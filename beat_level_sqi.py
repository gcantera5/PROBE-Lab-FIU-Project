import os
import glob
import json
import gzip
import shutil
from dataclasses import dataclass
from typing import Optional
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
from scipy.signal import cheby2, sosfiltfilt, find_peaks, resample, welch
from scipy.stats import skew

try:
    from dtw import dtw
except ImportError:
    dtw = None


# ============================================================
# CONFIG
# ============================================================

FS = 25  # FIU sampling frequency from our old code

CHANNEL_MAP = {
    "Unpolarized_A": {"Green": "c5",  "Red": "c2",  "IR": "c4"},
    "Unpolarized_B": {"Green": "c11", "Red": "c8",  "IR": "c10"},
    "Co-Polarized":  {"Green": "c13", "Red": "c12", "IR": "c15"},
    "Cross-Polarized": {"Green": "c19", "Red": "c18", "IR": "c21"},
}

# Folder names aren't fully consistent across days: some Day 3/4 folders say
# "Light" for the fair phantom, and some polarization labels end in a stray
# period ("Og. Pol." vs "Og. Pol"). Normalize so conditions group correctly.
SKIN_TONE_ALIASES = {"Light": "Fair"}


def normalize_skin_tone(skin):
    skin = skin.strip().capitalize()
    return SKIN_TONE_ALIASES.get(skin, skin)


def normalize_polarization_label(pol):
    return pol.strip().rstrip(".").strip()


WAVELENGTH_NM = {
    "Green": 525,
    "Red": 660,
    "IR": 940,
}


@dataclass
class BeatSQIConfig:
    """
    Config for the approved pipeline (per 10-second window):
      filter -> segment into beats (on the filtered signal)
             -> PI for each beat (from the RAW beat, same beat boundaries)
             -> normalize each beat (0-1)
             -> compare each beat to a template
             -> classify each beat -> roll up to a window label
                (bad if > bad_window_fraction_threshold of beats are bad)
             -> window PI = median of that window's beat PIs
    """
    fs: int = FS

    # Chebyshev Type II bandpass filter (used for beat segmentation/comparison)
    lowcut: float = 0.5
    highcut: float = 4.0
    filter_order: int = 4
    rs: int = 20

    # Window definition (10-second windows, 1-second stride)
    window_seconds: int = 10
    step_seconds: int = 1

    # Beat segmentation
    min_hr: float = 40
    max_hr: float = 180
    valley_prominence: Optional[float] = None
    # Valleys must be at least this fraction of the window's own pulse
    # period apart (period taken from the window's dominant frequency).
    # Stops the secondary notch in slower pulses from being counted as a
    # separate beat (Sept 2026: Slow recordings were being split in two).
    min_valley_spacing_fraction: float = 0.6

    # Beat normalization + template
    target_beat_length: int = 50
    n_template_beats: int = 12
    template_corr_threshold: float = 0.80

    # Beat classification thresholds.
    # sqi_lambda: Zia et al. (2020) value, used with their normalization
    # SQI = exp(-lambda * D / L), L = warping path length (see
    # compute_dtw_distance). This was also the pipeline's original value;
    # it had been changed to 0.30 on Sept 1 when the metrics briefly ran on
    # real (non-normalized) amplitudes, and was never set back afterwards.
    # (Sept 2026, after Rutendo's review.)
    #
    # MAD is no longer used as a rejection rule: it isn't an SQI in either
    # reference paper, and the old max_mad = 30 could never be reached on
    # 0-1 normalized beats. MAD is still saved as a descriptive column.
    sqi_lambda: float = 25.0
    min_template_sqi: float = 0.05
    min_corr: float = 0.80
    min_clipping_sqi: float = 0.80

    # Window rejection: bad if more than this fraction of a window's beats
    # are bad (starting value carried over from the original pipeline --
    # treat as tunable once this runs on real data)
    bad_window_fraction_threshold: float = 0.30


# ============================================================
# FILE HELPERS
# ============================================================

def unzip_file(filepath):
    """Unzip .json.gz files if needed."""
    if filepath.endswith(".gz"):
        new_path = filepath[:-3]
        if not os.path.exists(new_path):
            with gzip.open(filepath, "rb") as f_in, open(new_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        return new_path
    return filepath


def channel_name_map(clean_col_name):
    """Map cleaned channel names to hardware channel labels."""
    for pol, mapping in CHANNEL_MAP.items():
        for color, ch in mapping.items():
            if f"{pol}_{color}" == clean_col_name:
                return ch.upper()
    return "N/A"


def load_fiu_json(json_path, condition_info):
    """
    Load FIU JSON and return a cleaned dataframe with time + PPG channels.
    """
    json_path = unzip_file(json_path)

    with open(json_path, "r") as f:
        data = json.load(f)

    if "nirs4v1_adc24_32" in data:
        data_section = data["nirs4v1_adc24_32"]
    elif "nirs4v1_adc2" in data:
        data_section = data["nirs4v1_adc2"]
    elif "semi" in data:
        data_section = data["semi"]
    else:
        raise KeyError("Could not find expected FIU data key in JSON.")

    df = pd.DataFrame(data_section)
    time_col = df.get("ts", df.index)

    cleaned_data = {"time": time_col}

    for pol, mapping in CHANNEL_MAP.items():
        for color, channel in mapping.items():
            if channel in df.columns:
                cleaned_data[f"{pol}_{color}"] = df[channel]

    cleaned_df = pd.DataFrame(cleaned_data)

    for key, val in condition_info.items():
        cleaned_df[key] = val

    return cleaned_df


# ============================================================
# SIGNAL PREPROCESSING
# ============================================================

def bandpass_filter(ppg, config):
    """
    Chebyshev Type II bandpass filter
    """
    ppg = np.asarray(ppg, dtype=float)
    nyq = config.fs / 2.0

    low = config.lowcut / nyq
    high = config.highcut / nyq

    sos = cheby2(
        config.filter_order,
        config.rs,
        [low, high],
        btype="bandpass",
        output="sos"
    )

    return sosfiltfilt(sos, ppg)


def preprocess_ppg(ppg, config):
    """
    Preprocessing for one window: Chebyshev Type II bandpass filter only.
    Used to get a clean signal for beat segmentation and beat comparison.

    No separate median-filter baseline removal or zero-centering step
    (dropped per earlier discussion).
    """
    ppg = np.asarray(ppg, dtype=float)
    filtered = bandpass_filter(ppg, config)
    return np.nan_to_num(filtered, nan=0.0)


# ============================================================
# PERFUSION INDEX
# ============================================================

def compute_perfusion_index(signal):
    """
    Compute the perfusion index (PI) of ONE BEAT of the RAW (unfiltered)
    signal. Beat boundaries come from valley detection on the filtered
    window; the same indices are used to cut the raw window.

    Why per beat (Rutendo's review, Sept 2026): taking max/min over a whole
    10-second window folds baseline drift and artifacts into "AC". Doing it
    beat by beat keeps AC to one pulse, and the window PI is then the
    median across beats (robust to a few odd beats).

    Why raw, not filtered: the bandpass filter removes the DC baseline, so
    dividing by a near-zero "DC" sends PI into the hundreds/thousands of
    percent.

    DC = trough / minimum of the beat
    AC = peak - DC
    PI = (AC / DC) * 100
    """
    signal = np.asarray(signal, dtype=float)

    if len(signal) == 0:
        return np.nan, np.nan, np.nan

    dc_value = float(np.min(signal))
    peak_value = float(np.max(signal))
    ac_value = peak_value - dc_value

    if dc_value == 0:
        pi = np.nan
    else:
        pi = (ac_value / abs(dc_value)) * 100.0

    return ac_value, dc_value, pi


# ============================================================
# BEAT SEGMENTATION
# ============================================================

def estimate_pulse_period_samples(ppg, config):
    """
    Estimate the window's pulse period (in samples) from its dominant
    frequency inside the bandpass range. Returns None if it can't.
    """
    ppg = np.asarray(ppg, dtype=float)
    if len(ppg) < 2 * config.fs or np.std(ppg) == 0:
        return None

    freqs, power = welch(ppg - np.mean(ppg), fs=config.fs, nperseg=len(ppg))
    band = (freqs >= config.lowcut) & (freqs <= config.highcut)
    if not np.any(band):
        return None

    dominant_freq = freqs[band][np.argmax(power[band])]
    if dominant_freq <= 0:
        return None

    return config.fs / dominant_freq


def segment_beats_by_valleys(ppg, config):
    """
    Extract non-overlapping beats using valley-to-valley intervals, on the
    filtered window signal.

    Minimum valley spacing = the larger of
      - one beat at max_hr, and
      - min_valley_spacing_fraction x the window's own pulse period.
    The second rule keeps a pulse's secondary notch from being counted as
    its own valley (which was splitting every Slow-speed beat in two).
    """
    min_distance = int(config.fs * 60.0 / config.max_hr)

    period = estimate_pulse_period_samples(ppg, config)
    if period is not None:
        min_distance = max(
            min_distance,
            int(config.min_valley_spacing_fraction * period),
        )

    valleys, _ = find_peaks(
        -ppg,
        distance=min_distance,
        prominence=config.valley_prominence
    )

    raw_beats = []
    beat_indices = []

    for i in range(len(valleys) - 1):
        start = int(valleys[i])
        end = int(valleys[i + 1])

        duration = (end - start) / config.fs
        if duration <= 0:
            continue

        hr = 60.0 / duration

        if config.min_hr <= hr <= config.max_hr:
            raw_beats.append(ppg[start:end])
            beat_indices.append((start, end))

    return raw_beats, beat_indices


# ============================================================
# BEAT NORMALIZATION + TEMPLATE
# ============================================================

def safe_corr(x, y):
    """Safely compute correlation between two beats."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if np.std(x) == 0 or np.std(y) == 0:
        return 0.0

    return float(np.corrcoef(x, y)[0, 1])


def normalize_single_beat(beat, target_length):
    """
    Resample a beat to a uniform length, then min-max scale it to [0, 1].
    This is the "normalize the beat" step in the approved pipeline --
    beats are compared to the template AFTER this normalization.
    """
    beat = np.asarray(beat, dtype=float)
    resampled = resample(beat, target_length)

    beat_min = np.min(resampled)
    beat_range = np.max(resampled) - beat_min

    if beat_range > 0:
        resampled = (resampled - beat_min) / beat_range
    else:
        resampled = resampled - beat_min

    return resampled


def normalize_beats(raw_beats, config):
    """Normalize (resample + 0-1 scale) all beats in a window."""
    return np.array([
        normalize_single_beat(beat, config.target_beat_length)
        for beat in raw_beats
    ])


def make_clean_template(normalized_beats, config):
    """
    Build a clean template beat by averaging around 12 template-like
    beats, using the NORMALIZED beats (per the approved pipeline order:
    normalize each beat, then compare to the template).
    """
    if len(normalized_beats) == 0:
        raise ValueError("No beats available to build template.")

    n_initial = min(config.n_template_beats, len(normalized_beats))
    rough_template = np.mean(normalized_beats[:n_initial], axis=0)

    correlations = np.array([
        safe_corr(beat, rough_template)
        for beat in normalized_beats
    ])

    good_idx = np.where(correlations >= config.template_corr_threshold)[0]

    if len(good_idx) < min(3, len(normalized_beats)):
        good_idx = np.argsort(correlations)[::-1]

    selected_idx = good_idx[:min(config.n_template_beats, len(good_idx))]
    template = np.mean(normalized_beats[selected_idx], axis=0)

    return template, selected_idx


# ============================================================
# BEAT-LEVEL FEATURES
# ============================================================

def compute_dtw_distance(beat, template):
    """
    Compare beat to template with DTW, normalized the way Zia et al. (2020)
    do it: D(s,t) / L(s,t), where D is the summed distance along the warping
    path and L is the length of that path.

    Uses the "symmetric1" step pattern so every step on the path counts
    once (dtw-python's default, "symmetric2", double-weights diagonal steps
    and its normalizedDistance divides by N+M instead of the path length).
    """
    if dtw is None:
        raise ImportError("Please install dtw-python using: pip install dtw-python")

    alignment = dtw(beat, template, step_pattern="symmetric1", keep_internals=False)
    path_length = len(alignment.index1)
    return float(alignment.distance / path_length)


def compute_template_sqi(dtw_distance, config):
    """
    Zia et al. (2020) distance-based SQI: exp(-lambda * D / L), where
    dtw_distance is already D / L. Small distance = high SQI.
    """
    return float(np.exp(-config.sqi_lambda * dtw_distance))


def clipping_sqi(signal):
    """
    Estimate how much of a signal is clipped/saturated.
    1.0 = not clipped. Lower values = more clipping.
    Scale-invariant (a ratio relative to the signal's own min/max), so it
    gives the same result whether computed on the raw or normalized beat.
    """
    signal = np.asarray(signal, dtype=float)

    sig_range = np.max(signal) - np.min(signal)

    if sig_range == 0:
        return 0.0

    eps = 0.01 * sig_range

    near_max = signal >= (np.max(signal) - eps)
    near_min = signal <= (np.min(signal) + eps)

    clipped_fraction = np.mean(near_max | near_min)
    return float(1.0 - clipped_fraction)


def compute_beat_features(raw_beat, unfiltered_beat, normalized_beat, template, start_idx, end_idx, beat_number, config):
    """
    Compute beat-level features by comparing the NORMALIZED beat to the
    template (also built from normalized beats):
      - dtw_distance / template_sqi: Zia et al. (2020) distance-based SQI
      - correlation: Li & Clifford (2012) "linear resampling" SQI (beat
        resampled to a fixed length, then correlated with the template)
      - clipping_sqi: our own clipping check, loosely based on Li &
        Clifford's clipping detection
      - MAD, skewness: descriptive only (not used to reject beats)
    Plus the beat's perfusion index from the RAW (unfiltered) beat.
    (AD was dropped: with every beat resampled to target_beat_length
    points, AD = target_beat_length * MAD exactly, so it added nothing.)
    """
    difference = normalized_beat - template

    dtw_distance = compute_dtw_distance(normalized_beat, template)
    template_sqi = compute_template_sqi(dtw_distance, config)

    corr = max(0.0, safe_corr(normalized_beat, template))
    mad = float(np.mean(np.abs(difference)))
    beat_skewness = float(skew(normalized_beat, nan_policy="omit"))
    clip_sqi = clipping_sqi(raw_beat)
    beat_ac, beat_dc, beat_pi = compute_perfusion_index(unfiltered_beat)

    duration_sec = (end_idx - start_idx) / config.fs
    estimated_hr = 60.0 / duration_sec if duration_sec > 0 else np.nan

    return {
        "beat_number": beat_number,
        "beat_start_idx": start_idx,
        "beat_end_idx": end_idx,
        "beat_start_sec": start_idx / config.fs,
        "beat_end_sec": end_idx / config.fs,
        "duration_sec": duration_sec,
        "estimated_hr": estimated_hr,

        "dtw_distance": dtw_distance,
        "template_sqi": template_sqi,
        "correlation": corr,
        "MAD": mad,
        "skewness": beat_skewness,
        "clipping_sqi": clip_sqi,

        "beat_AC": beat_ac,
        "beat_DC": beat_dc,
        "beat_pi": beat_pi,
    }


def classify_beat(row, config):
    """
    Rule-based beat rejection (restored from the original pipeline).
    Also records why a beat was rejected.
    """
    reasons = []

    if row["template_sqi"] < config.min_template_sqi:
        reasons.append("low_template_sqi")

    if row["correlation"] < config.min_corr:
        reasons.append("low_correlation")

    if row["clipping_sqi"] < config.min_clipping_sqi:
        reasons.append("clipping")

    label = "bad" if reasons else "good"

    return label, ",".join(reasons)


# ============================================================
# BEAT COMPARISON -> WINDOW ROLL-UP
# ============================================================

def run_beat_comparison(filtered_window, raw_window, config):
    """
    Segment one filtered window into beats, normalize each beat, build a
    template from the normalized beats, compare every beat to that
    template, classify each beat, then roll the beat labels up into a
    window-level summary (bad if > bad_window_fraction_threshold of the
    beats are bad).

    The same beat boundaries are used to cut the RAW window, and PI is
    computed per raw beat. The window's PI is the median of its beat PIs
    (mean and SD of the beat PIs are saved too).

    Returns None if there aren't at least 3 beats in the window (matches
    the original pipeline's minimum for building a usable template).
    """
    raw_beats, beat_indices = segment_beats_by_valleys(filtered_window, config)

    if len(raw_beats) < 3:
        return None

    normalized_beats = normalize_beats(raw_beats, config)
    template, template_indices = make_clean_template(normalized_beats, config)

    rows = []

    for beat_number, (raw_beat, normalized_beat, (start_idx, end_idx)) in enumerate(
        zip(raw_beats, normalized_beats, beat_indices)
    ):
        row = compute_beat_features(
            raw_beat=raw_beat,
            unfiltered_beat=raw_window[start_idx:end_idx],
            normalized_beat=normalized_beat,
            template=template,
            start_idx=start_idx,
            end_idx=end_idx,
            beat_number=beat_number,
            config=config,
        )

        label, reasons = classify_beat(row, config)
        row["beat_label"] = label
        row["rejection_reasons"] = reasons
        row["used_for_template"] = beat_number in set(template_indices.tolist())

        rows.append(row)

    feature_table = pd.DataFrame(rows)

    percent_bad = float(np.mean(feature_table["beat_label"] == "bad"))

    window_label = (
        "bad_window"
        if percent_bad > config.bad_window_fraction_threshold
        else "good_window"
    )

    beat_pis = feature_table["beat_pi"].replace([np.inf, -np.inf], np.nan)

    summary = {
        "AC": float(feature_table["beat_AC"].median()),
        "DC": float(feature_table["beat_DC"].median()),
        "perfusion_index": float(beat_pis.median()),
        "mean_beat_pi": float(beat_pis.mean()),
        "std_beat_pi": float(beat_pis.std()),
        "num_beats": len(feature_table),
        "num_good_beats": int(np.sum(feature_table["beat_label"] == "good")),
        "num_bad_beats": int(np.sum(feature_table["beat_label"] == "bad")),
        "percent_bad": percent_bad,
        "mean_template_sqi": float(feature_table["template_sqi"].mean()),
        "mean_dtw_distance": float(feature_table["dtw_distance"].mean()),
        "mean_correlation": float(feature_table["correlation"].mean()),
        "mean_MAD": float(feature_table["MAD"].mean()),
        "mean_clipping_sqi": float(feature_table["clipping_sqi"].mean()),
        "mean_skewness": float(feature_table["skewness"].mean()),
        "window_label": window_label,
    }

    return {
        "feature_table": feature_table,
        "template": template,
        "normalized_beats": normalized_beats,
        "summary": summary,
    }


# ============================================================
# MAIN WINDOW-LEVEL PIPELINE
# ============================================================

def run_window_level_sqi(raw_window, config=None):
    """
    Run the approved pipeline on one 10-second PPG window:

        filter (for beat segmentation/comparison)
        -> segment into beats -> PI per RAW beat (window PI = median)
        -> normalize each beat
        -> compare each beat to a template -> classify each beat
        -> roll up to a window label (good_window / bad_window)

    Returns None if the window doesn't have enough beats to build a
    template (same >=3-beat minimum as before).
    """
    if config is None:
        config = BeatSQIConfig()

    raw_window = np.asarray(raw_window, dtype=float)
    filtered = preprocess_ppg(raw_window, config)

    beat_result = run_beat_comparison(filtered, raw_window, config)

    if beat_result is None:
        return None

    row = dict(beat_result["summary"])

    return {
        "filtered_signal": filtered,
        "feature_table": beat_result["feature_table"],
        "template": beat_result["template"],
        "normalized_beats": beat_result["normalized_beats"],
        "row": row,
    }


# ============================================================
# WINDOW-LEVEL WRAPPER
# ============================================================

def iter_windows(signal, config):
    """
    Sliding windows: 10-second windows with 1-second stride.
    """
    signal = np.asarray(signal, dtype=float)

    win_len = int(config.window_seconds * config.fs)
    step_len = int(config.step_seconds * config.fs)

    if len(signal) < win_len:
        return

    for start in range(0, len(signal) - win_len + 1, step_len):
        end = start + win_len
        yield start, end, signal[start:end]


def run_sqi_over_windows(signal, condition_info, channel_label, config=None):
    """
    Run the approved pipeline over 10-second FIU-style windows.

    Returns one DataFrame with one row per window (PI from the raw
    window, beat-comparison SQI values rolled up to the window, good/bad
    label, window/source info). Windows with fewer than 3 beats are
    skipped (same as the original pipeline).
    """
    if config is None:
        config = BeatSQIConfig()

    window_rows = []

    for win_start, win_end, win in iter_windows(signal, config):
        result = run_window_level_sqi(win, config)

        if result is None:
            continue

        row = result["row"]

        window_row = {
            "Channel": channel_label,
            "Hardware Channel": channel_name_map(channel_label),
            "WindowStartIdx": win_start,
            "WindowEndIdx": win_end,
            "WindowStartSec": win_start / config.fs,
            "WindowEndSec": win_end / config.fs,
            **row,
            **condition_info,
        }

        window_rows.append(window_row)

    return pd.DataFrame(window_rows)

def parse_day4_folder(folder_name):
    """
    Day 4 folder examples:
        Green Dark 0 Og. Pol
        Red Fair 0 Flip. Pol
        IR Dark 90 Un. Pol
        IR Fair 180 Og. Pol
    """
    parts = folder_name.split()

    if len(parts) < 4:
        raise ValueError(f"Could not parse Day 4 folder name: {folder_name}")

    wavelength = parts[0]
    skin = normalize_skin_tone(parts[1])
    orientation = parts[2]
    pol = normalize_polarization_label(" ".join(parts[3:]))

    return wavelength, skin, orientation, pol

def parse_day3_folder(folder_name):
    """
    Day 3 heartbeat folder examples:
        Green Dark
        Green Light
        Red Dark
        Red Light
    """
    parts = folder_name.split()

    if len(parts) != 2:
        raise ValueError(
            f"Could not parse Day 3 folder name: {folder_name}"
        )

    wavelength = parts[0].capitalize()
    skin = normalize_skin_tone(parts[1])

    return wavelength, skin

# ============================================================
# DAY 2: EXPERIMENT 1 PROCESSING
# ============================================================

def process_experiment1_complete(
    experiment_root="Experiment 1 Complete  copy",
    output_root="FIU_Beat_Level_SQI/Day_2/Experiment_1"
):
    """
    Process Experiment 1 folder structure:

    Experiment 1 Complete copy
        2.5 Dark Fast
        2.5 Dark Intermediate
        2.5 Dark Slow
        ...
    """
    config = BeatSQIConfig()
    os.makedirs(output_root, exist_ok=True)

    all_window_results = []
    all_summary_results = []


    condition_folders = [
        d for d in glob.glob(os.path.join(experiment_root, "*"))
        if os.path.isdir(d)
    ]
    

    for folder in condition_folders:
        folder_name = os.path.basename(folder)
        parts = folder_name.split()

        if len(parts) < 3:
            print(f"Skipping folder with unexpected name: {folder_name}")
            continue

        depth = parts[0] + "mm"
        skin = parts[1].capitalize()
        speed = parts[2].capitalize()

        condition_info = {
            "Day": "Day_2",
            "Experiment": "Experiment_1",
            "SkinTone": skin,
            "Speed": speed,
            "Depth": depth,
            "ExpectedBPM": {
                "Slow": 60,
                "Intermediate": 90,
                "Fast": 120
            }.get(speed),
            "Clamp": "Yes",
            "PolarizationPlacement": "Same",
            "ConditionFolder": folder_name,
        }

        # --------------------------------------------------
        # STEP 1: unzip any .json.gz files
        # --------------------------------------------------

        gz_files = glob.glob(os.path.join(folder, "*.json.gz"))

        for gz_file in gz_files:
            unzip_file(gz_file)

        # --------------------------------------------------
        # STEP 2: ONLY process .json files
        # --------------------------------------------------

        json_files = sorted(glob.glob(os.path.join(folder, "*.json")))

        for trial_number, json_path in enumerate( json_files, start=1):


            print(f"\nProcessing: {folder_name} / {os.path.basename(json_path)}")

            try:
                cleaned_df = load_fiu_json(json_path, condition_info)
            except Exception as e:
                print(f"Could not load {json_path}: {e}")
                continue

            for col in cleaned_df.columns:
                if not any(pol in col for pol in CHANNEL_MAP.keys()):
                    continue

                signal = cleaned_df[col].values

                try:
                    window_df = run_sqi_over_windows(
                        signal=signal,
                        condition_info=condition_info,
                        channel_label=col,
                        config=config
                    )
                    if len(window_df) > 0:
                        print("\n--------------------------------")
                        print(f"Condition: {folder_name}")
                        print(f"Channel: {col}")
                        print(f"Mean PI: {window_df['perfusion_index'].mean():.4f}%")
                        print(f"Mean DTW: {window_df['mean_dtw_distance'].mean():.4f}")
                        print(f"Mean Correlation: {window_df['mean_correlation'].mean():.4f}")
                        print(f"Mean MAD: {window_df['mean_MAD'].mean():.4f}")

                        good_windows = (window_df["window_label"] == "good_window").sum()
                        total_windows = len(window_df)
                        bad_windows = total_windows - good_windows

                        good_window_fraction = (
                            good_windows / total_windows
                            if total_windows > 0
                            else np.nan
                        )

                        print(f"Good Windows: {good_windows}/{total_windows}")
                        print("--------------------------------")

                        summary_row = {
                            "Day": "Day_2",
                            "Experiment": "Experiment_1",
                            "ConditionFolder": folder_name,
                            "Trial": trial_number,
                            "SourceFile": os.path.basename(json_path),
                            "Channel": col,
                            "Hardware Channel": channel_name_map(col),
                            "SkinTone": skin,
                            "Speed": speed,
                            "Depth": depth,
                            "ExpectedBPM": {
                                "Slow": 60,
                                "Intermediate": 90,
                                "Fast": 120
                            }.get(speed),
                            "Clamp": "Yes",
                            "PolarizationPlacement": "Same",
                            "GoodWindows": int(good_windows),
                            "BadWindows": int(bad_windows),
                            "TotalWindows": int(total_windows),
                            "GoodWindowFraction": float(
                                good_window_fraction
                            ),
                            "GoodWindowPercent": float(
                                good_window_fraction * 100.0
                            ),
                            "BadWindowPercent": float(
                                (1.0 - good_window_fraction) * 100.0
                            ),
                            "MeanPI": float(
                                window_df["perfusion_index"].mean()
                            ),
                            "MeanDTW": float(
                                window_df["mean_dtw_distance"].mean()
                            ),
                            "MeanCorrelation": float(
                                window_df["mean_correlation"].mean()
                            ),
                            "MeanMAD": float(
                                window_df["mean_MAD"].mean()
                            ),
                            "MeanTemplateSQI": float(
                                window_df["mean_template_sqi"].mean()
                            ),
                            "MeanClippingSQI": float(
                                window_df["mean_clipping_sqi"].mean()
                            ),
                            "MeanSkewness": float(
                                window_df["mean_skewness"].mean()
                            ),
                        }

                        all_summary_results.append(summary_row)

                except Exception as e:
                    print(f"SQI failed for {col}: {e}")
                    continue

                if len(window_df) > 0:
                    window_df["SourceFile"] = os.path.basename(json_path)
                    window_df["Trial"] = trial_number
                    all_window_results.append(window_df)

    if all_window_results:
        final_window_df = pd.concat(all_window_results, ignore_index=True)
        window_path = os.path.join(output_root, "day2_experiment1_all_window_sqi.csv")
        final_window_df.to_csv(window_path, index=False)
        print(f"\nSaved window-level SQI results to: {window_path}")

    if all_summary_results:
        summary_df = pd.DataFrame(all_summary_results)

        summary_path = os.path.join(
            output_root,
            "day2_experiment1_recording_summary.csv"
        )

        summary_df.to_csv(
            summary_path,
            index=False
        )

        print(
            f"Saved Day 2 Experiment 1 summary to: "
            f"{summary_path}"
        )

# ============================================================
# DAY 3: EXPERIMENT 2 HEARTBEAT PROCESSING
# ============================================================

def process_day3_experiment2(
    day3_root=(
        "Experiment 2 Test (Day 3) copy/"
        "Multilayered, 90 BPM, No Clamps & OG Polarization"
    ),
    output_root="FIU_Beat_Level_SQI/Day_3/Experiment_2"
):
    """
    Process the main Day 3 Experiment 2 pulsatile recordings.

    Folder structure:
        Multilayered, 90 BPM, No Clamps & OG Polarization
            Green Dark
            Green Light
            Red Dark
            Red Light

    Fixed experimental conditions:
        Depth = 3.5 mm
        Speed = Intermediate
        Expected BPM = 90
        Clamp = No
        Phantom = Multilayered
        Polarization = Original
    """
    config = BeatSQIConfig()
    os.makedirs(output_root, exist_ok=True)

    all_window_results = []
    all_summary_results = []

    if not os.path.exists(day3_root):
        print(f"Day 3 heartbeat folder not found: {day3_root}")
        return

    condition_folders = sorted(
        folder
        for folder in glob.glob(os.path.join(day3_root, "*"))
        if os.path.isdir(folder)
    )

    print(
        f"\nFound {len(condition_folders)} "
        f"Day 3 heartbeat condition folders."
    )

    for folder in condition_folders:
        folder_name = os.path.basename(folder)

        try:
            wavelength, skin = parse_day3_folder(folder_name)
        except ValueError as error:
            print(error)
            continue

        condition_info = {
            "Day": "Day_3",
            "Experiment": "Experiment_2",
            "ConditionFolder": folder_name,
            "Wavelength": wavelength,
            "WavelengthNm": WAVELENGTH_NM.get(wavelength),
            "SkinTone": skin,
            "Depth": "3.5mm",
            "Speed": "Intermediate",
            "ExpectedBPM": 90,
            "Clamp": "No",
            "PhantomType": "Multilayered",
            "PolarizationCondition": "Original",
        }

        # --------------------------------------------------
        # STEP 1: unzip compressed files if needed
        # --------------------------------------------------
        gz_files = sorted(
            glob.glob(os.path.join(folder, "*.json.gz"))
        )

        for gz_file in gz_files:
            unzip_file(gz_file)

        # --------------------------------------------------
        # STEP 2: process ONLY unzipped JSON files
        # --------------------------------------------------
        json_files = sorted(
            glob.glob(os.path.join(folder, "*.json"))
        )

        print(
            f"\n{folder_name}: found "
            f"{len(json_files)} unzipped JSON files."
        )

        for trial_number, json_path in enumerate(
            json_files,
            start=1
        ):
            source_file = os.path.basename(json_path)

            print(
                f"\nProcessing Day 3 / Experiment 2 / "
                f"{folder_name} / Trial {trial_number}"
            )

            try:
                cleaned_df = load_fiu_json(
                    json_path,
                    condition_info
                )
            except Exception as error:
                print(f"Could not load {json_path}: {error}")
                continue

            for col in cleaned_df.columns:
                if not any(
                    polarization in col
                    for polarization in CHANNEL_MAP.keys()
                ):
                    continue

                signal = cleaned_df[col].values

                try:
                    window_df = run_sqi_over_windows(
                        signal=signal,
                        condition_info=condition_info,
                        channel_label=col,
                        config=config
                    )
                except Exception as error:
                    print(f"SQI failed for {col}: {error}")
                    continue

                if len(window_df) == 0:
                    print(
                        f"No valid windows found for "
                        f"{folder_name} / {col}"
                    )
                    continue

                good_windows = (
                    window_df["window_label"] == "good_window"
                ).sum()

                total_windows = len(window_df)
                bad_windows = total_windows - good_windows

                good_window_fraction = (
                    good_windows / total_windows
                    if total_windows > 0
                    else np.nan
                )

                print("\n--------------------------------")
                print("Day: Day 3")
                print("Experiment: Experiment 2")
                print(f"Condition: {folder_name}")
                print(f"Trial: {trial_number}")
                print(f"Channel: {col}")
                print(
                    "Mean PI: "
                    f"{window_df['perfusion_index'].mean():.4f}%"
                )
                print(
                    "Mean DTW: "
                    f"{window_df['mean_dtw_distance'].mean():.4f}"
                )
                print(
                    "Mean Correlation: "
                    f"{window_df['mean_correlation'].mean():.4f}"
                )
                print(
                    "Mean MAD: "
                    f"{window_df['mean_MAD'].mean():.4f}"
                )
                print(
                    f"Good Windows: "
                    f"{good_windows}/{total_windows}"
                )
                print("--------------------------------")

                summary_row = {
                    "Day": "Day_3",
                    "Experiment": "Experiment_2",
                    "ConditionFolder": folder_name,
                    "Trial": trial_number,
                    "SourceFile": source_file,
                    "Channel": col,
                    "Hardware Channel": channel_name_map(col),
                    "Wavelength": wavelength,
                    "WavelengthNm": WAVELENGTH_NM.get(wavelength),
                    "SkinTone": skin,
                    "Depth": "3.5mm",
                    "Speed": "Intermediate",
                    "ExpectedBPM": 90,
                    "Clamp": "No",
                    "PhantomType": "Multilayered",
                    "PolarizationCondition": "Original",
                    "GoodWindows": int(good_windows),
                    "BadWindows": int(bad_windows),
                    "TotalWindows": int(total_windows),
                    "GoodWindowFraction": float(
                        good_window_fraction
                    ),
                    "GoodWindowPercent": float(
                        good_window_fraction * 100.0
                    ),
                    "BadWindowPercent": float(
                        (1.0 - good_window_fraction) * 100.0
                    ),
                    "MeanPI": float(
                        window_df["perfusion_index"].mean()
                    ),
                    "MeanDTW": float(
                        window_df["mean_dtw_distance"].mean()
                    ),
                    "MeanCorrelation": float(
                        window_df["mean_correlation"].mean()
                    ),
                    "MeanMAD": float(
                        window_df["mean_MAD"].mean()
                    ),
                    "MeanTemplateSQI": float(
                        window_df["mean_template_sqi"].mean()
                    ),
                    "MeanClippingSQI": float(
                        window_df["mean_clipping_sqi"].mean()
                    ),
                    "MeanSkewness": float(
                        window_df["mean_skewness"].mean()
                    ),
                }

                all_summary_results.append(summary_row)

                window_df["SourceFile"] = source_file
                window_df["Trial"] = trial_number
                all_window_results.append(window_df)

    # --------------------------------------------------
    # SAVE WINDOW RESULTS
    # --------------------------------------------------
    if all_window_results:
        final_window_df = pd.concat(
            all_window_results,
            ignore_index=True
        )

        window_path = os.path.join(
            output_root,
            "day3_experiment2_all_window_sqi.csv"
        )

        final_window_df.to_csv(
            window_path,
            index=False
        )

        print(
            f"\nSaved Day 3 window results to: "
            f"{window_path}"
        )

    # --------------------------------------------------
    # SAVE RECORDING SUMMARY
    # --------------------------------------------------
    if all_summary_results:
        summary_df = pd.DataFrame(all_summary_results)

        summary_path = os.path.join(
            output_root,
            "day3_experiment2_recording_summary.csv"
        )

        summary_df.to_csv(
            summary_path,
            index=False
        )

        print(
            f"Saved Day 3 recording summary to: "
            f"{summary_path}"
        )



def merge_shared_recordings(df, key_cols):
    """
    Some Day 4 recordings (IR 0-degree Og. Pol) are used in BOTH Experiment 2
    and Experiment 3. Keep them in each experiment's own file, but only once
    in the combined Day 4 file, labeled with both experiments.
    """
    experiments_per_file = df.groupby("SourceFile")["Experiment"].agg(
        lambda s: " & ".join(sorted(set(s)))
    )
    merged = df.copy()
    merged["Experiment"] = merged["SourceFile"].map(experiments_per_file)
    return merged.drop_duplicates(subset=key_cols).reset_index(drop=True)


def process_day4_experiments(
    day4_root="Experiment 2 & 3 (Day 4) copy",
    output_root="FIU_Beat_Level_SQI/Day_4"
):
    """
    Process Day 4 heartbeat/SQI data following the experiment handout.

    Day 4:
        Experiment 2 = wavelength / skin tone / polarization condition
        Experiment 3 = IR orientation / skin tone / polarization condition

    All Day 4 data:
        Depth = 3.5 mm
        Speed = Intermediate / 90 BPM
        Clamps = Yes
        Phantom = Multilayered
    """
    # -------------------------------
    # DEBUG: Verify folder structure
    # -------------------------------

    config = BeatSQIConfig()
    os.makedirs(output_root, exist_ok=True)

    all_window_results = []
    all_summary_results = []

    for exp_label in ["Experiment 2", "Experiment 3"]:

        exp_path = os.path.join(day4_root, exp_label)

        if not os.path.exists(exp_path):
            print(f"Skipping missing folder: {exp_path}")
            continue

        condition_folders = [
            d for d in glob.glob(os.path.join(exp_path, "*"))
            if os.path.isdir(d)
        ]

        for folder in condition_folders:

            folder_name = os.path.basename(folder)

            try:
                wavelength, skin, orientation, pol = parse_day4_folder(folder_name)
            except Exception as e:
                print(e)
                continue

            condition_info = {
                "Day": "Day_4",
                "Experiment": exp_label.replace(" ", "_"),
                "Wavelength": wavelength,
                "WavelengthNm": WAVELENGTH_NM.get(wavelength),
                "SkinTone": skin,
                "OrientationDegrees": int(orientation),
                "PolarizationCondition": pol,
                "Depth": "3.5mm",
                "Speed": "Intermediate",
                "ExpectedBPM": 90,
                "Clamp": "Yes",
                "PhantomType": "Multilayered",
                "ConditionFolder": folder_name,
            }

            # --------------------------------------------------
            # STEP 1: unzip .json.gz files if needed
            # --------------------------------------------------
            gz_files = glob.glob(os.path.join(folder, "*.json.gz"))

            for gz_file in gz_files:
                unzip_file(gz_file)

            # --------------------------------------------------
            # STEP 2: process ONLY unzipped .json files
            # --------------------------------------------------
            json_files = sorted(glob.glob(os.path.join(folder, "*.json")))

            for trial_number, json_path in enumerate(json_files, start=1):

                print(
                    f"\nProcessing: {exp_label} / "
                    f"{folder_name} / {os.path.basename(json_path)}"
                )

                try:
                    cleaned_df = load_fiu_json(json_path, condition_info)
                except Exception as e:
                    print(f"Could not load {json_path}: {e}")
                    continue

                for col in cleaned_df.columns:

                    if not any(pol_key in col for pol_key in CHANNEL_MAP.keys()):
                        continue

                    signal = cleaned_df[col].values

                    try:
                        window_df = run_sqi_over_windows(
                            signal=signal,
                            condition_info=condition_info,
                            channel_label=col,
                            config=config
                        )

                        if len(window_df) > 0:
                            print("\n--------------------------------")
                            print(f"Experiment: {exp_label}")
                            print(f"Condition: {folder_name}")
                            print(f"Channel: {col}")
                            print(f"Mean PI: " f"{window_df['perfusion_index'].mean():.4f}%")
                            print(f"Mean DTW: {window_df['mean_dtw_distance'].mean():.4f}")
                            print(f"Mean Correlation: {window_df['mean_correlation'].mean():.4f}")
                            print(f"Mean MAD: {window_df['mean_MAD'].mean():.4f}")

                            good_windows = (
                                window_df["window_label"] == "good_window"
                            ).sum()

                            total_windows = len(window_df)
                            bad_windows = total_windows - good_windows

                            good_window_fraction = (
                                good_windows / total_windows
                                if total_windows > 0
                                else np.nan
                            )

                            print(f"Good Windows: {good_windows}/{total_windows}")
                            print("--------------------------------")

                            summary_row = {
                                "Day": "Day_4",
                                "Experiment": exp_label.replace(" ", "_"),
                                "ConditionFolder": folder_name,
                                "Trial": trial_number,
                                "SourceFile": os.path.basename(json_path),
                                "Channel": col,
                                "Hardware Channel": channel_name_map(col),
                                "Wavelength": wavelength,
                                "WavelengthNm": WAVELENGTH_NM.get(wavelength),
                                "SkinTone": skin,
                                "OrientationDegrees": int(orientation),
                                "PolarizationCondition": pol,
                                "Depth": "3.5mm",
                                "Speed": "Intermediate",
                                "ExpectedBPM": 90,
                                "Clamp": "Yes",
                                "PhantomType": "Multilayered",
                                "GoodWindows": int(good_windows),
                                "BadWindows": int(bad_windows),
                                "TotalWindows": int(total_windows),
                                "GoodWindowFraction": float(
                                    good_window_fraction
                                ),
                                "GoodWindowPercent": float(
                                    good_window_fraction * 100.0
                                ),
                                "BadWindowPercent": float(
                                    (1.0 - good_window_fraction) * 100.0
                                ),
                                "MeanPI": float(
                                    window_df["perfusion_index"].mean()
                                ),
                                "MeanDTW": float(
                                    window_df["mean_dtw_distance"].mean()
                                ),
                                "MeanCorrelation": float(
                                    window_df["mean_correlation"].mean()
                                ),
                                "MeanMAD": float(
                                    window_df["mean_MAD"].mean()
                                ),
                                "MeanTemplateSQI": float(
                                    window_df["mean_template_sqi"].mean()
                                ),
                                "MeanClippingSQI": float(
                                    window_df["mean_clipping_sqi"].mean()
                                ),
                                "MeanSkewness": float(
                                    window_df["mean_skewness"].mean()
                                ),
                            }

                            all_summary_results.append(summary_row)


                    except Exception as e:
                        print(f"SQI failed for {col}: {e}")
                        continue

                    if len(window_df) > 0:
                        window_df["SourceFile"] = os.path.basename(json_path)
                        window_df["Trial"] = trial_number
                        all_window_results.append(window_df)


    # --------------------------------------------------
    # SAVE WINDOW-LEVEL RESULTS
    # --------------------------------------------------
    if all_window_results:
        final_window_df = pd.concat(
            all_window_results,
            ignore_index=True
        )

        # Save one combined Day 4 file
        window_path = os.path.join(
            output_root,
            "day4_all_window_sqi.csv"
        )

        merge_shared_recordings(
            final_window_df,
            key_cols=["SourceFile", "Channel", "WindowStartIdx"],
        ).to_csv(
            window_path,
            index=False
        )

        print(
            f"\nSaved Day 4 window-level SQI results to: "
            f"{window_path}"
        )

        # Save Experiment 2 and Experiment 3 separately
        for experiment_name in ["Experiment_2", "Experiment_3"]:

            experiment_window_df = final_window_df[
                final_window_df["Experiment"] == experiment_name
            ]

            if len(experiment_window_df) > 0:

                experiment_folder = os.path.join(
                    output_root,
                    experiment_name
                )

                os.makedirs(
                    experiment_folder,
                    exist_ok=True
                )

                experiment_window_path = os.path.join(
                    experiment_folder,
                    f"day4_{experiment_name.lower()}_window_sqi.csv"
                )

                experiment_window_df.to_csv(
                    experiment_window_path,
                    index=False
                )

                print(
                    f"Saved {experiment_name} window results to: "
                    f"{experiment_window_path}"
                )

    # --------------------------------------------------
    # SAVE RECORDING SUMMARY RESULTS
    # --------------------------------------------------
    if all_summary_results:
        summary_df = pd.DataFrame(all_summary_results)

        summary_path = os.path.join(
            output_root,
            "day4_recording_summary.csv"
        )

        merge_shared_recordings(
            summary_df,
            key_cols=["SourceFile", "Channel"],
        ).to_csv(
            summary_path,
            index=False
        )

        print(
            f"Saved Day 4 recording summary to: "
            f"{summary_path}"
        )

        # --------------------------------------------------
        # SAVE EXPERIMENT-SPECIFIC SUMMARIES
        # --------------------------------------------------

        for experiment_name in ["Experiment_2", "Experiment_3"]:

            experiment_summary_df = summary_df[
                summary_df["Experiment"] == experiment_name
            ]

            if len(experiment_summary_df) > 0:

                experiment_folder = os.path.join(
                    output_root,
                    experiment_name
                )

                os.makedirs(
                    experiment_folder,
                    exist_ok=True
                )

                summary_file = os.path.join(
                    experiment_folder,
                    f"{experiment_name.lower()}_summary.csv"
                )

                experiment_summary_df.to_csv(
                    summary_file,
                    index=False
                )

                print(
                    f"Saved {experiment_name} summary to: "
                    f"{summary_file}"
                )

# --------------------------------------------------
# Visual test for Day 2 Experiment 1
# --------------------------------------------------

def debug_one_file_one_channel(
    json_path,
    condition_info,
    channel_label="Cross-Polarized_IR",
):
    """
    Debug one recording and one channel so we can visually check:
    1. detected valleys/beats in the first window
    2. normalized beats + template for that window
    3. good/bad beat intervals within that window
    4. good/bad WINDOW labels over the whole recording
    """
    config = BeatSQIConfig()

    cleaned_df = load_fiu_json(json_path, condition_info)

    if channel_label not in cleaned_df.columns:
        print(f"{channel_label} not found.")
        print("Available columns:")
        print(cleaned_df.columns.tolist())
        return

    signal = cleaned_df[channel_label].values

    window_df = run_sqi_over_windows(
        signal=signal,
        condition_info=condition_info,
        channel_label=channel_label,
        config=config
    )

    if len(window_df) == 0:
        print("No result. Recording is shorter than one window, or no window had enough beats.")
        return

    print("\n===== DEBUG SUMMARY (window level) =====")
    print(f"Channel: {channel_label}")
    print(window_df[[
        "WindowStartSec",
        "WindowEndSec",
        "AC",
        "DC",
        "perfusion_index",
        "num_beats",
        "num_bad_beats",
        "mean_dtw_distance",
        "mean_correlation",
        "mean_MAD",
        "mean_template_sqi",
        "window_label",
    ]])

    good_windows = (window_df["window_label"] == "good_window").sum()
    total_windows = len(window_df)
    print(
        f"Good Windows: {good_windows}/{total_windows} "
        f"({100.0 * good_windows / total_windows:.1f}% good, "
        f"{100.0 * (total_windows - good_windows) / total_windows:.1f}% bad)"
    )

    # Detailed beat-level look at the first window that had enough beats
    win_len = int(config.window_seconds * config.fs)
    detail_result = None
    detail_start = None

    for win_start, win_end, win in iter_windows(signal, config):
        detail_result = run_window_level_sqi(win, config)
        if detail_result is not None:
            detail_start = win_start
            break

    if detail_result is None:
        print("No window had enough beats for a detailed beat-level plot.")
    else:
        filtered = detail_result["filtered_signal"]
        feature_table = detail_result["feature_table"]
        normalized_beats = detail_result["normalized_beats"]
        template = detail_result["template"]

        print(f"\n===== DEBUG SUMMARY (beat level, window starting at {detail_start / config.fs:.1f}s) =====")
        print(feature_table[[
            "beat_number",
            "dtw_distance",
            "correlation",
            "MAD",
            "beat_pi",
            "template_sqi",
            "clipping_sqi",
            "beat_label",
            "rejection_reasons",
        ]])

        # Plot 1: filtered window with detected valleys
        t_win = np.arange(len(filtered)) / config.fs
        valley_idxs = feature_table["beat_start_idx"].values.astype(int)

        plt.figure(figsize=(10, 4))
        plt.plot(t_win, filtered, linewidth=1.2)
        plt.scatter(
            valley_idxs / config.fs,
            filtered[valley_idxs],
            marker="v",
            s=50,
            label="Detected valleys"
        )
        plt.title(f"{channel_label}: Filtered Window with Detected Valleys")
        plt.xlabel("Time (s)")
        plt.ylabel("Filtered PPG")
        plt.legend()
        plt.grid(True, linestyle="--", alpha=0.4)
        plt.tight_layout()
        plt.show()

        # Plot 2: normalized beats + template
        plt.figure(figsize=(8, 4))
        for beat in normalized_beats:
            plt.plot(beat, alpha=0.25, linewidth=1)

        plt.plot(template, linewidth=3, label="Template beat")
        plt.title(f"{channel_label}: Normalized Beats + Template")
        plt.xlabel("Resampled sample")
        plt.ylabel("Normalized amplitude (0-1)")
        plt.legend()
        plt.grid(True, linestyle="--", alpha=0.4)
        plt.tight_layout()
        plt.show()

        # Plot 3: good/bad beat intervals within this window
        plt.figure(figsize=(10, 4))
        plt.plot(t_win, filtered, linewidth=1.2, color="black")

        for _, row in feature_table.iterrows():
            start = row["beat_start_idx"] / config.fs
            end = row["beat_end_idx"] / config.fs

            if row["beat_label"] == "good":
                plt.axvspan(start, end, color="green", alpha=0.18)
            else:
                plt.axvspan(start, end, color="red", alpha=0.30)

        plt.title(f"{channel_label}: Good vs Bad Beat Intervals")
        plt.xlabel("Time (s)")
        plt.ylabel("Filtered PPG")
        plt.grid(True, linestyle="--", alpha=0.4)

        plt.plot([], [], color="green", linewidth=8, alpha=0.4, label="Good beat")
        plt.plot([], [], color="red", linewidth=8, alpha=0.4, label="Bad beat")
        plt.legend()

        plt.tight_layout()
        plt.show()

    # Plot 4: good/bad WINDOW shading over the whole recording
    t_full = np.arange(len(signal)) / config.fs

    plt.figure(figsize=(10, 4))
    plt.plot(t_full, signal, linewidth=1.0, color="black", alpha=0.7)

    for _, row in window_df.iterrows():
        start = row["WindowStartSec"]
        end = row["WindowEndSec"]

        if row["window_label"] == "good_window":
            plt.axvspan(start, end, color="green", alpha=0.08)
        else:
            plt.axvspan(start, end, color="red", alpha=0.12)

    plt.title(f"{channel_label}: Good vs Bad Windows")
    plt.xlabel("Time (s)")
    plt.ylabel("Raw PPG")
    plt.grid(True, linestyle="--", alpha=0.4)

    plt.plot([], [], color="green", linewidth=8, alpha=0.3, label="Good window")
    plt.plot([], [], color="red", linewidth=8, alpha=0.3, label="Bad window")
    plt.legend()

    plt.tight_layout()
    plt.show()


# ============================================================
# MAIN BLOCK
# ============================================================

if __name__ == "__main__":

    # ------------------------------------------
    # Day 2: Experiment 1
    # ------------------------------------------
    process_experiment1_complete(
        experiment_root="Experiment 1 Complete  copy",
        output_root="FIU_Beat_Level_SQI/Day_2/Experiment_1"
    )

    # ------------------------------------------
    # Day 3: Experiment 2 heartbeat recordings
    # ------------------------------------------
    process_day3_experiment2(
        day3_root=(
            "Experiment 2 Test (Day 3) copy/"
            "Multilayered, 90 BPM, No Clamps & OG Polarization"
        ),
        output_root="FIU_Beat_Level_SQI/Day_3/Experiment_2"
    )

    # ------------------------------------------
    # Day 4: Experiments 2 and 3
    # ------------------------------------------
    process_day4_experiments(
        day4_root="Experiment 2 & 3 (Day 4) copy",
        output_root="FIU_Beat_Level_SQI/Day_4"
    )

# if __name__ == "__main__":
#
#     EXPERIMENT_1_FOLDER = "Experiment 1 Complete  copy"
#
#     process_experiment1_complete(
#         experiment_root=EXPERIMENT_1_FOLDER,
#         output_root="FIU_Beat_Level_SQI/Day_2/Experiment_1"
#     )

# if __name__ == "__main__":

#     test_json = "Experiment 1 Complete  copy/3.75 Fair Intermediate/2025-10-23T01-37-59-8516317d-a527-4f06-baaf-87ac9ffde0e7.json"

    # condition_info = {
    #     "Day": "Day_2",
    #     "Experiment": "Experiment_1",
    #     "SkinTone": "Fair",
    #     "Speed": "Intermediate",
    #     "Depth": "3.75mm",
    #     "ExpectedBPM": 90,
    #     "Clamp": "Yes",
    #     "PolarizationPlacement": "Same",
    #     "ConditionFolder": "3.75 Fair Intermediate",
    # }

#     debug_one_file_one_channel(
#         json_path=test_json,
#         condition_info=condition_info,
#         channel_label="Cross-Polarized_IR"
#     )