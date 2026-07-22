import os
import glob
import json
import gzip
import shutil
from dataclasses import dataclass
from typing import List, Tuple, Optional
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
from scipy.signal import cheby2, sosfiltfilt, find_peaks, medfilt, resample
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

WAVELENGTH_NM = {
    "Green": 525,
    "Red": 660,
    "IR": 940,
}


@dataclass
class BeatSQIConfig:
    fs: int = FS

    # GT recommended removing low-frequency components
    remove_baseline: bool = True
    baseline_kernel_seconds: float = 0.75

    # FIU-style filtering
    lowcut: float = 0.5
    highcut: float = 2.2
    filter_order: int = 4
    rs: int = 20

    # Beat segmentation
    min_hr: float = 40
    max_hr: float = 180
    valley_prominence: Optional[float] = None

    # Beat normalization/template
    target_beat_length: int = 50
    n_template_beats: int = 12
    template_corr_threshold: float = 0.80

    # SQI thresholds
    sqi_lambda: float = 25.0
    min_template_sqi: float = 0.05
    min_corr: float = 0.80
    max_mad: float = 0.60
    min_clipping_sqi: float = 0.80

    # Window rejection
    window_seconds: int = 10
    step_seconds: int = 1
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

def remove_low_frequency_baseline(ppg, config):
    """
    Removes slow baseline drift before beat detection.
    Removes low-frequency components.
    """
    ppg = np.asarray(ppg, dtype=float)

    if not config.remove_baseline:
        return ppg

    kernel = int(round(config.baseline_kernel_seconds * config.fs))

    if kernel < 3:
        return ppg

    if kernel % 2 == 0:
        kernel += 1

    baseline = medfilt(ppg, kernel_size=kernel)
    return ppg - baseline


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
    Full preprocessing:
    1. remove low-frequency baseline drift
    2. bandpass filter
    3. zero-center
    """
    no_baseline = remove_low_frequency_baseline(ppg, config)
    filtered = bandpass_filter(no_baseline, config)
    centered = filtered - np.nanmean(filtered)
    return np.nan_to_num(centered, nan=0.0)


# ============================================================
# BEAT SEGMENTATION
# ============================================================

def segment_beats_by_valleys(ppg, config):
    """
    Extract non-overlapping beats using valley-to-valley intervals.
    Beat-level workflow.
    """
    min_distance = int(config.fs * 60.0 / config.max_hr)

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
    Normalize amplitude and resample beat to uniform length.
    """
    beat = np.asarray(beat, dtype=float)
    beat = beat - np.mean(beat)

    std = np.std(beat)
    if std > 0:
        beat = beat / std

    return resample(beat, target_length)


def normalize_beats(raw_beats, config):
    """Normalize all beats to the same length."""
    return np.array([
        normalize_single_beat(beat, config.target_beat_length)
        for beat in raw_beats
    ])


def make_clean_template(norm_beats, config):
    """
    Build a clean template beat by averaging around 12 template-like beats.
    """
    if len(norm_beats) == 0:
        raise ValueError("No beats available to build template.")

    n_initial = min(config.n_template_beats, len(norm_beats))
    rough_template = np.mean(norm_beats[:n_initial], axis=0)

    correlations = np.array([
        safe_corr(beat, rough_template)
        for beat in norm_beats
    ])

    good_idx = np.where(correlations >= config.template_corr_threshold)[0]

    if len(good_idx) < min(3, len(norm_beats)):
        good_idx = np.argsort(correlations)[::-1]

    selected_idx = good_idx[:min(config.n_template_beats, len(good_idx))]
    template = np.mean(norm_beats[selected_idx], axis=0)

    return template, selected_idx


# ============================================================
# BEAT-LEVEL FEATURES
# ============================================================

def compute_dtw_distance(beat, template):
    """
    Compare beat to template with DTW.
    """
    if dtw is None:
        raise ImportError("Please install dtw-python using: pip install dtw-python")

    alignment = dtw(beat, template, keep_internals=False)
    return float(alignment.normalizedDistance)


def compute_template_sqi(dtw_distance, config):
    """
    Convert DTW distance into SQI.
    Small distance = high SQI.
    Large distance = low SQI.
    """
    return float(np.exp(-config.sqi_lambda * dtw_distance))


def clipping_sqi(raw_beat):
    """
    Estimate how much of a beat is clipped/saturated.
    1.0 = not clipped.
    Lower values = more clipping.
    """
    raw_beat = np.asarray(raw_beat, dtype=float)

    beat_range = np.max(raw_beat) - np.min(raw_beat)

    if beat_range == 0:
        return 0.0

    eps = 0.01 * beat_range

    near_max = raw_beat >= (np.max(raw_beat) - eps)
    near_min = raw_beat <= (np.min(raw_beat) + eps)

    clipped_fraction = np.mean(near_max | near_min)
    return float(1.0 - clipped_fraction)


def compute_beat_features(raw_beat, norm_beat, template, start_idx, end_idx, beat_number, config):
    """
    Compute beat-level SQI features:
    DTW, correlation, AD, MAD, skewness, clipping SQI, duration, HR.
    """
    difference = norm_beat - template

    dtw_distance = compute_dtw_distance(norm_beat, template)
    template_sqi = compute_template_sqi(dtw_distance, config)

    corr = max(0.0, safe_corr(norm_beat, template))
    mad = float(np.mean(np.abs(difference)))
    ad = float(np.sum(np.abs(difference)))
    beat_skewness = float(skew(norm_beat, nan_policy="omit"))
    clip_sqi = clipping_sqi(raw_beat)

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
        "AD": ad,
        "MAD": mad,
        "skewness": beat_skewness,
        "clipping_sqi": clip_sqi,
    }


# def classify_beat(row, config):
#     """
#     First-pass rule-based beat rejection.
#     These thresholds can be tuned after visual inspection.
#     """
#     is_bad = (
#         row["template_sqi"] < config.min_template_sqi
#         or row["correlation"] < config.min_corr
#         or row["MAD"] > config.max_mad
#         or row["clipping_sqi"] < config.min_clipping_sqi
#     )

#     return "bad" if is_bad else "good"

def classify_beat(row, config):
    """
    First-pass rule-based beat rejection.
    Also records why a beat was rejected.
    """
    reasons = []

    if row["template_sqi"] < config.min_template_sqi:
        reasons.append("low_template_sqi")

    if row["correlation"] < config.min_corr:
        reasons.append("low_correlation")

    if row["MAD"] > config.max_mad:
        reasons.append("high_MAD")

    if row["clipping_sqi"] < config.min_clipping_sqi:
        reasons.append("clipping")

    label = "bad" if reasons else "good"

    return label, ",".join(reasons)


# ============================================================
# MAIN BEAT-LEVEL PIPELINE
# ============================================================

def run_beat_level_sqi(ppg, config=None):
    """
    Run beat-level SQI on one PPG segment/window.
    """
    if config is None:
        config = BeatSQIConfig()

    filtered = preprocess_ppg(ppg, config)

    raw_beats, beat_indices = segment_beats_by_valleys(filtered, config)

    if len(raw_beats) < 3:
        return None

    norm_beats = normalize_beats(raw_beats, config)

    template, template_indices = make_clean_template(norm_beats, config)

    rows = []

    for beat_number, (raw_beat, norm_beat, (start_idx, end_idx)) in enumerate(
        zip(raw_beats, norm_beats, beat_indices)
    ):
        row = compute_beat_features(
            raw_beat=raw_beat,
            norm_beat=norm_beat,
            template=template,
            start_idx=start_idx,
            end_idx=end_idx,
            beat_number=beat_number,
            config=config
        )

        #row["beat_label"] = classify_beat(row, config)
        label, reasons = classify_beat(row, config)
        row["beat_label"] = label
        row["rejection_reasons"] = reasons
        row["used_for_template"] = beat_number in set(template_indices.tolist())

        rows.append(row)

    feature_table = pd.DataFrame(rows)

    percent_bad = np.mean(feature_table["beat_label"] == "bad")

    window_label = (
        "bad_window"
        if percent_bad > config.bad_window_fraction_threshold
        else "good_window"
    )

    summary = {
        "num_beats": len(feature_table),
        "num_good_beats": int(np.sum(feature_table["beat_label"] == "good")),
        "num_bad_beats": int(np.sum(feature_table["beat_label"] == "bad")),
        "percent_bad": float(percent_bad),
        "mean_template_sqi": float(feature_table["template_sqi"].mean()),
        "mean_dtw_distance": float(feature_table["dtw_distance"].mean()),
        "mean_correlation": float(feature_table["correlation"].mean()),
        "mean_MAD": float(feature_table["MAD"].mean()),
        "window_label": window_label,
    }

    return {
        "filtered_signal": filtered,
        "raw_beats": raw_beats,
        "normalized_beats": norm_beats,
        "template": template,
        "feature_table": feature_table,
        "summary": summary,
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
    Run beat-level SQI over 10-second FIU-style windows.
    """
    if config is None:
        config = BeatSQIConfig()

    window_rows = []
    beat_rows = []

    for win_start, win_end, win in iter_windows(signal, config):
        result = run_beat_level_sqi(win, config)

        if result is None:
            continue

        summary = result["summary"]

        window_row = {
            "Channel": channel_label,
            "Hardware Channel": channel_name_map(channel_label),
            "WindowStartIdx": win_start,
            "WindowEndIdx": win_end,
            "WindowStartSec": win_start / config.fs,
            "WindowEndSec": win_end / config.fs,
            **summary,
            **condition_info,
        }

        window_rows.append(window_row)

        beats = result["feature_table"].copy()
        beats["Channel"] = channel_label
        beats["Hardware Channel"] = channel_name_map(channel_label)

        # convert beat indices from local window index to full-signal index
        beats["beat_start_idx_global"] = beats["beat_start_idx"] + win_start
        beats["beat_end_idx_global"] = beats["beat_end_idx"] + win_start
        beats["beat_start_sec_global"] = beats["beat_start_idx_global"] / config.fs
        beats["beat_end_sec_global"] = beats["beat_end_idx_global"] / config.fs

        for key, val in condition_info.items():
            beats[key] = val

        beat_rows.append(beats)

    window_df = pd.DataFrame(window_rows)

    if beat_rows:
        beat_df = pd.concat(beat_rows, ignore_index=True)
    else:
        beat_df = pd.DataFrame()

    return window_df, beat_df

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
    skin = parts[1].capitalize()
    orientation = parts[2]
    pol = " ".join(parts[3:])

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
    skin = parts[1].capitalize()

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
    all_beat_results = []
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
                    window_df, beat_df = run_sqi_over_windows(
                        signal=signal,
                        condition_info=condition_info,
                        channel_label=col,
                        config=config
                    )
                    if len(window_df) > 0:
                        print("\n--------------------------------")
                        print(f"Condition: {folder_name}")
                        print(f"Channel: {col}")
                        print(f"Mean DTW: {window_df['mean_dtw_distance'].mean():.4f}")
                        print(f"Mean Correlation: {window_df['mean_correlation'].mean():.4f}")
                        print(f"Mean MAD: {window_df['mean_MAD'].mean():.4f}")

                        good_windows = (window_df["window_label"] == "good_window").sum()
                        total_windows = len(window_df)

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
                            "TotalWindows": int(total_windows),
                            "GoodWindowFraction": float(
                                good_window_fraction
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
                        }

                        all_summary_results.append(summary_row)
                        
                except Exception as e:
                    print(f"SQI failed for {col}: {e}")
                    continue

                if len(window_df) > 0:
                    window_df["SourceFile"] = os.path.basename(json_path)
                    window_df["Trial"] = trial_number
                    all_window_results.append(window_df)

                if len(beat_df) > 0:
                    beat_df["SourceFile"] = os.path.basename(json_path)
                    beat_df["Trial"] = trial_number
                    all_beat_results.append(beat_df)

    if all_window_results:
        final_window_df = pd.concat(all_window_results, ignore_index=True)
        window_path = os.path.join(output_root, "day2_experiment1_all_window_sqi.csv")
        final_window_df.to_csv(window_path, index=False)
        print(f"\nSaved window-level SQI results to: {window_path}")

    if all_beat_results:
        final_beat_df = pd.concat(all_beat_results, ignore_index=True)
        beat_path = os.path.join(output_root, "day2_experiment1_all_beat_sqi.csv")
        final_beat_df.to_csv(beat_path, index=False)
        print(f"Saved beat-level SQI results to: {beat_path}")
    
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
    all_beat_results = []
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
                    window_df, beat_df = run_sqi_over_windows(
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
                    "TotalWindows": int(total_windows),
                    "GoodWindowFraction": float(
                        good_window_fraction
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
                }

                all_summary_results.append(summary_row)

                window_df["SourceFile"] = source_file
                window_df["Trial"] = trial_number
                all_window_results.append(window_df)

                if len(beat_df) > 0:
                    beat_df["SourceFile"] = source_file
                    beat_df["Trial"] = trial_number
                    all_beat_results.append(beat_df)

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
    # SAVE BEAT RESULTS
    # --------------------------------------------------
    if all_beat_results:
        final_beat_df = pd.concat(
            all_beat_results,
            ignore_index=True
        )

        beat_path = os.path.join(
            output_root,
            "day3_experiment2_all_beat_sqi.csv"
        )

        final_beat_df.to_csv(
            beat_path,
            index=False
        )

        print(
            f"Saved Day 3 beat results to: "
            f"{beat_path}"
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
    all_beat_results = []
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
                        window_df, beat_df = run_sqi_over_windows(
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
                            print(f"Mean DTW: {window_df['mean_dtw_distance'].mean():.4f}")
                            print(f"Mean Correlation: {window_df['mean_correlation'].mean():.4f}")
                            print(f"Mean MAD: {window_df['mean_MAD'].mean():.4f}")

                            good_windows = (
                                window_df["window_label"] == "good_window"
                            ).sum()

                            total_windows = len(window_df)
                            
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
                                "TotalWindows": int(total_windows),
                                "GoodWindowFraction": float(
                                    good_window_fraction
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
                            }

                            all_summary_results.append(summary_row)


                    except Exception as e:
                        print(f"SQI failed for {col}: {e}")
                        continue

                    if len(window_df) > 0:
                        window_df["SourceFile"] = os.path.basename(json_path)
                        window_df["Trial"] = trial_number
                        all_window_results.append(window_df)

                    if len(beat_df) > 0:
                        beat_df["SourceFile"] = os.path.basename(json_path)
                        beat_df["Trial"] = trial_number
                        all_beat_results.append(beat_df)


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

        final_window_df.to_csv(
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
    # SAVE BEAT-LEVEL RESULTS
    # --------------------------------------------------
    if all_beat_results:
        final_beat_df = pd.concat(
            all_beat_results,
            ignore_index=True
        )

        # Save one combined Day 4 file
        beat_path = os.path.join(
            output_root,
            "day4_all_beat_sqi.csv"
        )

        final_beat_df.to_csv(
            beat_path,
            index=False
        )

        print(
            f"Saved Day 4 beat-level SQI results to: "
            f"{beat_path}"
        )

        # Save Experiment 2 and Experiment 3 separately
        for experiment_name in ["Experiment_2", "Experiment_3"]:

            experiment_beat_df = final_beat_df[
                final_beat_df["Experiment"] == experiment_name
            ]

            if len(experiment_beat_df) > 0:

                experiment_folder = os.path.join(
                    output_root,
                    experiment_name
                )

                os.makedirs(
                    experiment_folder,
                    exist_ok=True
                )

                experiment_beat_path = os.path.join(
                    experiment_folder,
                    f"day4_{experiment_name.lower()}_beat_sqi.csv"
                )

                experiment_beat_df.to_csv(
                    experiment_beat_path,
                    index=False
                )

                print(
                    f"Saved {experiment_name} beat results to: "
                    f"{experiment_beat_path}"
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

        summary_df.to_csv(
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
    1. detected valleys
    2. normalized beats + template
    3. good/bad beat labels
    """
    config = BeatSQIConfig()

    cleaned_df = load_fiu_json(json_path, condition_info)

    if channel_label not in cleaned_df.columns:
        print(f"{channel_label} not found.")
        print("Available columns:")
        print(cleaned_df.columns.tolist())
        return

    signal = cleaned_df[channel_label].values

    result = run_beat_level_sqi(signal, config=config)

    if result is None:
        print("No result. Not enough beats were detected.")
        return

    filtered = result["filtered_signal"]
    feature_table = result["feature_table"]
    norm_beats = result["normalized_beats"]
    template = result["template"]

    print("\n===== DEBUG SUMMARY =====")
    print(f"Channel: {channel_label}")
    print(result["summary"])
    print(feature_table[[
        "beat_number",
        "dtw_distance",
        "correlation",
        "MAD",
        "template_sqi",
        "clipping_sqi",
        "beat_label",
        "rejection_reasons"
    ]])

    # Plot 1: filtered PPG with detected valleys
    t = np.arange(len(filtered)) / config.fs
    valley_idxs = feature_table["beat_start_idx"].values.astype(int)

    plt.figure(figsize=(10, 4))
    plt.plot(t, filtered, linewidth=1.2)
    plt.scatter(
        valley_idxs / config.fs,
        filtered[valley_idxs],
        marker="v",
        s=50,
        label="Detected valleys"
    )
    plt.title(f"{channel_label}: Filtered PPG with Detected Valleys")
    plt.xlabel("Time (s)")
    plt.ylabel("Filtered PPG")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    plt.show()

    # Plot 2: normalized beats + template
    plt.figure(figsize=(8, 4))
    for beat in norm_beats:
        plt.plot(beat, alpha=0.25, linewidth=1)

    plt.plot(template, linewidth=3, label="Template beat")
    plt.title(f"{channel_label}: Normalized Beats + Template")
    plt.xlabel("Normalized sample")
    plt.ylabel("Normalized amplitude")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    plt.show()

    # Plot 3: good/bad beat intervals
    plt.figure(figsize=(10, 4))
    plt.plot(t, filtered, linewidth=1.2, color="black")

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

    # fake legend handles
    plt.plot([], [], color="green", linewidth=8, alpha=0.4, label="Good beat")
    plt.plot([], [], color="red", linewidth=8, alpha=0.4, label="Bad beat")
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