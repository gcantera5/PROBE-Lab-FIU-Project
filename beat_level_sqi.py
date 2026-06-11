import os
import glob
import json
import gzip
import shutil
from dataclasses import dataclass
from typing import List, Tuple, Optional

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


# ============================================================
# EXPERIMENT 1 PROCESSING
# ============================================================

def process_experiment1_complete(
    experiment_root="Experiment 1 Complete  copy",
    output_root="FIU_Beat_Level_SQI/Experiment_1"
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
            "Day": "Day_1",
            "SkinTone": skin,
            "Speed": speed,
            "Depth": depth,
            "Experiment": "Experiment_2",
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

        json_files = glob.glob(os.path.join(folder, "*.json"))

        for json_path in json_files:

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

                        print(f"Good Windows: {good_windows}/{total_windows}")
                        print("--------------------------------")
                        
                except Exception as e:
                    print(f"SQI failed for {col}: {e}")
                    continue

                if len(window_df) > 0:
                    window_df["SourceFile"] = os.path.basename(json_path)
                    all_window_results.append(window_df)

                if len(beat_df) > 0:
                    beat_df["SourceFile"] = os.path.basename(json_path)
                    all_beat_results.append(beat_df)

    if all_window_results:
        final_window_df = pd.concat(all_window_results, ignore_index=True)
        window_path = os.path.join(output_root, "experiment_1_all_window_sqi.csv")
        final_window_df.to_csv(window_path, index=False)
        print(f"\nSaved window-level SQI results to: {window_path}")

    if all_beat_results:
        final_beat_df = pd.concat(all_beat_results, ignore_index=True)
        beat_path = os.path.join(output_root, "experiment_1_all_beat_sqi.csv")
        final_beat_df.to_csv(beat_path, index=False)
        print(f"Saved beat-level SQI results to: {beat_path}")


# ============================================================
# MAIN BLOCK
# ============================================================

if __name__ == "__main__":

    # Change this if your folder name/path is different
    EXPERIMENT_1_FOLDER = "Experiment 1 Complete  copy"

    process_experiment1_complete(
        experiment_root=EXPERIMENT_1_FOLDER,
        output_root="FIU_Beat_Level_SQI/Experiment_1"
    )