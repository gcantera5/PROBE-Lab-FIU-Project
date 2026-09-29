"""
One-page analysis report for a single measurement (one recording x one channel).

Built for Mark's request (Sept 2026):
  1. Measurement condition (skin tone, depth, flow speed/BPM, wavelength,
     polarization, ...)
  2. PI for the measurement (mean and SD across windows)
  3. Bad window rate for each SQI (plus the overall rate)
  4. A reference PPG waveform: one 10-second window and its average beat

Usage (from the FIU Project folder):
    python measurement_report.py <recording id or part of the file name> <channel>

    e.g.  python measurement_report.py 02-09-29 Co-Polarized_IR

Outputs go to FIU_Reports/:
    <day>_<condition>_<channel>_<file id>.png / .pdf   the report page
    report_numbers.csv                                  one row per report
                                                        (appended), so reports
                                                        can be compared later

The report re-runs the same pipeline functions as beat_level_sqi.py (same
config), so its numbers match the saved window CSVs.
"""

import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd

import beat_level_sqi as sqi


# ------------------------------------------------------------
# Where things live
# ------------------------------------------------------------

SUMMARY_FILES = [
    "FIU_Beat_Level_SQI/Day_2/Experiment_1/day2_experiment1_recording_summary.csv",
    "FIU_Beat_Level_SQI/Day_3/Experiment_2/day3_experiment2_recording_summary.csv",
    "FIU_Beat_Level_SQI/Day_4/day4_recording_summary.csv",
]

DATA_ROOTS = [
    "Experiment 1 Complete  copy",
    "Experiment 2 Test (Day 3) copy",
    "Experiment 2 & 3 (Day 4) copy",
]

OUTPUT_DIR = "FIU_Reports"

# The three checks that can reject a beat, and the reason tag
# classify_beat() writes for each one.
SQI_CHECKS = [
    ("Template SQI", "low_template_sqi"),
    ("Correlation", "low_correlation"),
    ("Clipping", "clipping"),
]

# Colors (reference data-viz palette)
BLUE = "#2a78d6"        # main waveform / bars
GRAY = "#b8b6b0"        # individual beats, secondary marks
INK = "#0b0b0b"         # primary text
INK_2 = "#52514e"       # secondary text
GRID = "#e4e2dc"
GOOD = "#0ca30c"        # status: good window
BAD = "#d03b3b"         # status: bad window


# ------------------------------------------------------------
# Finding the measurement
# ------------------------------------------------------------

def find_measurement(file_id, channel):
    """
    Look up one recording x channel in the saved recording summaries, so the
    report uses the exact condition labels the pipeline saved.
    Returns (condition_row, json_path).
    """
    rows = []
    for path in SUMMARY_FILES:
        if os.path.exists(path):
            df = pd.read_csv(path)
            rows.append(df[df["SourceFile"].str.contains(file_id, regex=False)
                           & (df["Channel"] == channel)])

    matches = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()

    if len(matches) == 0:
        raise ValueError(f"No saved summary row for file '{file_id}' + channel '{channel}'. "
                         f"Run beat_level_sqi.py first, and check the channel name.")
    if matches["SourceFile"].nunique() > 1:
        raise ValueError(f"'{file_id}' matches more than one recording: "
                         f"{matches['SourceFile'].unique().tolist()}")

    condition = matches.iloc[0]

    json_paths = []
    for root in DATA_ROOTS:
        json_paths += glob.glob(os.path.join(root, "**", condition["SourceFile"]), recursive=True)
    if not json_paths:
        raise FileNotFoundError(f"Couldn't find {condition['SourceFile']} under {DATA_ROOTS}")

    return condition, json_paths[0]


# ------------------------------------------------------------
# Computing the numbers
# ------------------------------------------------------------

def analyze_measurement(json_path, channel, config):
    """
    Run the pipeline window by window and keep what the report needs:
    per-window PI and labels, per-SQI window failures, beat HRs, and the
    per-window detail needed to plot a reference window.
    """
    df = sqi.load_fiu_json(json_path, {})
    signal = df[channel].values.astype(float)

    windows = []
    details = {}
    beat_hrs = []

    for start, end, raw_window in sqi.iter_windows(signal, config):
        result = sqi.run_window_level_sqi(raw_window, config)
        if result is None:
            continue

        beats = result["feature_table"]
        row = {
            "start_sec": start / config.fs,
            "perfusion_index": result["row"]["perfusion_index"],
            "mean_correlation": result["row"]["mean_correlation"],
            "window_label": result["row"]["window_label"],
            "num_beats": result["row"]["num_beats"],
        }

        # A window is "bad by check X" if more than the window threshold of
        # its beats fail check X -- same roll-up rule as the overall label,
        # applied to one check at a time.
        for _, tag in SQI_CHECKS:
            fails = beats["rejection_reasons"].str.contains(tag, regex=False)
            row[f"bad_{tag}"] = fails.mean() > config.bad_window_fraction_threshold

        windows.append(row)
        details[start] = (raw_window, result)
        beat_hrs.extend(beats["estimated_hr"].tolist())

    windows = pd.DataFrame(windows)
    windows["start_idx"] = (windows["start_sec"] * config.fs).round().astype(int)

    return signal, windows, details, np.array(beat_hrs)


def pick_reference_window(windows):
    """
    Reference window = a TYPICAL good window: among the good windows, the one
    whose mean beat correlation is closest to the median of the good windows.
    If there are no good windows, fall back to the best (highest-correlation)
    window, and say so on the report.
    """
    good = windows[windows["window_label"] == "good_window"]
    if len(good) > 0:
        target = good["mean_correlation"].median()
        pick = good.iloc[(good["mean_correlation"] - target).abs().argmin()]
        return int(pick["start_idx"]), "typical good window"

    pick = windows.iloc[windows["mean_correlation"].argmax()]
    return int(pick["start_idx"]), "no good windows -- best available window shown"


# ------------------------------------------------------------
# Report page
# ------------------------------------------------------------

def _style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(GRID)
    ax.tick_params(colors=INK_2, labelsize=8, length=0)
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def _fmt(value, spec=".2f", suffix=""):
    return "--" if pd.isna(value) else f"{value:{spec}}{suffix}"


def condition_table(condition, channel, json_path, signal, config, measured_bpm):
    pol_group, led = channel.rsplit("_", 1)
    led_nm = sqi.WAVELENGTH_NM.get(led)

    def get(col):
        return condition[col] if col in condition.index and not pd.isna(condition[col]) else None

    rows = [
        ("Day / experiment", f"{str(get('Day')).replace('_', ' ')}  /  "
                             f"{str(get('Experiment')).replace('Experiment_', 'Exp ').replace(' & Exp ', ' & ')}"),
        ("Skin tone", get("SkinTone")),
        ("Vessel depth", get("Depth")),
        ("Flow speed", get("Speed")),
        ("BPM (labeled / measured)", f"{get('ExpectedBPM')}  /  {measured_bpm:.0f}"),
        ("Phantom wavelength", f"{get('Wavelength')} ({get('WavelengthNm'):.0f} nm)"
                                if get("Wavelength") else "n/a (Experiment 1)"),
        ("LED channel", f"{led} ({led_nm} nm)"),
        ("Polarization channel", pol_group.replace("_", " ")),
        ("Polarization setup", get("PolarizationCondition") or get("PolarizationPlacement")),
        ("Device orientation", f"{get('OrientationDegrees'):.0f} deg" if get("OrientationDegrees") is not None else None),
        ("Clamp / phantom", " / ".join(str(x) for x in [get("Clamp"), get("PhantomType")] if x)),
        ("Hardware channel", get("Hardware Channel")),
        ("Recording", f"{os.path.basename(json_path)[:19]}  (trial {get('Trial')}, {len(signal) / config.fs:.0f} s)"),
    ]
    return [(k, str(v)) for k, v in rows if v not in (None, "", "None")]


def build_report(condition, json_path, channel, out_dir=OUTPUT_DIR, config=None):
    config = config or sqi.BeatSQIConfig()
    signal, windows, details, beat_hrs = analyze_measurement(json_path, channel, config)

    if len(windows) == 0:
        raise ValueError("No window in this measurement had enough beats to analyze.")

    # ---- numbers ----
    pi = windows["perfusion_index"]
    measured_bpm = float(np.nanmedian(beat_hrs))
    rates = [("Overall", 100 * (windows["window_label"] == "bad_window").mean())]
    rates += [(name, 100 * windows[f"bad_{tag}"].mean()) for name, tag in SQI_CHECKS]

    ref_start, ref_note = pick_reference_window(windows)
    ref_raw, ref_result = details[ref_start]
    ref_label = windows.loc[windows["start_idx"] == ref_start, "window_label"].iloc[0]

    # ---- page ----
    fig = plt.figure(figsize=(11, 8.5), facecolor="white")
    gs = GridSpec(3, 3, figure=fig, height_ratios=[1.0, 0.95, 1.1],
                  hspace=0.6, wspace=0.34, left=0.08, right=0.97, top=0.87, bottom=0.1)

    title = f"{condition['ConditionFolder'].strip()}  |  {channel.replace('_', ' ')}"
    fig.text(0.08, 0.955, title, fontsize=15, weight="bold", color=INK)
    fig.text(0.08, 0.928, f"Measurement report  -  {len(windows)} windows analyzed "
             f"(10 s, 1 s step)  -  window is bad if > {config.bad_window_fraction_threshold:.0%} of its beats fail",
             fontsize=9, color=INK_2)

    # 1. Conditions
    ax = fig.add_subplot(gs[0, 0:2]); ax.axis("off")
    ax.set_title("1. Measurement condition", loc="left", fontsize=11, weight="bold", color=INK)
    rows = condition_table(condition, channel, json_path, signal, config, measured_bpm)
    half = (len(rows) + 1) // 2
    for col, chunk in enumerate([rows[:half], rows[half:]]):
        for i, (k, v) in enumerate(chunk):
            y = 0.9 - i * 0.135
            ax.text(col * 0.52, y, k, fontsize=8.5, color=INK_2, transform=ax.transAxes, va="top")
            ax.text(col * 0.52 + 0.255, y, v, fontsize=8.5, color=INK, transform=ax.transAxes, va="top")

    # 2. PI headline
    ax = fig.add_subplot(gs[0, 2]); ax.axis("off")
    ax.set_title("2. Perfusion index", loc="left", fontsize=11, weight="bold", color=INK)
    ax.text(0, 0.72, f"{_fmt(pi.mean())} %", fontsize=26, weight="bold", color=INK, transform=ax.transAxes)
    ax.text(0, 0.56, f"mean  ±  {_fmt(pi.std())} % SD  (across windows)", fontsize=9, color=INK_2, transform=ax.transAxes)
    ax.text(0, 0.40, f"median {_fmt(pi.median())} %   range {_fmt(pi.min())} - {_fmt(pi.max())} %",
            fontsize=9, color=INK_2, transform=ax.transAxes)
    ax.text(0, 0.18, "Window PI = median of its beat PIs\n(raw signal, AC/DC x 100, DC = beat trough)",
            fontsize=7.5, color=INK_2, transform=ax.transAxes, va="top")

    # 3. Bad window rates
    ax = fig.add_subplot(gs[1, 0])
    names = [n for n, _ in rates][::-1]
    vals = [v for _, v in rates][::-1]
    bars = ax.barh(names, vals, color=[BLUE] * (len(vals) - 1) + [INK_2], height=0.55)
    for bar, v in zip(bars, vals):
        inside = v > 85
        ax.text(v - 2 if inside else v + 2, bar.get_y() + bar.get_height() / 2, f"{v:.0f}%",
                va="center", ha="right" if inside else "left", fontsize=9,
                color="white" if inside else INK)
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of windows bad", fontsize=8, color=INK_2)
    _style(ax); ax.grid(axis="y", visible=False)
    ax.set_title("3. Bad window rate", loc="left", fontsize=11, weight="bold", color=INK)

    # PI + window label over the recording
    ax = fig.add_subplot(gs[1, 1:3])
    ax.plot(windows["start_sec"], pi, color=BLUE, linewidth=1.6)
    ymin, ymax = np.nanmin(pi), np.nanmax(pi)
    pad = 0.1 * (ymax - ymin if ymax > ymin else 1)
    strip = ymin - pad * 2.5
    colors = np.where(windows["window_label"] == "good_window", GOOD, BAD)
    ax.scatter(windows["start_sec"], np.full(len(windows), strip), c=colors, marker="s", s=14, linewidths=0)
    ax.axvspan(ref_start / config.fs, ref_start / config.fs + config.window_seconds, color=BLUE, alpha=0.08)
    ax.set_ylim(strip - pad, ymax + pad)
    ax.set_xlabel("window start (s)", fontsize=8, color=INK_2)
    ax.set_ylabel("PI (%)", fontsize=8, color=INK_2)
    _style(ax)
    ax.set_title("PI per window  (strip below = window label;  shaded = reference window)",
                 loc="left", fontsize=9, color=INK_2)
    # legend chips for the strip, so good/bad isn't color-only
    ax.scatter([], [], c=GOOD, marker="s", s=20, label="good window")
    ax.scatter([], [], c=BAD, marker="s", s=20, label="bad window")
    ax.legend(loc="upper right", fontsize=7.5, frameon=False, ncol=2)

    # 4. Reference waveform: raw window, filtered window, average beat
    t = np.arange(len(ref_raw)) / config.fs
    beats = ref_result["feature_table"]

    ax = fig.add_subplot(gs[2, 0])
    ax.plot(t, ref_raw, color=INK_2, linewidth=1.2)
    _style(ax)
    ax.set_xlabel("time in window (s)", fontsize=8, color=INK_2)
    ax.set_ylabel("raw ADC counts", fontsize=8, color=INK_2)
    ax.set_title("4. Reference window: raw", loc="left", fontsize=11, weight="bold", color=INK)

    ax = fig.add_subplot(gs[2, 1])
    filt = ref_result["filtered_signal"]
    ax.plot(t, filt, color=BLUE, linewidth=1.4)
    for _, b in beats.iterrows():
        ax.axvspan(b["beat_start_sec"], b["beat_end_sec"],
                   color=GOOD if b["beat_label"] == "good" else BAD, alpha=0.10, linewidth=0)
    v = beats["beat_start_idx"].astype(int).tolist() + [int(beats["beat_end_idx"].iloc[-1])]
    ax.scatter(np.array(v) / config.fs, filt[v], color=INK, marker="v", s=14, zorder=3)
    _style(ax)
    n_good = int((beats["beat_label"] == "good").sum())
    ax.set_xlabel("time in window (s)", fontsize=8, color=INK_2)
    ax.set_title(f"filtered + beats ({n_good}/{len(beats)} good)", loc="left", fontsize=9, color=INK_2)

    ax = fig.add_subplot(gs[2, 2])
    x = np.linspace(0, 100, ref_result["normalized_beats"].shape[1])
    for nb in ref_result["normalized_beats"]:
        ax.plot(x, nb, color=GRAY, linewidth=0.8)
    ax.plot(x, ref_result["template"], color=BLUE, linewidth=2.2, label="average beat (template)")
    ax.plot([], [], color=GRAY, linewidth=0.8, label="individual beats")
    ax.legend(loc="upper right", fontsize=7, frameon=False)
    _style(ax)
    ax.set_ylim(-0.05, 1.25)
    ax.set_xlabel("% of beat", fontsize=8, color=INK_2)
    ax.set_ylabel("normalized (0-1)", fontsize=8, color=INK_2)
    ax.set_title("average beat", loc="left", fontsize=9, color=INK_2)

    fig.text(0.08, 0.035,
             f"Reference window: {ref_start / config.fs:.0f}-{ref_start / config.fs + config.window_seconds:.0f} s "
             f"({ref_note}; labeled {ref_label.replace('_', ' ')}).  "
             f"A beat fails if template SQI < {config.min_template_sqi}, correlation < {config.min_corr}, "
             f"or clipping SQI < {config.min_clipping_sqi}.",
             fontsize=7.5, color=INK_2)
    fig.text(0.08, 0.015, f"Source: {os.path.basename(json_path)}", fontsize=7.5, color=INK_2)

    # ---- save ----
    os.makedirs(out_dir, exist_ok=True)
    file_id = os.path.basename(json_path)[11:19]
    stem = f"{condition['Day']}_{condition['ConditionFolder'].strip().replace(' ', '_').replace('.', '')}_{channel}_{file_id}"
    png = os.path.join(out_dir, stem + ".png")
    pdf = os.path.join(out_dir, stem + ".pdf")
    fig.savefig(png, dpi=150)
    fig.savefig(pdf)
    plt.close(fig)

    numbers = {
        "Report": stem,
        "SourceFile": os.path.basename(json_path),
        "Channel": channel,
        **{k: condition[k] for k in condition.index
           if k in ["Day", "Experiment", "ConditionFolder", "Trial", "SkinTone", "Depth", "Speed",
                    "ExpectedBPM", "Wavelength", "WavelengthNm", "PolarizationCondition",
                    "PolarizationPlacement", "OrientationDegrees", "Clamp", "PhantomType",
                    "Hardware Channel"]},
        "MeasuredBPM": measured_bpm,
        "Windows": len(windows),
        "PI_mean": pi.mean(), "PI_sd": pi.std(), "PI_median": pi.median(),
        **{f"BadWindowPct_{n.split(' ')[0]}": v for n, v in rates},
        "ReferenceWindowStartSec": ref_start / config.fs,
    }
    numbers_path = os.path.join(out_dir, "report_numbers.csv")
    table = pd.DataFrame([numbers])
    if os.path.exists(numbers_path):
        old = pd.read_csv(numbers_path)
        table = pd.concat([old[old["Report"] != stem], table], ignore_index=True)
    table.to_csv(numbers_path, index=False)

    return png, pdf, numbers


def generate_measurement_report(file_id, channel, out_dir=OUTPUT_DIR):
    condition, json_path = find_measurement(file_id, channel)
    png, pdf, numbers = build_report(condition, json_path, channel, out_dir)
    print(f"Saved {png}\nSaved {pdf}")
    return png, pdf, numbers


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    generate_measurement_report(sys.argv[1], sys.argv[2])
