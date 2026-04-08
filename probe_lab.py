import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, butter, filtfilt
import os


"""
- do CVS pipeline for mean and std data 
- have the data saved seperately for stat analysis
- lab field data??? 
- rewrite names to be participant #s -- no ID check
"""

# ---------------------------------
# Helper Functions
# ---------------------------------

def lowpass_filter(signal, cutoff_hz, fs, order=4):
    """Apply a low-pass Butterworth filter (replicates MATLAB lowpass())."""
    b, a = butter(order, cutoff_hz / (fs / 2), btype='low')
    return filtfilt(b, a, signal)

def compute_PI(peaks, troughs):
    """Compute Perfusion Index as (AC/DC)*100, sign-corrected."""
    AC = peaks + troughs
    PI = -(AC / troughs) * 100
    return PI

# ---------------------------------
# Load Data
# ---------------------------------

# Replace with your actual file path (CSV or Excel)
file_path = 'Updated Channels Excel/JuliaDataNewChannels.xlsx - Sheet1.csv'
data = pd.read_csv(file_path)
data = data.apply(pd.to_numeric, errors='coerce')

# Extract participant name based on filename convention 
# (e.g., ___iDataNewChannels.xlsx)
base_name = os.path.basename(file_path)

# Extract Name on file
participant_name = base_name.split("DataNewChannels")[0]  # Extract 
participant_name = participant_name.replace(" ", "").strip()

print(f"Processing data for participant: {participant_name}")

time = data['time'].values
cross940 = data['c16'].values
co940 = data['c22'].values
cross655 = data['c14'].values
co655 = data['c20'].values

fs = 50   # Sampling frequency (Hz)
cutoff = 5  # Low-pass cutoff frequency (Hz)

# ---------------------------------
# Process Each Channel
# ---------------------------------
"""
Apply a 5 Hz low-pass filter with 50 Hz sampling to remove high-frequency noise and keep cardiac content.
The Python helper uses scipy.signal.butter + filtfilt to replicate MATLAB's zero-phase lowpass()
"""

def process_channel(signal, fs, cutoff, prominence, distance, label, intervals):
    """Filter signal, find peaks/troughs, compute PI per interval, and plot (MATLAB-style)."""
    isolated = lowpass_filter(signal, cutoff, fs)

    # MATLAB equivalent to detect troughs (invert signal)
    locs_troughs, _ = find_peaks(-isolated, prominence=prominence, distance=distance)
    locs_peaks, _ = find_peaks(isolated, prominence=prominence, distance=distance)

    # Convert to values and positions
    troughs = -isolated[locs_troughs]  # amplitude values of troughs
    peaks = isolated[locs_peaks]       # amplitude values of peaks

    results = []

    for idx_list, interval_label in intervals:
        # --- ensure idx_list is a 1D array of valid integer indices ---
        idx_array = np.atleast_1d(np.array(idx_list, dtype=int).flatten())

        # convert to NumPy integer array for clean indexing
        idx_array = np.asarray(idx_array, dtype=np.int64)

        # safely slice peaks/troughs by position index (like MATLAB 14:29)
        sel_troughs = troughs[idx_array]
        sel_locs_troughs = locs_troughs[idx_array]
        sel_peaks = peaks[idx_array]

        # Compute PI for this interval (same as MATLAB)
        PI = compute_PI(sel_peaks, sel_troughs)
        results.append(PI)

        # --- Plot filtered signal with troughs and peaks ---
        plt.figure(figsize=(10, 4))
        plt.plot(time, isolated, label=f"{label} Filtered")
        plt.plot(time[sel_locs_troughs], -sel_troughs, "x", label="Troughs")
        plt.plot(time[locs_peaks[idx_array]], sel_peaks, "o", label="Peaks")

        # --- Zoom in to the region of interest ---
        # Get time range for the current interval
        t_min = time[sel_locs_troughs[0]] - 2   # small margin
        t_max = time[sel_locs_troughs[-1]] + 2
        plt.xlim([t_min, t_max])

        segment = isolated[sel_locs_troughs[0]-50 : sel_locs_troughs[-1]+50]
        plt.ylim([np.min(segment) - 500, np.max(segment) + 500])

        plt.title(f"{label}: {interval_label}")
        plt.xlabel("Time (s)")
        plt.ylabel("Filtered Signal")
        plt.legend()
        plt.tight_layout()
        #plt.show()

        # --- Define output directory ---
        output_dir = os.path.join("Generated_Plots", participant_name)
        os.makedirs(output_dir, exist_ok=True)  # create folder if not exists

        # --- Generate unique filename for each channel/interval ---
        filename = f"{label.replace(' ', '_')}_{interval_label.replace(' ', '_').replace('(', '').replace(')', '').replace('–','-')}.png"
        save_path = os.path.join(output_dir, filename)

        # --- Save only if file doesn't already exist ---
        if not os.path.exists(save_path):
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Saved plot: {save_path}")
        else:
            print(f"Skipped (already exists): {save_path}")
        
        plt.show()

    # Combine all intervals (replicates MATLAB concatenation + mean/std)
    all_PI = np.concatenate(results)
    mean_PI = np.mean(all_PI)
    std_PI = np.std(all_PI)

    print(f"{label} — Mean PI: {mean_PI:.3f}, Std PI: {std_PI:.3f}")

    return results, mean_PI, std_PI

# Define intervals (index-based like MATLAB)
"""
same sample indices to pick stable 10-s signal windows.
np.arange replicates MATLAB's start:stop, and np.concatenate replicates [a; b]
"""
intervals = [
    (np.arange(14, 30), "Interval 1 (105–115 s)"),
    (np.concatenate([np.arange(94, 96), np.arange(98, 107)]), "Interval 2 (143–153 s)"),
    (np.arange(110, 126), "Interval 3 (180–190 s)")
]


# ---------------------------------
# Run for All Channels
# ---------------------------------

cross940_results, mean_940cr, std_940cr = process_channel(cross940, fs, cutoff, 25, 10, "940 Cross", intervals)
cross655_results, mean_655cr, std_655cr = process_channel(cross655, fs, cutoff, 18, 10, "655 Cross", intervals)
co940_results, mean_940co, std_940co = process_channel(co940, fs, cutoff, 40, 8, "940 Co", intervals)
co655_results, mean_655co, std_655co = process_channel(co655, fs, cutoff, 70, 11, "655 Co", intervals)



# ---------------------------------
# Save Participant Summary (Mean + Std)
# ---------------------------------

summary_data = {
    "Participant": [participant_name],
    "940 Cross Mean": [mean_940cr], "940 Cross Std": [std_940cr],
    "655 Cross Mean": [mean_655cr], "655 Cross Std": [std_655cr],
    "940 Co Mean": [mean_940co], "940 Co Std": [std_940co],
    "655 Co Mean": [mean_655co], "655 Co Std": [std_655co],
}

summary_df = pd.DataFrame(summary_data)

# Save in the same participant folder under Generated_Plots
summary_path = os.path.join("Generated_Plots", participant_name, "summary.csv")

# create CVS for statistical analysis
if os.path.exists(summary_path):
    existing_df = pd.read_csv(summary_path)
    updated_df = pd.concat([existing_df, summary_df], ignore_index=True)
    updated_df.to_csv(summary_path, index=False)
    print(f"Updated summary CSV: {summary_path}")
else:
    summary_df.to_csv(summary_path, index=False)
    print(f"Saved new summary CSV: {summary_path}")

# ---------------------------------
# Plot Comparison (PI Across Channels)
# ---------------------------------

# Create a simple sample index axis for each concatenated PI array
pi_x = np.arange(len(np.concatenate(cross940_results)))

plt.figure(figsize=(8, 6))
plt.step(pi_x, np.concatenate(cross940_results), where='mid', label="940 Cross PI", color='r')
plt.step(pi_x, np.concatenate(cross655_results), where='mid', label="655 Cross PI", color='b')
plt.step(pi_x, np.concatenate(co940_results), where='mid', label="940 Co PI", color='m')
plt.step(pi_x, np.concatenate(co655_results), where='mid', label="655 Co PI", color='k')

plt.xlabel("Sample Index (relative PI points)")
plt.ylabel("Perfusion Index")
plt.title("Perfusion Index Across Channels (Interval 3 Example)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()


# ---------------------------------
# Signal to Noise Metic
# ---------------------------------

def pick_good_segments(signal, fs, sqi_func, win_sec=10, hop_sec=2, n_segments=3):
    """
    Scan the signal and return indices of n_segments best 10s windows by SQI.
    signal : 1D numpy array
    fs : sampling frequency (Hz)
    sqi_func : function that takes (segment, fs) -> score
    win_sec: window size
    hop_sec: how far to move the window before computing the next segments quality.
    n_segments : number of good windows to select
    """
    win = int(win_sec * fs)
    hop = int(hop_sec * fs)
    N = len(signal)
    scores = []
    
    # Compute SQI for each window
    for start in range(0, N - win + 1, hop):
        seg = signal[start:start + win]
        score = sqi_func(seg, fs)
        scores.append((start, start + win, score))
    
    # Sort windows by descending quality
    # Sorts all the windows so the best ones (highest SQI) come first.
    scores.sort(key=lambda x: x[2], reverse=True)
    
    # Goes through the ranked list from best → worst.

    selected = []
    for start, end, score in scores:
        # ensure chosen windows don't overlap with existing ones
        if all(abs(start - s[0]) > win for s in selected):
            selected.append((start, end, score))
        # Stops once 3 good segments are picked.
        if len(selected) >= n_segments:
            break
    """
    returns three separate 10-second windows spread 
    through the 2 min recording, not three overlapping ones
    """
    return selected