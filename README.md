# FIU Pulse Oximetry Analysis Pipeline

## Overview

This project processes photoplethysmography (PPG) data collected from a phantom-based pulse oximetry study. The goal is to evaluate how polarization, skin tone, vessel depth, flow speed, wavelength, and device orientation affect signal quality and perfusion index (PI).

This pipeline provides a structured and reproducible workflow for:

- loading and organizing raw experimental data  
- cleaning and extracting relevant PPG signals  
- evaluating signal quality using sliding window metrics  
- selecting the best-quality segment of each signal  
- computing perfusion index from reliable data  
- generating plots and summary tables for analysis  

---

## What This Pipeline Does (Intuition)

PPG signals are not always clean. Some parts of the signal look good and show clear pulsatile behavior, while other parts may be noisy, flat, or unstable.

Instead of trusting the entire signal, this pipeline:

1. breaks the signal into smaller windows  
2. evaluates the quality of each window  
3. selects the best window  
4. computes metrics only from that window  

This ensures that results are based on high-quality signal segments, making comparisons across experimental conditions more reliable.

---

## Channel Map

CHANNEL_MAP = {
    "Unpolarized_A": {"Green": "c5",  "Red": "c2",  "IR": "c4"},
    "Unpolarized_B": {"Green": "c11", "Red": "c8",  "IR": "c10"},
    "Co-Polarized":  {"Green": "c13", "Red": "c12", "IR": "c15"},
    "Cross-Polarized": {"Green": "c19", "Red": "c18", "IR": "c21"},
}

## Explanation of Signal Processing Pipeline

## Overview

This section explains the reasoning behind each step of the signal processing pipeline used in the FIU pulse oximetry phantom experiments.

The goal is to ensure that all computed metrics, especially perfusion index (PI), are based on high-quality portions of the signal rather than noisy or unreliable data.

---

## Why We Use Sliding Windows

PPG signals are not stable over time. Even in controlled phantom experiments, the signal can vary due to:

- flow instability  
- sensor positioning  
- noise in the system  
- drift in the signal baseline  

Because of this, analyzing the entire signal at once can lead to inaccurate results.

Instead, the signal is divided into smaller segments (windows), and each segment is evaluated independently. This allows us to identify which part of the signal is the most reliable.

---

## Why We Bandpass Filter

The bandpass filter isolates the frequency range where heart-related pulsations occur.

Typical heart rate frequencies fall between approximately:

- 0.5 Hz (30 bpm)  
- 2.2 Hz (132 bpm)  

Everything outside this range is considered noise.

The Chebyshev Type II filter is used because:

- it strongly removes unwanted frequencies  
- it provides good control over stopband attenuation  
- when applied in zero-phase mode, it does not shift the signal in time  

This step ensures that the remaining signal reflects true physiological pulsations.

---

## Why We Zero-Center the Signal

PPG signals often contain a DC offset, meaning the signal is shifted upward or downward.

Zero-centering subtracts the mean from the signal so that it is centered around zero.

This is important because:

- it makes the signal easier to analyze  
- it improves peak detection  
- it ensures consistency across different signals  

---

## Why We Detrend

Over time, PPG signals can slowly drift due to:

- temperature changes  
- sensor movement  
- electronic noise  

Detrending removes this slow baseline drift so that only the pulsatile component remains.

This helps prevent the drift from interfering with peak and trough detection.

---

## Why We Use Signal Quality Metrics

Not every window contains useful data.

To determine which window is best, we compute:

- standard deviation  
- skewness  

### Standard Deviation

This measures how much the signal varies.

- low standard deviation → flat or weak signal  
- high standard deviation → strong pulsatile variation  

### Skewness

This measures the shape of the waveform.

A good PPG waveform is not perfectly symmetric, so skewness helps identify signals that resemble real pulse shapes.

---

## How We Select the Best Window

A window is considered good if:

- standard deviation is greater than or equal to 0.02  
- skewness is greater than 0  

If no window meets both conditions, the window with the highest standard deviation is selected.

This ensures that:

- completely flat signals are avoided  
- the most informative segment is always used  

---

## Why We Invert the Signal

PPG sensors measure reflected light intensity.

When blood volume increases during a pulse:

- more light is absorbed  
- less light is detected  

This causes the signal to appear inverted relative to the actual physiological event.

Inverting the signal ensures that:

- peaks correspond to systolic pulses  
- the waveform matches expected physiological behavior  

---

## Why We Detect Peaks and Troughs

Perfusion index relies on measuring the difference between:

- peaks (maximum signal values)  
- troughs (minimum signal values)  

These correspond to:

- systolic phase (peak)  
- diastolic phase (trough)  

Accurate detection is critical for computing a reliable PI.

---

## What Perfusion Index Represents

Perfusion index measures the strength of the pulsatile signal relative to the baseline.

Higher PI means:

- stronger pulsatile signal  
- better blood flow representation  

Lower PI means:

- weaker signal  
- more noise or poor signal quality  

---

## Why We Compute PI Only on the Best Window

If PI were computed using the entire signal, it would include:

- noisy segments  
- flat regions  
- unstable portions  

This would reduce accuracy and consistency.

By computing PI only on the best window:

- results are more reliable  
- comparisons across conditions are fair  
- noise impact is minimized  

---

## Summary

The pipeline is designed to:

- isolate the cleanest part of each signal  
- remove noise and drift  
- extract meaningful physiological features  
- compute perfusion index from high-quality data  

This approach improves both the accuracy and robustness of the analysis across all experimental conditions.
