# FIU Pulse Oximetry Analysis Pipeline

# FIU Pulse Oximetry Analysis Pipeline

The objective of this study is to evaluate how physiological parameters and polarization conditions, including unpolarized, co polarized, and cross polarized light, affect the perfusion index and overall signal quality measured by our device. Using a controlled phantom setup allows us to systematically vary factors such as skin color, optical depth, and flow speed while minimizing the variability typically seen in human testing. This approach makes it possible to isolate the effect of each parameter on PPG detection and to assess the device’s limitations under different optical and physical conditions. By analyzing performance across skin tones and depths, we aim to identify which polarization state provides the most reliable and consistent PPG measurements. Together, these investigations will validate the device’s detection capabilities and guide improvements for future human testing.

---

## Overview

This pipeline is designed to:

- load and organize raw experimental data
- clean and extract relevant PPG signals
- evaluate signal quality using sliding window metrics
- select the best-quality segment of each signal
- compute perfusion index from reliable data
- generate plots and summary tables for downstream analysis

Rather than trusting the entire signal equally, the pipeline identifies the most reliable segment of each recording and computes metrics from that portion only. This improves consistency across experimental conditions and reduces the effect of noise, drift, or unstable signal segments.

---

## Study Design

The phantom study was built around a factorial design that varied key optical and physiological parameters:

- 3 skin tones
- 3 phantom flow speeds
- 3 vessel depths
- 3 polarization conditions

This gives a total of 81 theoretical conditions. However, because the device records unpolarized, co-polarized, and cross-polarized conditions simultaneously, the effective number of testing samples is reduced to 27. 

### Experiment 1 Variables
- skin tones: fair, medium, dark
- flow speeds: slow (60 BPM), intermediate (90 BPM), fast (120 BPM)
- vessel depths: 2.5 mm, 3.5 mm, 5 mm
- polarization conditions recorded simultaneously

### Experiment 2 Variables
- phantom wavelengths: 525 nm, 660 nm, 940 nm
- fixed vessel depth: 3.5 mm
- fixed speed: intermediate, 90 BPM
- multilayered phantom
- polarization conditions recorded and compared

### Experiment 3 Variables
- 940 nm condition
- pulse oximeter orientation rotated to 0 degrees, 90 degrees, and 180 degrees
- original polarization and unpolarized conditions compared under rotation testing

---

## Limitations

We recognize that the phantom cell model cannot fully replicate the complexity of human anatomy or physiological variability. However, demonstrating that the device can generate and detect PPG signals within this controlled system provides valuable insight into its fundamental capabilities and potential limitations. The results from this model will help identify areas where the device performs reliably and where further refinement may be needed before proceeding to human testing.

---

## Experimental Workflow

### Day 1: Calibration
The first day focused on calibrating the device to the phantom and familiarizing the Ramella Lab with the system. During calibration, it was observed that when the polarizers were placed on the top row, the channels were often very close to the noise floor, making the synthetic PPG signal more difficult to detect. Data were collected for fair and medium skin tone phantoms, and some files were repeated. 

### Day 2: Experiment 1
The second day focused on structured data collection across skin tone, vessel depth, and flow speed conditions. During this phase, the team also began investigating the effect of polarization by comparing signal quality with and without polarizers. This was motivated by the observation that channel 18 consistently produced a stronger signal even when the polarization was flipped. To improve repeatability, the Ramella Lab printed custom clamps to better regulate contact pressure during measurements. 

### Day 3: Experiment 2
The third day focused on wavelength-based experiments using multilayered phantoms with a fixed vessel depth of 3.5 mm and an intermediate speed of 90 BPM. Measurements were collected without clamps at first, and although consistent data could still be obtained, it became clear that clamps were important for maintaining consistency across experiments. Additional testing was also performed to investigate whether an offset was present in the data by turning off the photodiodes and LEDs to see whether the channels would properly zero out. Polarization was switched again to investigate whether signal behavior was tied to polarization condition or to channel-specific dependence. 

### Day 4: Experiment 2 and Experiment 3
The fourth day continued wavelength testing under more controlled conditions using clamps. These experiments again used multilayered phantoms, a vessel depth of 3.5 mm, and an intermediate speed of 90 BPM. Additional trials were performed to test both flipped polarization and device rotation, especially for the 940 nm channel. Rotation testing was conducted at 0 degrees, 90 degrees, and 180 degrees for fair and dark skin tone conditions.
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
