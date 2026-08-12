# FIU Pulse Oximetry Analysis Pipeline

The goal of this project is to evaluate how different experimental conditions affect the quality of the PPG signals measured by our pulse oximetry device. In particular, we are interested in understanding how factors such as skin tone, vessel depth, flow speed, wavelength, and polarization condition affect signal quality and perfusion index (PI).

Using a controlled phantom setup allows us to change these conditions one at a time while reducing some of the variability that would normally come with human testing. This gives us a better way to understand how the device behaves under different optical and physical conditions before moving toward larger human studies.

A major focus of the analysis is comparing unpolarized, co-polarized, and cross-polarized measurements across the different phantom conditions. By looking at signal quality at both the window and individual beat level, we can better understand which conditions consistently produce reliable PPG signals and where the device may have limitations.

---

## Overview

The current pipeline is designed to:

- load and organize the FIU phantom recordings
- process only uncompressed `.json` recordings to avoid analyzing duplicate `.json.gz` files
- extract the PPG channels of interest
- filter and clean the PPG signals
- divide each recording into smaller windows
- detect individual pulse beats within each window
- create a representative beat template
- compare individual beats to that template
- calculate multiple signal quality metrics for each beat
- classify beats based on their signal quality
- determine the quality of each window based on its beats
- calculate perfusion index (PI) on a beat-by-beat basis
- save beat-level, window-level, and recording-level results to CSV files

Instead of assuming that an entire recording has the same signal quality, the pipeline looks more closely at what is happening throughout the recording. A recording may contain sections with very clean PPG waveforms and other sections that are noisy or unstable. Looking at smaller windows and individual beats gives us a much more detailed picture of the signal.

---

# Study Design

The phantom study was designed around several optical and physiological parameters:

- 3 skin tones
- 3 phantom flow speeds
- 3 vessel depths
- 3 polarization conditions

This gives a total of 81 theoretical conditions. However, because the device records unpolarized, co-polarized, and cross-polarized conditions simultaneously, the effective number of testing samples is reduced to 27.

## Experiment 1 Variables

Experiment 1 focused on changing:

- skin tone: fair, medium, dark
- flow speed: slow (60 BPM), intermediate (90 BPM), fast (120 BPM)
- vessel depth: 2.5 mm, 3.5 mm, 5 mm
- polarization condition

These recordings allow us to look at how changes in the phantom itself affect the quality of the PPG signal.

## Experiment 2 Variables

Experiment 2 focused more specifically on wavelength.

The main conditions included:

- phantom wavelengths: 525 nm, 660 nm, and 940 nm
- fixed vessel depth: 3.5 mm
- fixed flow speed: intermediate (90 BPM)
- multilayered phantom
- different polarization conditions

## Experiment 3 Variables

Experiment 3 focused on the effect of device orientation.

The main conditions included:

- 940 nm wavelength
- pulse oximeter orientations of 0°, 90°, and 180°
- original polarization and unpolarized conditions
- comparisons across rotation conditions

---

# Experimental Workflow

## Day 1: Calibration

The first day mainly focused on calibrating the device to the phantom setup and becoming familiar with how the system behaved.

During calibration, we noticed that when the polarizers were placed on the top row, some of the channels were very close to the noise floor. This made the synthetic PPG signal harder to detect.

Data were collected for fair and medium skin tone phantoms, and some recordings were repeated while the setup was being adjusted.

Because these recordings were mainly used for calibration and troubleshooting, they are not currently part of the main beat-level SQI analysis.

---

## Day 2: Experiment 1

Day 2 focused on Experiment 1 and included structured data collection across:

- skin tone
- vessel depth
- flow speed
- polarization condition

During this experiment, we also started looking more closely at the effect of polarization. This was partly motivated by seeing that certain channels consistently produced stronger signals even when the polarization configuration was changed.

Custom clamps were also introduced to help regulate contact pressure and improve consistency between measurements.

The heartbeat recordings from Day 2 are processed using the beat-level SQI pipeline.

Results are saved under:

```text
FIU_Beat_Level_SQI/Day_2/Experiment_1/
```

---

## Day 3: Experiment 2

Day 3 focused on wavelength testing using multilayered phantoms.

The main heartbeat recordings used:

- 3.5 mm vessel depth
- intermediate speed (90 BPM)
- no clamps
- original polarization configuration

Although we were able to collect consistent data without the clamps, it became clear that having controlled contact pressure would be important for making comparisons across experiments more consistent.

Additional calibration testing was also performed during Day 3. This included blackout/offset testing and switching the polarization configuration to investigate whether differences in the signals were related to the polarization itself or to specific channels.

For the beat-level SQI analysis, we only process the **multilayered heartbeat recordings** from Day 3.

Results are saved under:

```text
FIU_Beat_Level_SQI/Day_3/Experiment_2/
```

---

## Why Day 3 Calibration Data Are Not Included in the Beat-Level Analysis

The calibration recordings from Day 3 were collected for a different purpose than the heartbeat recordings.

These recordings were used to:

- check whether the channels properly zero out
- investigate possible offsets
- test blackout conditions
- investigate channel dependence
- test changes in polarization placement

The beat-level SQI pipeline assumes that the signal contains repeated pulse cycles that can be segmented and compared.

Because the calibration recordings were not necessarily collected to represent normal periodic heartbeat signals, applying the same heartbeat analysis to them would not always be meaningful.

For that reason, the current Day 3 beat-level analysis focuses specifically on the multilayered heartbeat recordings.

---

## Day 4: Experiments 2 and 3

Day 4 continued Experiment 2 wavelength testing under more controlled conditions and also included Experiment 3 orientation testing.

These recordings again used multilayered phantoms with:

- 3.5 mm vessel depth
- intermediate speed (90 BPM)
- controlled contact using clamps

Additional trials looked at flipped polarization and device rotation, particularly for the 940 nm condition.

Rotation testing was performed at:

- 0°
- 90°
- 180°

The pipeline separates Experiment 2 and Experiment 3 so that the results can be analyzed independently.

Experiment-specific results are saved under:

```text
FIU_Beat_Level_SQI/Day_4/Experiment_2/
```

and:

```text
FIU_Beat_Level_SQI/Day_4/Experiment_3/
```

Combined Day 4 results are also saved directly inside the `Day_4` folder.

---

# Channel Map

The following channel map is used to connect each polarization and wavelength condition to its corresponding device channel.

```python
CHANNEL_MAP = {
    "Unpolarized_A": {
        "Green": "c5",
        "Red": "c2",
        "IR": "c4"
    },
    "Unpolarized_B": {
        "Green": "c11",
        "Red": "c8",
        "IR": "c10"
    },
    "Co-Polarized": {
        "Green": "c13",
        "Red": "c12",
        "IR": "c15"
    },
    "Cross-Polarized": {
        "Green": "c19",
        "Red": "c18",
        "IR": "c21"
    },
}
```

---

# Signal Processing Pipeline

The goal of the signal processing pipeline is to take the raw PPG recordings and determine how reliable the signal is at a much smaller scale.

The overall structure of the analysis is:

```text
Recording
    ↓
Windows
    ↓
Individual Beats
    ↓
Beat-Level SQI Metrics + Perfusion Index
```

This lets us look at signal quality at multiple levels instead of assigning one value to an entire recording.

---

## 1. Signal Preprocessing

Before detecting beats, the PPG signal is cleaned and filtered.

PPG recordings can contain:

- baseline drift
- electronic noise
- high-frequency noise
- DC offsets
- slow changes that are unrelated to the pulse

Filtering helps isolate the part of the signal that contains the periodic pulse waveform.

---

## 2. Bandpass Filtering

The bandpass filter keeps the frequency range where we expect the heartbeat signal to occur while removing frequencies outside that range.

The expected heart rate range is approximately:

- 0.5 Hz = 30 BPM
- 2.2 Hz = 132 BPM

This range includes the phantom speeds used in the experiments:

- 60 BPM
- 90 BPM
- 120 BPM

Removing frequencies outside this range helps make the pulse waveform easier to detect.

---

## 3. Zero-Centering and Detrending

PPG signals can have a baseline offset or slowly drift over time.

Zero-centering removes the mean from the signal so that the waveform is centered around zero.

Detrending removes slow changes in the baseline that could interfere with pulse detection.

Together, these steps make the pulsatile part of the signal easier to analyze.

---

# Window-Level Analysis

## Why We Use Windows

Signal quality is not always consistent throughout an entire recording.

For example, one part of a recording may have a very clear pulse waveform while another part may contain noise or instability.

Instead of treating the entire recording as equally reliable, the pipeline divides it into smaller windows.

Each window can then be evaluated separately.

This makes it easier to identify:

- clean portions of the recording
- noisy portions
- unstable pulse signals
- changes in signal quality over time

---

# Beat-Level Analysis

## Beat Segmentation

After the signal has been divided into windows, individual pulse beats are detected.

The beats are segmented from valley to valley so that each segment represents approximately one complete pulse cycle.

This gives us individual waveforms that can be compared instead of only looking at the average behavior of the entire window.

---

## Representative Beat Template

Once the individual beats in a window are detected, the beats are normalized and used to create a representative beat template.

The template gives us an estimate of what a typical beat looks like within that section of the recording.

Each detected beat can then be compared against this template.

A beat that looks very similar to the template is more likely to represent a consistent PPG pulse.

A beat that looks very different may contain noise, distortion, or an unstable waveform.

---

# Beat-Level Signal Quality Metrics

Rather than relying on only one measurement, the pipeline uses multiple SQI metrics to describe different parts of beat quality.

The main metrics currently include:

- Dynamic Time Warping (DTW)
- correlation
- Mean Absolute Deviation (MAD)
- template-based signal quality information
- beat acceptance
- perfusion index (PI)

---

## Dynamic Time Warping (DTW)

Dynamic Time Warping measures how different the shape of an individual beat is from the representative beat template.

It allows for small differences in timing while still comparing the overall morphology of the waveforms.

In general:

- smaller DTW distance = beat is more similar to the template
- larger DTW distance = beat differs more from the template

This makes DTW useful for identifying beats with unusual or distorted shapes.

---

## Correlation

Correlation measures how closely the shape of an individual beat follows the representative template.

In general:

- higher correlation = stronger similarity
- lower correlation = weaker similarity

A high correlation suggests that the beat follows the general waveform shape expected for that window.

---

## Mean Absolute Deviation (MAD)

MAD measures the average difference between an individual beat and the representative template.

In general:

- lower MAD = beat is closer to the template
- higher MAD = beat differs more from the template

While correlation focuses more on whether two waveforms follow a similar shape, MAD gives us information about how far apart they are.

---

# How the SQI Metrics Work Together

The current pipeline uses a **rule-based approach** rather than a machine-learning classifier.

This means that the SQI metrics are calculated separately and evaluated using predefined thresholds.

The pipeline does not currently train a model to combine the SQIs into one learned signal quality score.

Using individual thresholds gives us a more interpretable starting point because we can see exactly which signal-quality criteria a beat does or does not meet.

The current thresholds should still be considered preliminary. As more of the FIU data are analyzed, we can look at the distributions of these metrics and determine whether the thresholds should be adjusted.

A future version of the pipeline could use these SQIs as features for a machine-learning model that learns how to classify good and bad beats.

---

# Beat-Level Perfusion Index

Perfusion index has now been added as an additional **beat-level metric**.

Instead of only calculating PI from a larger portion of the recording, the updated pipeline calculates PI for each individual detected beat.

Perfusion index represents the strength of the pulsatile part of the signal relative to the baseline signal.

Conceptually:

```text
PI = (AC / DC) × 100
```

where:

- **AC** represents the pulsatile change in the PPG signal
- **DC** represents the underlying baseline intensity

A larger PI generally represents a stronger pulsatile component relative to the baseline.

Calculating PI beat-by-beat gives us more information than calculating only one PI value for an entire recording.

For example, within the same recording we can now investigate whether:

- PI stays relatively consistent between beats
- PI changes when signal quality decreases
- different polarization conditions produce different PI distributions
- PI changes across skin tone, wavelength, depth, or flow speed

PI is currently included as an **additional measurement** rather than being used by itself to determine whether a beat is good or bad.

This allows us to study how PI relates to the other SQI measurements before deciding whether it should eventually contribute to beat classification.

---

# Output CSV Files

The pipeline saves the analysis at three main levels:

```text
Beat Level
    ↓
Window Level
    ↓
Recording Level
```

Each CSV serves a different purpose.

---

## 1. Beat-Level SQI CSV

Example:

```text
day2_experiment1_all_beat_sqi.csv
```

This is the most detailed output produced by the pipeline.

Each row represents an individual detected beat.

The beat-level CSV contains measurements that allow us to look at the quality of individual pulse cycles, including the calculated SQI metrics and beat-level perfusion index.

This file can be used to:

- compare individual beats
- look at distributions of SQI metrics
- compare PI across beats
- investigate why certain beats fail
- compare signal quality across experimental conditions
- eventually create features for machine-learning analysis

---

## 2. Window-Level SQI CSV

Example:

```text
day2_experiment1_all_window_sqi.csv
```

This file summarizes the signal quality of each analyzed window.

Instead of looking at every individual beat, it gives us a higher-level view of how well each section of the recording performed.

This makes it useful for seeing:

- how many beats were detected
- how many beats passed the SQI criteria
- how signal quality changes throughout a recording
- whether certain sections of a recording are consistently better than others

---

## 3. Recording Summary CSV

Example:

```text
day2_experiment1_recording_summary.csv
```

This file gives us the highest-level summary of the analysis.

It connects the SQI results back to the experimental conditions.

Depending on the experiment, this includes information such as:

- day
- experiment
- source recording
- skin tone
- vessel depth
- flow speed
- expected BPM
- wavelength
- polarization condition
- device orientation
- hardware channel
- overall signal quality information

This file is especially useful when comparing results across experimental conditions.

---

# Current Output Folder Structure

The processed results are organized by day and experiment.

```text
FIU_Beat_Level_SQI/
│
├── Day_2/
│   └── Experiment_1/
│       ├── day2_experiment1_all_beat_sqi.csv
│       ├── day2_experiment1_all_window_sqi.csv
│       └── day2_experiment1_recording_summary.csv
│
├── Day_3/
│   └── Experiment_2/
│       ├── day3_experiment2_all_beat_sqi.csv
│       ├── day3_experiment2_all_window_sqi.csv
│       └── day3_experiment2_recording_summary.csv
│
└── Day_4/
    │
    ├── Experiment_2/
    │   ├── day4_experiment_2_beat_sqi.csv
    │   ├── day4_experiment_2_window_sqi.csv
    │   └── experiment_2_summary.csv
    │
    ├── Experiment_3/
    │   ├── day4_experiment_3_beat_sqi.csv
    │   ├── day4_experiment_3_window_sqi.csv
    │   └── experiment_3_summary.csv
    │
    ├── day4_all_beat_sqi.csv
    ├── day4_all_window_sqi.csv
    └── day4_recording_summary.csv
```

Day 4 includes both experiment-specific files and combined Day 4 files so that we can either analyze each experiment separately or look at all Day 4 recordings together.

---

# Avoiding Duplicate Recordings

Some of the original FIU folders contain both:

```text
recording.json
```

and:

```text
recording.json.gz
```

These can represent the same recording in compressed and uncompressed formats.

To avoid accidentally analyzing the same recording twice, the current processing pipeline only searches for and processes the uncompressed `.json` recordings.

The `.json.gz` files are therefore not included in the heartbeat analysis.

---

# Current Analysis Goal

The current version of the pipeline gives us a way to evaluate PPG quality at several different levels:

```text
Experimental Condition
        ↓
Recording
        ↓
Window
        ↓
Individual Beat
        ↓
SQI Metrics + Perfusion Index
```

This lets us move beyond simply asking whether an entire recording looks good or bad.

Instead, we can investigate questions such as:

- Are certain polarization conditions producing more consistent beats?
- Does signal quality change across skin tones?
- Does vessel depth affect beat morphology?
- Does flow speed affect beat detection?
- Do certain wavelengths produce stronger or more consistent signals?
- Does device orientation affect signal quality?
- How does perfusion index change across these conditions?
- Is PI related to whether a beat passes or fails the SQI criteria?

The goal is to use these measurements to better understand where the device performs consistently and where signal quality begins to break down.

---

# Limitations

The phantom model cannot fully recreate the complexity of human anatomy or physiological variability.

However, the controlled phantom setup gives us an important way to test the device under repeatable conditions.

It allows us to isolate specific variables and determine whether the device can consistently detect PPG signals before introducing the additional variability that comes with human testing.

The current SQI thresholds are also preliminary.

They provide a starting point for separating more consistent beats from lower-quality beats, but they have not yet been treated as final validated thresholds.

As more of the dataset is analyzed, these values can be revisited using the actual distributions of the SQI metrics.

---

# Future Work

The current beat-level pipeline gives us a foundation for more detailed signal-quality analysis.

Possible next steps include:

- analyze the distributions of each SQI metric across the full dataset
- refine the current SQI thresholds using the experimental data
- compare SQI distributions between skin tones
- compare signal quality across polarization conditions
- compare signal quality across wavelengths
- compare PI between good and bad beats
- determine whether PI should contribute to beat-quality classification
- evaluate whether certain SQIs are more informative than others
- combine multiple SQIs using a machine-learning model
- investigate which polarization and wavelength combinations consistently provide the strongest PPG signals
- use the phantom results to guide future human testing

Ultimately, the goal is to build a more complete understanding of how the device performs across different optical and physical conditions and use that information to improve the reliability of future PPG measurements.