> Background reading for this project -- study design, experiment variables, and a full walkthrough of how and why the pipeline works the way it does.
>
> If you just want to set this up and run it yourself, see [README.md](README.md) instead.

---

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
- divide each recording into 10-second windows
- filter each window and detect individual pulse beats within the filtered window
- calculate a perfusion index (PI) for each beat from the raw, unfiltered signal, and summarize it per window (median across beats)
- normalize each beat individually and build a representative beat template from the normalized beats
- compare individual beats to that template
- calculate multiple signal quality metrics for each beat
- classify beats based on their signal quality, then roll those beat labels up into a good/bad label for the window
- save window-level and recording-level results to CSV files

Instead of assuming that an entire recording has the same signal quality, the pipeline looks more closely at what is happening throughout the recording. A recording may contain sections with very clean PPG waveforms and other sections that are noisy or unstable. Looking at smaller windows and individual beats gives us a much more detailed picture of the signal.

---

# Recent Pipeline Updates

**August 2026:** Following a review of the processing order, a few changes were made to `beat_level_sqi.py`:

- The bandpass filter's upper cutoff was widened from 2.2 Hz (132 BPM) to 4 Hz (240 BPM).
- Beat quality metrics (DTW distance, correlation, MAD, template SQI) are now calculated **before** beats are rescaled in height, using beats that have only been lined up in time. Previously, beats were rescaled in height first, which could hide real quality differences tied to a beat's amplitude.
- The height-rescaling step itself changed from a mean/standard-deviation rescaling to a minimum/maximum rescaling (0 to 1), and now runs *after* quality has been determined rather than before, since it's used for downstream purposes rather than the quality decision itself.
- Because quality metrics are now computed in real amplitude units instead of rescaled ones, the current SQI thresholds are expected to need re-tuning (see Limitations).

**September 2026:** After more feedback from Mark and Rutendo, the pipeline went through another round of changes -- a bigger one this time.

We tried a window-level-only version for a bit: filter the window, calculate PI, normalize the whole window at once, and classify it directly with no beat segmentation. Running it on real data came back almost 100% "good," which traced back to two problems -- the window-level check was only looking at clipping (a pretty lenient check on its own), and PI was blowing up to hundreds/thousands of percent because it was being calculated on the filtered signal.

Mark and Rutendo clarified that we still needed real beat-to-beat comparison, and that normalization was supposed to happen per beat, not per window. So the pipeline moved back toward beat-level comparison -- closer to how it worked before the August update above, but with a few real changes:

- Beats are now normalized (resampled to a fixed length, then min-max scaled to 0-1) **before** being compared to the template, not after. This reverses the August change above -- the DTW distance, correlation, MAD, and template SQI described later in this README are now calculated on the normalized beats, not the real-amplitude ones. Clipping SQI is the one exception, since it's a ratio and gives the same answer either way.
- Perfusion index moved from a beat-level metric to a window-level one, and is now calculated from the **raw, unfiltered** window signal instead of the filtered one. Filtering removes the DC baseline that PI needs, which is what was causing PI to blow up in the window-level-only attempt above.
- Beat quality is still decided beat-by-beat, but the beat-level CSV is no longer saved on its own -- each window's beat results are now summarized directly into the window-level CSV (mean DTW distance, mean correlation, mean MAD, percent of bad beats, etc.) instead of writing out one row per beat.
- The SQI thresholds (`sqi_lambda`, `min_template_sqi`, `max_mad`) were briefly retuned using real per-beat percentile data, since they were calibrated for the old real-amplitude scale and barely triggered on the new normalized one. The retuned values barely changed the overall good/bad split, so we reverted back to the original thresholds and documented the reasoning in the code instead of carrying non-standard values forward.

**September 2026 (Rutendo's review):** Rutendo went through the code against the two reference papers (Li & Clifford 2012; Zia et al. 2020) and raised a few issues, which led to these changes:

- **PI is now calculated per beat, not per window.** Taking the max and min of a whole 10-second window lets baseline drift or an artifact inflate "AC". PI is now calculated for each beat on the raw signal (using the same beat boundaries found on the filtered signal), and the window's PI is the median across its beats. The mean and SD of the beat PIs are saved too.
- **The DTW-based SQI now follows Zia et al. exactly:** SQI = exp(-λ · D/L), where L is the length of the warping path. Before, we used `dtw-python`'s `normalizedDistance`, which divides by the two beats' lengths added together (N + M) rather than L, and uses a step pattern that double-counts diagonal steps.
- **λ is back to 25 (Zia et al.'s value).** 25 was the pipeline's original value too; it was changed to 0.30 on Sept 1 when the metrics briefly ran on real amplitudes, and it was never changed back once beats were normalized again. So 0.30 was never tuned on FIU data.
- **AD was removed.** Every beat is resampled to 50 points, so AD = 50 × MAD exactly and it added nothing.
- **MAD is no longer used to reject beats.** It isn't an SQI in either reference paper (it was something we added), and the old `max_mad = 30` could never be reached on 0-to-1 normalized beats. MAD is still saved as a descriptive column.

**September 2026 (results review):** Going through the saved results turned up a few more issues:

- **Beat segmentation fix.** In the Slow (60 BPM label) recordings, each pulse has a deep trough followed by a shallow secondary notch, and the valley detector was counting both, splitting every beat in two. Slow windows were only 1.2% good vs. ~25% at the other speeds. Valleys now have to be at least 60% of the window's own pulse period apart (period taken from the window's dominant frequency), on top of the old max-HR spacing. Slow went from 1.2% to 30.6% good; Fast and Intermediate barely changed (25.2 → 25.4%, 24.7 → 24.9%).
- **Condition labels normalized.** "Light" is now mapped to "Fair" for skin tone, and the trailing period in "Og. Pol." is dropped, so conditions group correctly across folders.
- **Day 4 combined files no longer double-count.** The IR 0° Og. Pol recordings are used in both Experiment 2 and Experiment 3. They still appear in each experiment's own files, but only once in the combined Day 4 files, labeled "Experiment_2 & Experiment_3".
- **Old beat-level CSVs archived.** The Aug 31 `*_beat_sqi.csv` files (old pipeline) were moved to `FIU_Beat_Level_SQI/_archive_old_beat_sqi_aug31/` so they aren't mistaken for current results.

Still open (needs Mark/Rutendo):

- **Clipping SQI doesn't measure clipping at 25 Hz.** For a typical beat it works out to exactly 1 − 2/N (N = samples in the beat), because it's only counting the single max and min sample. So the 0.80 cutoff just rejects beats shorter than 10 samples, and it fails ~21% of beats that correlate ≥ 0.95 with the template. In a 22-file sample, removing it takes good windows from ~36% to ~61%. Options: replace it with a real saturation check, or drop it.
- **Actual pump rates don't match the labels.** The spectra show ~72 / 105 / 130 BPM for Slow / Intermediate / Fast, not 60 / 90 / 120 (the 25 Hz sampling rate checks out against the timestamps).
- **Time gaps inside recordings.** About half of the Day 2 files have a jump in their timestamps (often > 30 s), and the loader joins the pieces end to end, so windows crossing the join mix two separate stretches of recording.

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
Individual Beats (segmented from the filtered window)
    ↓
Beat-Level SQI Metrics  +  Beat-Level PI (from the raw beat)
    ↓
Window-Level Good/Bad Label
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
- 4.0 Hz = 240 BPM

This range comfortably covers the phantom speeds used in the experiments:

- 60 BPM
- 90 BPM
- 120 BPM

The upper cutoff was widened from 2.2 Hz (132 BPM) to 4 Hz (240 BPM) so the filter keeps more of the natural shape of each pulse instead of smoothing it away. The narrower cutoff was closer to only the fundamental heart-rate frequency, which could flatten out real shape differences between a clean beat and a distorted one before the beat-level SQI metrics get a chance to see them.

Removing frequencies outside this range still helps make the pulse waveform easier to detect.

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

To keep a pulse's secondary notch from being counted as its own valley, valleys must be at least 60% of the window's pulse period apart. The period comes from the window's dominant frequency.

This gives us individual waveforms that can be compared instead of only looking at the average behavior of the entire window.

---

## Normalizing Beats Before Comparison

Beats don't all last the same amount of time -- a fast heartbeat produces a shorter beat than a slow one -- and they don't all sit at the same height either, depending on things like contact pressure or channel. Before beats can be compared to each other or to a template, each one is resampled onto the same fixed-length timeline and then min-max scaled so its values run from 0 to 1.

This is a change from how the pipeline worked earlier: beats used to be lined up in time only, keeping their real height, and height differences were exactly what the quality metrics were meant to catch. After talking it through with Mark and Rutendo, this per-beat normalization turned out to be what was actually intended -- so both the timing *and* the height of each beat are standardized before anything gets compared.

---

## Representative Beat Template

Once the individual beats in a window have been normalized (as described above), they're used to build a representative beat template -- built from the normalized versions, not the real-amplitude ones.

The template is built by averaging the beats that already look most alike, so a few odd or noisy beats can't drag it off course. The result is an estimate of what a typical *normalized* beat looks like within that section of the recording.

Each detected beat can then be compared against this template.

A beat that looks very similar to the template is more likely to represent a consistent PPG pulse. A beat that looks very different may contain noise, distortion, or an unstable waveform. Clipping is checked separately, straight off the beat's real (non-normalized) amplitude -- see Clipping SQI below.

---

# Beat-Level Signal Quality Metrics

Rather than relying on only one measurement, the pipeline uses multiple SQI metrics to describe different parts of beat quality. Except for clipping, these are all calculated on the *normalized* version of the beat (see "Normalizing Beats Before Comparison" above) -- so they're comparing shape and relative amplitude within a standardized 0-to-1 range, not raw signal height.

The metrics used to accept or reject a beat are:

- DTW-based template SQI (Zia et al. 2020)
- correlation after linear resampling (Li & Clifford 2012)
- clipping (our own check, loosely based on Li & Clifford's clipping detection)

MAD and skewness are also calculated and saved, but only as descriptive values -- they are not used to reject beats.

Perfusion index (PI) is calculated for each beat too, from the raw signal -- see "Perfusion Index" below.

**What we took from each paper vs. what we changed.** Li & Clifford calculate three correlation-based SQIs (direct matching, correlation after linear resampling, and correlation after DTW alignment) plus a clipping check. We only use the linear-resampling correlation. For DTW we use Zia et al.'s distance-based SQI instead of Li & Clifford's DTW correlation. Our clipping formula is our own.

---

## Dynamic Time Warping (DTW)

Dynamic Time Warping measures how different an individual beat's shape is from the representative beat template, using the normalized versions of both.

Following Zia et al. (2020), the DTW distance D is divided by the length L of the warping path, and turned into a template SQI with SQI = exp(-λ · D/L), with λ = 25.

It allows for small differences in timing while still comparing the overall morphology of the waveforms.

In general:

- smaller DTW distance = beat is more similar to the template
- larger DTW distance = beat differs more from the template

This makes DTW useful for identifying beats with unusual or distorted shapes, even after amplitude differences between beats have been normalized away.

---

## Correlation

Correlation measures how closely the shape of an individual beat follows the representative template. Correlation isn't affected by how tall or short a beat is -- it only cares about shape -- so it gives the same answer whether it's computed on the raw beat or the normalized one.

In general:

- higher correlation = stronger similarity
- lower correlation = weaker similarity

A high correlation suggests that the beat follows the general waveform shape expected for that window.

---

## Mean Absolute Deviation (MAD) -- descriptive only

MAD isn't from either reference paper, so it is no longer used to reject beats; it's kept as a descriptive column. (AD, the summed version, was removed because AD = 50 × MAD exactly once beats are resampled to 50 points.)

MAD measures the average difference between an individual beat and the representative template, sample by sample, using the normalized versions of both. Unlike correlation, MAD is still affected by height -- it's just the *normalized* height now, so it's picking up on shape differences that survive normalization rather than raw amplitude differences.

In general:

- lower MAD = beat is closer to the template
- higher MAD = beat differs more from the template

While correlation focuses more on whether two waveforms follow a similar shape, MAD gives us information about how far apart they are, sample by sample.

---

## Clipping SQI

Clipping SQI estimates how much of a beat is sitting at (or very near) its own minimum or maximum value, which is usually a sign of signal clipping or saturation.

This one is calculated on the beat's real, non-normalized amplitude -- but since it's a ratio relative to the beat's own min and max, it gives the same result whether or not the beat has been normalized, so it didn't need to change with the rest of the pipeline.

In general:

- closer to 1.0 = little or no clipping
- lower values = more of the beat is clipped

---

# How the SQI Metrics Work Together

The current pipeline uses a **rule-based approach** rather than a machine-learning classifier.

This means that the SQI metrics are calculated separately and evaluated using predefined thresholds.

The pipeline does not currently train a model to combine the SQIs into one learned signal quality score.

Using individual thresholds gives us a more interpretable starting point because we can see exactly which signal-quality criteria a beat does or does not meet.

The current thresholds should still be considered preliminary. As more of the FIU data are analyzed, we can look at the distributions of these metrics and determine whether the thresholds should be adjusted.

A future version of the pipeline could use these SQIs as features for a machine-learning model that learns how to classify good and bad beats.

---

# Perfusion Index

Perfusion index is calculated for each **beat**, from the raw, unfiltered signal, using the same beat boundaries found on the filtered signal. Each window's PI is then the **median** of its beat PIs (the mean and SD of the beat PIs are saved as `mean_beat_pi` and `std_beat_pi`).

An earlier version calculated one PI per window from the max and min of the whole 10-second window. That lets baseline drift or a single artifact inflate "AC", so it was changed after Rutendo's review.

Perfusion index represents the strength of the pulsatile part of the signal relative to the baseline signal.

Conceptually:

```text
PI = (AC / DC) × 100
```

where:

- **AC** represents the pulsatile change in the PPG signal
- **DC** represents the underlying baseline intensity

A larger PI generally represents a stronger pulsatile component relative to the baseline.

PI is calculated on the raw signal specifically because our bandpass filter removes the DC baseline on purpose -- if PI were calculated on the filtered signal, DC would be artificially close to zero, and dividing by a near-zero number is what was sending PI up to hundreds or thousands of percent in an earlier version of the pipeline (see "Recent Pipeline Updates" above).

For example, within the same recording we can now investigate whether:

- PI stays relatively consistent between windows
- PI changes when window quality decreases
- different polarization conditions produce different PI distributions
- PI changes across skin tone, wavelength, depth, or flow speed

PI is currently included as an **additional measurement** rather than being used by itself to determine whether a window is good or bad.

This allows us to study how PI relates to the other SQI measurements before deciding whether it should eventually contribute to window classification.

One open item: a small number of recordings still show extreme PI values, and they all trace back to one specific channel (Unpolarized_A_Green / hardware channel C5), where the raw signal itself has near-zero amplitude. Switching to per-beat PI shrank this a lot (recording-channel combos with mean PI > 100% went from 116 to 56, and the worst case from ~29,000% to ~480%), and all 56 remaining are on C5. This looks more like a hardware/sensor issue on that channel than a problem with the PI calculation itself, but it's still an open question for Mark and Rutendo (see Limitations).

---

# Output CSV Files

The pipeline saves the analysis at two main levels:

```text
Window Level
    ↓
Recording Level
```

Each window's beat-level results (DTW distance, correlation, MAD, clipping, percent of bad beats, etc.) are summarized directly into the window-level CSV rather than being saved as a separate beat-by-beat file -- see "Recent Pipeline Updates" above.

Each CSV serves a different purpose.

---

## 1. Window-Level SQI CSV

Example:

```text
day2_experiment1_all_window_sqi.csv
```

This is the most detailed output the pipeline currently saves. Each row represents one 10-second window and includes the window's perfusion index (median of its beat PIs, plus their mean and SD), its beat-level results summarized into means (mean DTW distance, mean correlation, mean MAD, mean template SQI, mean clipping SQI, mean skewness), the number of good/bad beats, the percent of bad beats, and the resulting good/bad window label.

This makes it useful for seeing:

- how many beats were detected in each window, and how many passed the SQI criteria
- how signal quality changes throughout a recording, window by window
- whether certain sections of a recording are consistently better than others
- how a window's PI relates to its other SQI values

---

## 2. Recording Summary CSV

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
- number of good/bad windows and the percent good/bad
- mean PI, mean DTW distance, mean correlation, mean MAD, mean template SQI, mean clipping SQI, and mean skewness across the recording's windows

This file is especially useful when comparing results across experimental conditions.

---

# Current Output Folder Structure

The processed results are organized by day and experiment.

```text
FIU_Beat_Level_SQI/
│
├── Day_2/
│   └── Experiment_1/
│       ├── day2_experiment1_all_window_sqi.csv
│       └── day2_experiment1_recording_summary.csv
│
├── Day_3/
│   └── Experiment_2/
│       ├── day3_experiment2_all_window_sqi.csv
│       └── day3_experiment2_recording_summary.csv
│
└── Day_4/
    │
    ├── Experiment_2/
    │   ├── day4_experiment_2_window_sqi.csv
    │   └── experiment_2_summary.csv
    │
    ├── Experiment_3/
    │   ├── day4_experiment_3_window_sqi.csv
    │   └── experiment_3_summary.csv
    │
    ├── day4_all_window_sqi.csv
    └── day4_recording_summary.csv
```

Day 4 includes both experiment-specific files and combined Day 4 files so that we can either analyze each experiment separately or look at all Day 4 recordings together.

Older `*_beat_sqi.csv` files from a previous version of the pipeline have been moved to `FIU_Beat_Level_SQI/_archive_old_beat_sqi_aug31/`. They're no longer regenerated, since beat-level detail is now summarized into the window-level CSV instead.

The combined Day 4 files (`day4_all_window_sqi.csv`, `day4_recording_summary.csv`) list recordings shared by Experiments 2 and 3 only once, with Experiment = "Experiment_2 & Experiment_3".

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

They provide a starting point for separating more consistent beats from lower-quality beats, but they have not yet been treated as final validated thresholds. After Rutendo's review, λ now follows Zia et al. (25, with their D/L normalization), but `min_template_sqi`, `min_corr`, and `min_clipping_sqi` are still starting values rather than values tuned on FIU data.

Two specific open items from the current results:

- A small number of recordings still show extreme PI values, all on one channel -- Unpolarized_A_Green (hardware channel C5) -- where the raw signal itself has near-zero amplitude. This looks like a hardware/sensor issue on that channel rather than a problem with the PI calculation, but it hasn't been confirmed.
- Clipping SQI values across real beats cluster pretty tightly (roughly 0.55-0.94, median around 0.82) right around our 0.80 cutoff, for both good and bad beats. It may be worth revisiting whether that threshold, or the clipping formula itself, needs adjusting so it separates good and bad beats more clearly.

As more of the dataset is analyzed, the thresholds can be revisited again using the actual distributions of the SQI metrics.

---

# Future Work

The current beat-level pipeline gives us a foundation for more detailed signal-quality analysis.

Possible next steps include:

- analyze the distributions of each SQI metric across the full dataset
- refine the current SQI thresholds using the experimental data
- investigate the Unpolarized_A_Green (C5) PI outlier issue further -- confirm whether it's a hardware/sensor problem and decide whether those windows should be excluded or flagged
- revisit the clipping SQI threshold/formula given how tightly real clipping values cluster around the current 0.80 cutoff
- compare SQI distributions between skin tones
- compare signal quality across polarization conditions
- compare signal quality across wavelengths
- compare PI between good and bad windows
- determine whether PI should contribute to window-quality classification
- evaluate whether certain SQIs are more informative than others
- combine multiple SQIs using a machine-learning model
- investigate which polarization and wavelength combinations consistently provide the strongest PPG signals
- use the phantom results to guide future human testing

Ultimately, the goal is to build a more complete understanding of how the device performs across different optical and physical conditions and use that information to improve the reliability of future PPG measurements.