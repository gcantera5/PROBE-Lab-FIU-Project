# Lab Update – SQI Pipeline Changes

Quick recap of what's changed on the signal-quality-index side since the last update, based on the back-and-forth with Mark and Rutendo.

## Why we changed things

The original pipeline was doing beat-by-beat comparison against a template, but after talking it through with Mark and Rutendo we tried moving to a simpler window-level approach instead — filter the signal, calculate PI, normalize the whole 10-second window, run a couple of SQI checks on that window, and classify it as good or bad directly, no beat segmentation at all.

When I ran that version on real data, basically everything came back as "good." Dug into it and found two problems: the window-level classification was only checking clipping, which is a pretty lenient check on its own, so almost nothing failed it. And separately, PI was blowing up into the hundreds/thousands of percent on some windows because I was calculating it off the *filtered* signal — the bandpass filter strips out the DC baseline, so once that's gone the PI math (AC/DC) doesn't mean anything anymore.

## Getting realigned with Mark and Rutendo

Flagged both issues to them, and that's when Mark clarified two things that changed the plan: we still need beat-to-beat comparison happening (not just one check on the whole window), and normalization should happen per beat, not per window — he'd been picturing per-beat normalization the whole time, which is actually closer to the original pipeline than the window-only version we'd built.

So the pipeline landed here (Mark signed off on this version):

- Filter the window (Chebyshev II bandpass) — this filtered version is only used for beat segmentation and beat comparison
- Calculate PI from the **raw, unfiltered** window signal, not the filtered one
- Segment the filtered window into individual beats (valley-to-valley)
- Normalize each beat individually (0–1 scale)
- Compare each normalized beat to a template beat (DTW distance, correlation, MAD, clipping)
- Classify each beat as good/bad
- Roll that up into a percent-bad-beats number for the window, and classify the whole window as good/bad off of that (currently: bad if more than 30% of its beats are bad)
- Export everything to CSV per window, then a recording-level summary with %good/%bad windows

Mark did note he's still a little unsure about calculating PI on an unfiltered window vs. a "cleaner" one, but didn't have a clean way to do that without detrending messing with the filter, so we agreed to move forward with this and revisit later if it becomes an issue.

## Implementing it + what we're seeing on real data

Got this version implemented and running across Day 2/3/4 recordings. Current aggregate result is **33.12% good windows / 66.88% bad windows**, and it's fully reproducible — reran it and got the exact same split, so the pipeline itself is deterministic, nothing random going on there.

A couple things came out of digging into the real results:

**PI outliers:** There's a small set of windows (about 116) with extreme PI values. Traced almost all of them (105/116) to one specific channel — Unpolarized_A_Green (hardware channel C5) — where the raw signal itself has near-zero amplitude in those windows. This looks less like a pipeline bug and more like a hardware/sensor issue on that specific channel, but wanted to flag it before assuming that.

**Threshold retuning:** Since beat comparison switched from comparing raw-amplitude beats to comparing normalized (0–1) beats, some of our classification thresholds (sqi_lambda, min_template_sqi, max_mad) were tuned for the old raw-amplitude scale and were basically dead weight on the new scale — they weren't triggering across the whole dataset. Pulled real percentile distributions from ~43,000 actual beats and retuned those three thresholds to match. Reran with the new numbers and the overall good/bad split barely moved (33.03%/66.97% vs. 33.12%/66.88%), so given how small the difference was, decided it wasn't worth carrying non-standard threshold values and reverted back to the original numbers. Documented the whole retune-then-revert reasoning in the code in case we want to revisit it.

## Open questions for Mark and Rutendo

- Is the Unpolarized_A_Green (C5) PI issue something on the hardware/sensor side that we should just exclude or flag, or is there a fix we should be trying on the processing side?
- Our clipping_sqi values are clustering pretty tightly (roughly 0.55–0.94, median around 0.82) right around our 0.80 cutoff, for both good and bad beats. Worth rethinking whether that threshold or the underlying clipping formula needs adjusting so it's actually separating good from bad rather than sitting right in the middle of the distribution.
