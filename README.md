# FIU Pulse Oximetry Analysis Pipeline

Signal-quality-index (SQI) pipeline for the FIU phantom PPG data -- filters each recording, splits it into 10-second windows, segments and classifies individual beats, and rolls that up into window-level and recording-level good/bad summaries.

For background on the study design, experiment variables, and a full walkthrough of how and why the pipeline works the way it does, see [PROJECT_INFO.md](PROJECT_INFO.md).

---

## Requirements

- Python 3
- `numpy`, `pandas`, `scipy`, `matplotlib`, `dtw-python`

```bash
pip install numpy pandas scipy matplotlib dtw-python
```

---

## Setup

1. Clone/pull this repo.
2. Make sure the FIU phantom data folders (`Experiment 1 Complete  copy`, `Experiment 2 Test (Day 3) copy`, `Experiment 2 & 3 (Day 4) copy`) are sitting alongside `beat_level_sqi.py`.
3. Only uncompressed `.json` recordings get processed -- `.json.gz` files are unzipped automatically if needed, and duplicates are skipped so nothing gets analyzed twice.

---

## Running the pipeline

```bash
python beat_level_sqi.py
```

This runs all three data sets -- Day 2 Experiment 1, Day 3 Experiment 2, and Day 4 Experiments 2 & 3 -- and saves results under `FIU_Beat_Level_SQI/`.

To run just one instead of all three, comment out the others in the `if __name__ == "__main__":` block at the bottom of `beat_level_sqi.py`, or call the function directly:

```python
from beat_level_sqi import process_experiment1_complete

process_experiment1_complete(
    experiment_root="Experiment 1 Complete  copy",
    output_root="FIU_Beat_Level_SQI/Day_2/Experiment_1"
)
```

---

## Output

Each data set saves two CSVs under `FIU_Beat_Level_SQI/`:

- a window-level CSV -- one row per 10-second window
- a recording-level summary CSV -- one row per recording, with %good/%bad windows

See [PROJECT_INFO.md](PROJECT_INFO.md) for what's actually in each file, and for the full pipeline write-up.
