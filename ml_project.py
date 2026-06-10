import os
import glob
import pandas as pd
import numpy as np

from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.model_selection import train_test_split


CLEANED_ROOT = "FIU_Cleaned_Data/Day_1"

STD_THRESHOLD = 0.02
SKEW_THRESHOLD = 0.0


# ============================================================
# LOAD ALL SQI FILES
# ============================================================

sqi_files = glob.glob(os.path.join(CLEANED_ROOT, "**", "sqi.csv"), recursive=True)

if len(sqi_files) == 0:
    raise FileNotFoundError(f"No sqi.csv files found under {CLEANED_ROOT}")

all_dfs = []

for file in sqi_files:
    df = pd.read_csv(file)
    df["source_file"] = file
    all_dfs.append(df)

data = pd.concat(all_dfs, ignore_index=True)

print("\nLoaded SQI files:", len(sqi_files))
print("Combined dataset shape:", data.shape)
print("\nColumns:")
print(data.columns.tolist())


# ============================================================
# CREATE WEAK GOOD/BAD LABELS
# ============================================================

def make_label(row):
    if row["Std"] >= STD_THRESHOLD and row["Skewness"] > SKEW_THRESHOLD:
        return 1   # good
    else:
        return 0   # bad

data["label"] = data.apply(make_label, axis=1)

print("\nLabel counts:")
print(data["label"].value_counts())
print("\nLabel meaning: 1 = good, 0 = bad")


# ============================================================
# SELECT FEATURES
# ============================================================

feature_cols = [
    "Mean",
    "Std",
    "Skewness",
    "Channel",
    "Hardware Channel",
    "SkinTone",
    "Speed",
    "Depth",
    "Experiment",
]

# Keep only columns that actually exist
feature_cols = [c for c in feature_cols if c in data.columns]

model_df = data[feature_cols + ["label", "source_file"]].copy()

# Drop rows with missing feature values
model_df = model_df.dropna()

print("\nUsing features:")
print(feature_cols)


# ============================================================
# ONE-HOT ENCODE CATEGORICAL FEATURES
# ============================================================

categorical_cols = model_df.select_dtypes(include=["object"]).columns.tolist()

# Do not encode source_file as a model feature
if "source_file" in categorical_cols:
    categorical_cols.remove("source_file")

model_encoded = pd.get_dummies(
    model_df,
    columns=categorical_cols,
    drop_first=False
)

X = model_encoded.drop(columns=["label", "source_file"])
y = model_encoded["label"]

# ============================================================
# TRAIN/TEST SPLIT
# ============================================================

# Basic split for first test
X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.25,
    random_state=42,
    stratify=y
)

# ============================================================
# TRAIN RANDOM FOREST
# ============================================================

model = RandomForestClassifier(
    n_estimators=200,
    random_state=42,
    class_weight="balanced"
)

model.fit(X_train, y_train)

# ============================================================
# EVALUATE MODEL
# ============================================================

preds = model.predict(X_test)

print("\nConfusion Matrix:")
print(confusion_matrix(y_test, preds))

print("\nClassification Report:")
print(classification_report(y_test, preds, target_names=["bad", "good"]))

# ============================================================
#  FEATURE IMPORTANCE
# ============================================================

importances = pd.DataFrame({
    "feature": X.columns,
    "importance": model.feature_importances_
}).sort_values("importance", ascending=False)

print("\nTop 20 Most Important Features:")
print(importances.head(20))


# ============================================================
# SAVE COMBINED DATASET
# ============================================================

data.to_csv("day1_combined_sqi_with_labels.csv", index=False)
importances.to_csv("day1_feature_importances.csv", index=False)

print("\nSaved:")
print("day1_combined_sqi_with_labels.csv")
print("day1_feature_importances.csv")