# 🏥 Fetal Healthcare Data — Project README

![Dashboard](<img width="1919" height="720" alt="Screenshot 2026-01-22 185013" src="https://github.com/user-attachments/assets/35b65b15-9cc4-4c1f-be65-1d7ea1b4d7d0" />
) 
---

## 📌 Project Overview

This project analyses **Cardiotocography (CTG)** examination data to classify fetal health status into three categories: **Normal**, **Suspect**, and **Pathological**. The Excel workbook contains raw clinical data, derived calculated fields, a risk scoring engine, and an interactive Excel dashboard — all built to support prenatal monitoring and early detection of fetal distress.

The dashboard (see screenshot) provides a visual summary of key health indicators across the three case types, enabling clinicians and analysts to quickly identify patterns and anomalies.

---

## 📂 Workbook Structure

The Excel file contains **4 sheets**:

| Sheet | Purpose |
|---|---|
| `fetal_health` | Core dataset — 2,126 records × 22 clinical features |
| `data set` | Extended version of the dataset with added derived columns (risk scoring, composite index, label mapping) |
| `CAL` | Pivot/calculation layer powering the dashboard charts |
| `dash` | Excel dashboard with interactive visuals |

---

## 📊 Dataset Description

### Source Sheet: `fetal_health`
- **Total Records:** 2,126
- **Features:** 22 columns

### Class Distribution (Target Variable: `fetal_health`)

| Class Code | Label | Count | % of Total |
|---|---|---|---|
| 1 | Normal Case | 1,655 | 77.8% |
| 2 | Suspect Case | 295 | 13.9% |
| 3 | Pathological Case | 176 | 8.3% |

> ⚠️ The dataset is **imbalanced** — Normal cases dominate. Any machine learning model trained on this data should account for class imbalance (e.g., using SMOTE, class weights, or stratified sampling).

---

## 🔬 Feature Descriptions

### Cardiotocography (CTG) Signals

| Column Name | Description |
|---|---|
| `baseline value` | Fetal Heart Rate (FHR) baseline in beats per minute (bpm) |
| `accelerations` | Number of accelerations per second |
| `fetal_movement` | Number of fetal movements per second |
| `uterine_contractions` | Number of uterine contractions per second |
| `light_decelerations` | Number of light decelerations (LDs) per second |
| `severe_decelerations` | Number of severe decelerations (SDs) per second |
| `prolongued_decelerations` | Number of prolonged decelerations (PDs) per second |

### Short-Term Variability

| Column Name | Description |
|---|---|
| `abnormal_short_term_variability` | Percentage of time with abnormal short-term variability |
| `mean_value_of_short_term_variability` | Mean value of short-term variability |

### Long-Term Variability

| Column Name | Description |
|---|---|
| `percentage_of_time_with_abnormal_long_term_variability` | % of time with abnormal long-term variability |
| `mean_value_of_long_term_variability` | Mean value of long-term variability |

### Histogram Features (FHR Distribution)

| Column Name | Description |
|---|---|
| `histogram_width` | Width of the FHR histogram |
| `histogram_min` | Minimum value of FHR histogram |
| `histogram_max` | Maximum value of FHR histogram |
| `histogram_number_of_peaks` | Number of histogram peaks |
| `histogram_number_of_zeroes` | Number of histogram zeroes |
| `histogram_mode` | Histogram mode |
| `histogram_mean` | Histogram mean |
| `histogram_median` | Histogram median |
| `histogram_variance` | Histogram variance |
| `histogram_tendency` | Histogram tendency (left=-1, symmetric=0, right=1) |

### Target Variable

| Column Name | Description |
|---|---|
| `fetal_health` | **Target class** — 1 = Normal, 2 = Suspect, 3 = Pathological |

---

## ⚙️ Derived & Engineered Columns (`data set` sheet)

The `data set` sheet extends the raw data with three engineered columns:

### 1. `FETAL_HEALTH2` — Human-Readable Label
Maps the numeric `fetal_health` code to a string label:
```
=IF(fetal_health=1, "Normal_Case",
  IF(fetal_health=2, "Suspect_Case",
    IF(fetal_health=3, "Pathological_Case")))
```

### 2. `RISK SCORING` — Rule-Based Risk Flag
A conditional risk classification based on key clinical thresholds:
```
=IF(severe_decelerations > 0.5,       "High Risk",
  IF(prolongued_decelerations > 0.3,  "Medium Risk",
    IF(accelerations < 0.1,           "Low variability",
      IF(abnormal_short_term_variability > 20, "High Risk", ""))))
```

### 3. `composite index` — Weighted Risk Score
A numeric composite score aggregating the most clinically significant features:
```
= (0.30 × severe_decelerations)
+ (0.20 × prolongued_decelerations)
+ (0.20 × abnormal_short_term_variability)
+ (0.15 × percentage_of_time_with_abnormal_long_term_variability)
+ (0.15 × (1/accelerations + 1/fetal_movement))
```

**Weight rationale:**
- Severe & prolonged decelerations carry the highest weight (clinical priority indicators)
- Short-term variability accounts for 20% (key distress marker)
- Long-term variability and movement/acceleration each contribute 15%

---

## 📈 Dashboard Overview (`dash` sheet)

The dashboard visualises aggregated statistics from the `CAL` pivot sheet and presents the following panels:

### KPI Cards (Top Row)
| Metric | Value |
|---|---|
| Severe Deceleration (avg) | 3.3E-06 |
| Mean Acceleration (avg) | 0.0032 |
| Histogram Variance (avg) | 18.91 |
| Mean Short-Term Variability (avg) | 3.3E-06 |

### Charts & Visuals

| Panel | Chart Type | Description |
|---|---|---|
| Fetal Health Demography | Donut Chart | Case distribution: 1,646 Normal / 292 Suspect / 175 Pathological |
| Abnormal Short-Term Variability | Bar Chart | Comparison across Normal, Pathological, Suspect cases |
| Mean Long-Term Variability | Line Chart | Trend of long-term variability by health class |
| % Time with Abnormal Short-Term Variability | Donut Chart | Proportion breakdown by health class |
| Histogram Variability | Scatter Plot | Variance spread across case types |
| Prolonged Deceleration | Line Chart | Average prolonged deceleration by class |
| Severe Deceleration | Line Chart | Average severe deceleration by class |
| Light Deceleration | Line Chart | Average light deceleration by class |
| Average Fetal Movement | Line Chart | Fetal movement trends by class |
| Average Acceleration | Line Chart | Acceleration averages by class |

### Slicers / Filters
The dashboard includes interactive slicers for:
- `Normal_Case`
- `Pathological_Case`
- `Suspect_Case`

---

## 🔄 Data Processing Workflow

```
Raw CTG Data (fetal_health sheet)
         │
         ▼
   data set sheet
   ┌─────────────────────────────────┐
   │ + FETAL_HEALTH2 (label mapping) │
   │ + RISK SCORING (rule-based)     │
   │ + composite index (weighted)    │
   └─────────────────────────────────┘
         │
         ▼
   CAL sheet (Pivot Tables & Aggregations)
         │
         ▼
   dash sheet (Interactive Dashboard)
```

---

## 🚀 How to Use

1. **Open the workbook** in Microsoft Excel (2016 or later recommended for full slicer/chart support).
2. **Navigate to the `dash` sheet** to explore the interactive dashboard.
3. **Use the slicers** on the left panel (Normal / Pathological / Suspect) to filter all charts simultaneously.
4. **Explore raw data** in the `fetal_health` or `data set` sheets for column-level analysis.
5. **Check the `CAL` sheet** to review underlying pivot table calculations feeding the dashboard.

---

## 📋 Key Clinical Insights

- **77.8% of cases are Normal** — the dataset reflects real-world prevalence but introduces class imbalance challenges for modelling.
- **Severe decelerations** are the strongest single indicator of fetal distress and carry the highest weight in the composite risk index.
- **Pathological cases** show markedly higher `abnormal_short_term_variability` and lower `mean_value_of_long_term_variability` compared to Normal cases.
- **Prolonged decelerations** are rare but clinically critical; they receive a dedicated "Medium Risk" flag in the risk scoring logic.

---

## 📁 File Information

| Property | Value |
|---|---|
| File Name | `fetal_health_-_Copy_AutoRecovered__L.xlsx` |
| Format | Microsoft Excel Workbook (.xlsx) |
| Total Sheets | 4 |
| Primary Dataset Rows | 2,126 |
| Primary Dataset Columns | 22 |
| Extended Dataset Columns | 26 (inc. 4 derived columns) |

---

## 👤 Author
Sunday Abdulsalam
Linkedin: https://www.linkedin.com/in/sundayabdulsalam/


---

