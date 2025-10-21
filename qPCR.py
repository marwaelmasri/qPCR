# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 11:17:28 2025

@author: marwa

This script analyses and plots qPCR data
"""

import pandas as pd
import numpy as np
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt


# Open excel.
data = pd.read_excel(r"C:\Users\marwa\OneDrive\Documents\2025-10-15_ME_IPSC_MSN.xlsx",
                     sheet_name="Results")

# Define Column names.
header = data.iloc[43,:]

# Subset the data to useful rows.
qpcr=data.iloc[44:160,:]

# Set header as column names.
qpcr.columns=header

# Also filter out wells you want to exclude if you find them to be outliers.
# Alternatively filter out i excel before import.

# Reset index
qpcr = qpcr.reset_index(drop=True)


# Create new column based on substring match.
qpcr["celltype"] = np.where(qpcr["Sample Name"].str.contains("IPSC"), "iPSC",
                            np.where(qpcr["Sample Name"].str.contains("MSN"), "MSN",
                            qpcr["Sample Name"]))

# Check qpcr.dtypes.
# CT is object, change to float.
# First remove undetermined.
qpcr_values = qpcr[qpcr["CT"] != "Undetermined"].copy()
qpcr_values["CT"] = qpcr_values["CT"].astype(float)
# check qpcr_values.dtypes to make sure.

# Removing Outliers.
# Get Q1 and Q3.
# Filter rows where Target Name = target, celltype = cell.
# Can add condition if needed.

def get_quartiles(row, df):
    # Filter rows matching current row's Target Name and celltype
    mask = (df["Target Name"] == row["Target Name"]) & (df["celltype"] == row["celltype"]) & (df["Sample Name"] == row["Sample Name"])
    subset = df.loc[mask, "CT"]
    
    if subset.empty:
        return pd.Series({"Q1": None, "Q3": None})
    
    return pd.Series({
        "Q1": subset.quantile(0.25),
        "Q3": subset.quantile(0.75)
    })

# Apply Function across all rows.
qpcr_values[["Q1", "Q3"]] = qpcr_values.apply(get_quartiles, df=qpcr_values, axis=1)


# Get IQR.
qpcr_values["IQR"] = qpcr_values["Q3"] - qpcr_values["Q1"]

# Multiple IQR by 1.5.
qpcr_values["OutlierIf"] = qpcr_values["IQR"]*1.5

# Filter out outliers, assign to new variable.
qpcr_OL_removed = qpcr_values[(qpcr_values["CT"] > (qpcr_values["Q1"] - qpcr_values["OutlierIf"]))  &
 (qpcr_values["CT"] < (qpcr_values["Q3"] + qpcr_values["OutlierIf"]))].copy()


# HK CT averages.
def get_avg(row, df, HK):
    # Filter where target == HK (fixed), and celltype matches the current row
    mask = (df["Target Name"] == HK) & (df["celltype"] == row["celltype"]) & (df["Sample Name"] == row["Sample Name"])
    subset = df.loc[mask, "CT"]

    # Return mean or blank if no matches
    return subset.mean() if not subset.empty else ""

# Apply for EIF4A2.
qpcr_OL_removed["EIF4A2"] = qpcr_OL_removed.apply(lambda r: get_avg(r, qpcr_OL_removed, "EIF4A2"), axis=1)


# Apply for UBC.
qpcr_OL_removed["UBC"] = qpcr_OL_removed.apply(lambda r: get_avg(r, qpcr_OL_removed, "UBC"), axis=1)


# Geomean of HK.
qpcr_OL_removed["geomean"] = stats.gmean(qpcr_OL_removed[["EIF4A2", "UBC"]], axis=1)

# dct
qpcr_OL_removed["dct"] = qpcr_OL_removed["CT"] - qpcr_OL_removed["geomean"]

# dct.con with a function, con is the control you want to nomralise to.
def get_dctcon(row, df, con):
    # Filter where target == HK (fixed), and celltype matches the current row
    mask = (df["Target Name"] == row["Target Name"]) & (df["celltype"] == con)
    subset = df.loc[mask, "dct"]

    # Return mean or blank if no matches
    return subset.mean() if not subset.empty else ""

# Apply normalisation to control for all rows.
qpcr_OL_removed["dct.con"] = qpcr_OL_removed.apply(lambda r: get_dctcon(r, qpcr_OL_removed, "iPSC"), axis=1)


# ddct
qpcr_OL_removed["ddct"] = qpcr_OL_removed["dct"] - qpcr_OL_removed["dct.con"]


# 2^-ddct
qpcr_OL_removed["2^-ddct"] = 2**(-qpcr_OL_removed["ddct"])



# export csv
qpcr_OL_removed.to_csv(r"C:\Users\marwa\OneDrive\Documents\qPCR-preFANS.csv", index=False)

# To plot:
# remove HK
final = qpcr_OL_removed[(qpcr_OL_removed["Target Name"] != "EIF") & (qpcr_OL_removed["Target Name"] != "UBC")]

# plot
sns.barplot(
    data=final,
    x="Target Name",
    y="2^-ddct",
    hue="celltype",
    palette=["#E8E6E6", "#9D67E6"],
    dodge=True,
    edgecolor="black",
    linewidth=0.5, 
    errorbar="sd",
    order=["S100B", "BCL", "PPP", "RBFOX3"]
)

# Overlay data points
sns.stripplot(
    data=final,
    x="Target Name",
    y="2^-ddct",
    hue="celltype",
    dodge=True,
    palette=["#000000", "#000000"],  # black
    size=4,
    jitter=True,                     
    alpha=0.6                        # transparency
)

plt.xlabel("Marker")

# Adjust legend to avoid duplication
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles[:2], labels[:2], title="Cell Type",
           bbox_to_anchor=(1.02, 1),
           loc="upper left",
           borderaxespad=0)

plt.tight_layout()
plt.show()

