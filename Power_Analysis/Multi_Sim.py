# -*- coding: utf-8 -*-
"""
Created on Tue May 20 23:16:26 2025

@author: kphud
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
from statsmodels.formula.api import ols
from statsmodels.stats.outliers_influence import summary_table
import os

np.random.seed(123)

# Setup
nsim = 1000
N = 72
temps = np.tile([14, 18, 22, 26, 30, 34], N // 6)
ids = [str(i+1) for i in range(N)]

# Constants
MSMRSD = 0.015 * np.sqrt(25)
CortSD = 2.85

# Storage
MSMR_models = []
Cort_models = []
Cort_MSMR_models = []

# Simulation loop
for _ in tqdm(range(nsim), desc="Simulations"):
    sex = np.random.choice(["M", "F"], size=N)
    mass = np.array([
        np.random.gamma(shape=(6.8**2)/(0.4**2), scale=(0.4**2)/6.8) if s == 'M' else
        np.random.gamma(shape=(20**2)/(2**2), scale=(2**2)/20) for s in sex
    ])
    
    MSMR_det = mass**-0.25 * np.exp(temps/10) * (0.2/6.04)
    Cort_det = mass**-0.25 * np.exp(temps/10) * (7.8/6.04)
    
    MSMR = np.random.gamma((MSMR_det**2)/(MSMRSD**2), (MSMRSD**2)/MSMR_det)
    Cort = np.random.gamma((Cort_det**2)/(CortSD**2), (CortSD**2)/Cort_det)

    df = pd.DataFrame({
        'ID': ids,
        'Temp': temps,
        'sex': sex,
        'Mass': mass,
        'MSMR': MSMR,
        'Cort': Cort
    })

    msmr_model = ols("np.log(MSMR) ~ np.log(Mass) + Temp", data=df).fit()
    cort_model = ols("np.log(Cort) ~ np.log(Mass) + Temp", data=df).fit()
    cort_msmr_model = ols("np.log(Cort) ~ np.log(MSMR)", data=df).fit()

    MSMR_models.append(msmr_model)
    Cort_models.append(cort_model)
    Cort_MSMR_models.append(cort_msmr_model)

# Save the last simulation
df.to_csv("SimulatedData.csv", index=False)

# Stats Summary
def average_coef_stats(models, predictors):
    means = {p: [] for p in predictors}
    rsq = []

    for model in models:
        rsq.append(model.rsquared)
        for p in predictors:
            means[p].append(model.params[p])
    
    means_avg = {p: np.mean(means[p]) for p in predictors}
    se = {p: np.std(means[p], ddof=1) for p in predictors}
    return means_avg, se, np.mean(rsq)

# Extract results
MSMR_avg, MSMR_se, MSMR_rsq = average_coef_stats(MSMR_models, ['Intercept', 'np.log(Mass)', 'Temp'])
Cort_avg, Cort_se, Cort_rsq = average_coef_stats(Cort_models, ['Intercept', 'np.log(Mass)', 'Temp'])
CortMSMR_avg, CortMSMR_se, CortMSMR_rsq = average_coef_stats(Cort_MSMR_models, ['Intercept', 'np.log(MSMR)'])

# Build stat table
stat_table = pd.DataFrame([
    [MSMR_avg['np.log(Mass)'], MSMR_se['np.log(Mass)'], '', '', MSMR_avg['Intercept'], MSMR_rsq],
    [MSMR_avg['Temp'], MSMR_se['Temp'], '', '', '', ''],
    [Cort_avg['np.log(Mass)'], Cort_se['np.log(Mass)'], '', '', Cort_avg['Intercept'], Cort_rsq],
    [Cort_avg['Temp'], Cort_se['Temp'], '', '', '', ''],
    [CortMSMR_avg['np.log(MSMR)'], CortMSMR_se['np.log(MSMR)'], '', '', CortMSMR_avg['Intercept'], CortMSMR_rsq],
], columns=["Estimate", "SE Est.", "T value", "p value", "Intercept", "R squared"], 
   index=["MSMR ~ mass", "MSMR ~ Temp", "Cort ~ mass", "Cort ~ Temp", "Cort ~ MSMR"])

# Round and format
stat_table = stat_table.round(3)
stat_table.loc[["MSMR ~ Temp", "Cort ~ Temp"], "Intercept"] = ""
stat_table.loc[["MSMR ~ Temp", "Cort ~ Temp"], "R squared"] = ""

# Save to CSV
os.makedirs("Power_Analysis", exist_ok=True)
stat_table.to_csv("StatsTab.csv")

# Plotting (Summary overlays)
def plot_model_lines(x, y, models, x_label, y_label, x_log=False, y_log=False, ylims=None):
    plt.figure(figsize=(5, 4))
    x_vals = np.linspace(min(x), max(x), 100)
    for model in models:
        intercept, slope = model.params[:2]
        y_vals = intercept + slope * np.log(x_vals) if x_log else intercept + slope * x_vals
        plt.plot(np.log(x_vals) if x_log else x_vals, y_vals, color='blue', alpha=0.02)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    if ylims:
        plt.ylim(ylims)
    plt.grid(False)
    plt.tight_layout()

plot_model_lines(df["Mass"], df["MSMR"], MSMR_models, "ln Body Mass (g)", "ln MSMR", x_log=True, y_log=True, ylims=(-7, 0))
plot_model_lines(df["Temp"], df["MSMR"], MSMR_models, "Temperature (C)", "ln MSMR", y_log=True)

plot_model_lines(df["Mass"], df["Cort"], Cort_models, "ln Body Mass (g)", "ln Cort", x_log=True, y_log=True, ylims=(-5, 2))
plot_model_lines(df["Temp"], df["Cort"], Cort_models, "Temperature (C)", "ln Cort", y_log=True)

plot_model_lines(df["MSMR"], df["Cort"], Cort_MSMR_models, "ln MSMR", "ln Cort", x_log=True, y_log=True)

plt.show()
