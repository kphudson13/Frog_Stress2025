# -*- coding: utf-8 -*-
"""
Created on Tue May 20 16:47:14 2025

@author: kphud
"""

# pip install pandas scipy statsmodels tqdm matplotlib seaborn


import numpy as np
import pandas as pd
from scipy.stats import gamma
import statsmodels.api as sm
from statsmodels.formula.api import ols
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns

np.random.seed(123)

nsim = 100
numvec = np.arange(6,251,3)

MSMRMassPowerVec = []
MSMRTempPowerVec = []
CortMassPowerVec = []
CortTempPowerVec = []
CortMSMRPowerVec = []

MSMRSD = 0.015 * np.sqrt(25)
CortSD = 2.85

# Simulation loop
for N in tqdm(numvec, desc="Running simulations"):
    MSMRMasspval = []
    MSMRTemppval = []
    CortMasspval = []
    CortTemppval = []
    CortMSMRpval = []
    
    for _ in range(nsim):
        temps = np.tile([14, 18, 22, 26, 30, 34], int(np.ceil(N/6)))[:N]
        sex = np.random.choice(['M', 'F'], size=N, replace=True)

        mass = np.array([
            gamma.rvs(a=(6.8**2)/(0.4**2), scale=(0.4**2)/6.8) if s == 'M' 
            else gamma.rvs(a=(20**2)/(2**2), scale=(2**2)/20)
            for s in sex
        ])

        MSMR_det = mass**-0.25 * np.exp(temps/10) * (0.2/6.04)
        Cort_det = mass**-0.25 * np.exp(temps/10) * (7.8/6.04)

        MSMR = gamma.rvs(a=(MSMR_det**2)/(MSMRSD**2), scale=(MSMRSD**2)/MSMR_det)
        Cort = gamma.rvs(a=(Cort_det**2)/(CortSD**2), scale=(CortSD**2)/Cort_det)

        df = pd.DataFrame({
            'Temp': temps,
            'sex': sex,
            'Mass': mass,
            'MSMR': MSMR,
            'Cort': Cort
        })

        df['logMass'] = np.log(df['Mass'])
        df['logMSMR'] = np.log(df['MSMR'])
        df['logCort'] = np.log(df['Cort'])

        MSMR_model = ols('logMSMR ~ logMass + Temp', data=df).fit()
        Cort_model = ols('logCort ~ logMass + Temp', data=df).fit()
        Cort_MSMR_model = ols('logCort ~ logMSMR', data=df).fit()

        MSMRMasspval.append(MSMR_model.pvalues['logMass'])
        MSMRTemppval.append(MSMR_model.pvalues['Temp'])
        CortMasspval.append(Cort_model.pvalues['logMass'])
        CortTemppval.append(Cort_model.pvalues['Temp'])
        CortMSMRpval.append(Cort_MSMR_model.pvalues['logMSMR'])

    # Calculate power
    MSMRMassPowerVec.append(np.mean(np.array(MSMRMasspval) < 0.05))
    MSMRTempPowerVec.append(np.mean(np.array(MSMRTemppval) < 0.05))
    CortMassPowerVec.append(np.mean(np.array(CortMasspval) < 0.05))
    CortTempPowerVec.append(np.mean(np.array(CortTemppval) < 0.05))
    CortMSMRPowerVec.append(np.mean(np.array(CortMSMRpval) < 0.05))

# Convert power lists to arrays
MSMRMassPowerVec = np.array(MSMRMassPowerVec)
MSMRTempPowerVec = np.array(MSMRTempPowerVec)
CortMassPowerVec = np.array(CortMassPowerVec)
CortTempPowerVec = np.array(CortTempPowerVec)
CortMSMRPowerVec = np.array(CortMSMRPowerVec)

# --- Plotting ---
plt.figure(figsize=(15, 10))
plot_labels = ['A', 'B', 'C', 'D', 'E']
titles = ['MSMR ~ Mass', 'MSMR ~ Temp', 'Cort ~ Mass', 'Cort ~ Temp', 'Cort ~ MSMR']
ydata = [MSMRMassPowerVec, MSMRTempPowerVec[:10], CortMassPowerVec,
         CortTempPowerVec[:10], CortMSMRPowerVec[:20]]
xdata = [numvec, numvec[:10], numvec, numvec[:10], numvec[:20]]

for i in range(5):
    plt.subplot(2, 3, i+1)
    sns.scatterplot(x=xdata[i], y=ydata[i], s=40)
    sns.regplot(x=xdata[i], y=ydata[i], lowess=True, scatter=False, color='red')
    plt.axhline(0.8, color='blue', linestyle='--')
    plt.xlabel("Sample size")
    plt.ylabel("Statistical power")
    plt.title(titles[i])
    plt.text(0.01, 0.95, plot_labels[i], transform=plt.gca().transAxes, fontsize=14, weight='bold')

plt.tight_layout()
plt.show()
