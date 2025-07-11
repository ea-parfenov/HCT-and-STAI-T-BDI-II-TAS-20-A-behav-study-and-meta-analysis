"""
Calculate effect sizes of the folowwing studies
1. Rae et al. 2019
2. Critchley et al. 2023
3. Arslanova et al., 2023
4. Quadt et al., 2021
5. Leganes-Fonteneau et al., 2019
6. Bokk and Foster, 2022
7. Botan et al., 2021
"""

import pandas as pd
import numpy as np
import pyreadstat
from numpy.random import randint
from scipy.stats import norm, spearmanr
from joblib import Parallel, delayed

# Open Critchley et al., 2023 data
critchley = pd.read_csv('Critchley_2023.csv')
# Drop participants without no HCT and STAI-T scores
critchley = critchley.dropna(subset=['HBTacc', 'STAI_Y2'])
critchley_healthy = critchley[critchley['Primary_psychiatric_diagnosis'] == 'Comparison']
critchley_patients = critchley[critchley['Primary_psychiatric_diagnosis'] != 'Comparison']
# Open Rae et al., 2019 data
rae = pd.read_csv('Rae_2019.csv')
rae_healthy = rae[rae['Group'] == 2]
rae_patients = rae[rae['Group'] == 1]
# Open Arslanova et al., 2023 data
arslanova = pd.read_csv('Arslanova_2023.csv').dropna()
# Open Quadt et al., 2021 data
quadt = pd.read_csv('Quadt_2021.csv').dropna()
# Open Leganes-Fonteneau et al., 2019 data
leganes = pd.read_csv('Leganes_2019.csv')
# Open Bokk and Foster, 2022
bokk = pd.read_csv('Bokk_2022.csv')
# open Botan et al., 2021
botan, meta = pyreadstat.read_sav("Botan_2021.sav")
botan = botan[['Gender', 'Tracking', 'Anxiety_Trait']].dropna()


# Spearman correlation with CI
def get_boot_r(x, y):
    idx = randint(0, y.shape[0], y.shape[0])
    return spearmanr(y[idx], x[idx])[0]


def get_jack_r(x, y, i):
    y = np.delete(y, i, axis=0)
    x = np.delete(x, i, axis=0)
    return spearmanr(x, y)[0]


def spearman_ci(x, y, n_samples=10000, alpha=0.05):
    """
    Calculate confidence interval at specified alpha.

    Parameters
    ----------
    x : Array of floats
        First variable.
    y : Array of floats
        Second variable.
    n_samples : int, optional
        Number of bootstrap permutations to run. The default is 10000.
    alpha : float, optional
        Alpha level for confidence interval. The default is 0.05.

    Returns
    -------
    None.

    """
    alphas = np.array([alpha/2, 1-alpha/2])
    x = np.array(x)
    y = np.array(y)
    orig = spearmanr(x, y)[0]
    boot = sorted(Parallel(n_jobs=8, verbose=0)(delayed(get_boot_r)(x, y) for i in range(n_samples)))
    z0 = norm.ppf((1.0*np.sum(boot < orig, axis=0)) / n_samples)
    jstat = Parallel(n_jobs=8, verbose=0)(delayed(get_jack_r)(x, y, i) for i in range(len(y)))
    jmean = np.mean(jstat, axis=0)
    a = np.sum((jmean - jstat)**3, axis=0) / (6.0 * np.sum((jmean - jstat)**2, axis=0)**1.5)
    zs = z0 + norm.ppf(alphas).reshape(alphas.shape+(1,)*z0.ndim)
    avals = norm.cdf(z0 + zs/(1-a*zs))
    nvals = np.around((n_samples-1)*avals).astype('int')
    return boot[nvals[0]], boot[nvals[1]]


def spearman(x, y, x_lab, y_lab):
    """
    Calculate spearman correlation between two variables with 95% confidence
    interval.
    """
    # All values
    r, p = spearmanr(x, y)
    ci = spearman_ci(x, y)
    n = len(x)

    print(f'{x_lab} vs {y_lab} correlation')
    print(f'rho = {r:.2f}, p = {p:.3f}, ci = {ci[0]:.2f} {ci[1]:.2f}, n = {n}', '\n')


# Calculate correlation for Critchley et al., 2023
for df in [critchley_healthy, critchley_patients, critchley]:
    spearman(df['HBTacc'], df['STAI_Y2'], 'IAcc', 'STAI-T')
    print("% Female :", len(df.loc[df['sex_m1f2o3'] == 2]) / len(df['sex_m1f2o3'])*100)

# Calculate correlation for Rae et al., 2019
for df in [rae_healthy, rae_patients, rae]:
    spearman(df['Interoceptive_accuracy_tracking'], df['STAI_trait'], 'IAcc', 'STAI-T')

# Calculate correlation for Arslanova et al., 2023
spearman(arslanova['HBacc'], arslanova['TAS2'], 'IAcc', 'TAS-20')
# Calculate correlation for Quadt et al., 2023
spearman(quadt['T0_Track_meanacc'], quadt['T0_TAS_total'], 'IAcc', 'TAS-20')
spearman(quadt['T0_Track_meanacc'], quadt['T0_STAI_T'], 'IAcc', 'STAI-T')
# Calculate correlation for Leganes-Fonteneau et al., 2019
spearman(leganes['HT_SCORE'], leganes['TAS'], 'IAcc', 'TAS-20')
# Calculate correlation for Bokk and Forest, 2022
spearman(bokk['Hearbeat Score'], bokk['STAI'], 'IAcc', 'STAI-T')
# Calculate correlation for Boran et al., 2021
spearman(botan['Tracking'], botan['Anxiety_Trait'], 'IAcc', 'STAI-T')
print("% Female", len(botan.loc[botan['Gender'] == 1]) / len(botan['Gender'])*100)
