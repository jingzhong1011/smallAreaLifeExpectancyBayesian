import numpy as np
import pandas as pd
import pytensor.tensor as pt
from pandas import DataFrame
from scipy import stats
import pymc as pm
print(f'Running on PyMC v{pm.__version__}')
import arviz as az
import matplotlib.pyplot as plt
%config InlineBackend.figure_format = 'retina'
az.style.use("arviz-darkgrid")
# Boosting sampling using nutpie
import nutpie
# import jax

DeathCount = pd.read_csv('DeathCount_v5.csv', encoding = 'latin1')
# Only estimate TAIWAN REGION data
# That means Kingmen, Lianjiang would not be concluded
Data = DeathCount[(DeathCount['SEX'] == 1) & (DeathCount['Urbanicity_AST2023'] != 9)]
# lookup
lookup_urban = dict(zip(Data.Urbanicity_AST2023.unique(), range(len(Data.Urbanicity_AST2023.unique()))))
lookup_urban = Data.Urbanicity_AST2023.replace(lookup_urban).values

lookup_moi = dict(zip(Data.MOI.unique(), range(len(Data.MOI.unique()))))
lookup_moi = Data.MOI.replace(lookup_moi).values

lookup_agp = dict(zip(Data.AGP.unique(), range(len(Data.AGP.unique()))))
lookup_agp = Data.AGP.replace(lookup_agp).values

lookup_year = dict(zip(Data.YEAR.unique(), Data.YEAR.unique() - 2010))
lookup_year = Data.YEAR.replace(lookup_year).values
# To make sure the hierarchy working
lookup_moi_urban = Data[['MOI', 'Urbanicity_AST2023']].drop_duplicates()
lookup_moi_urban = lookup_moi_urban['Urbanicity_AST2023'] - 1
lookup_moi_urban = lookup_moi_urban.values

NB = pd.read_csv('NB_mat.csv', encoding = 'latin1')
NB_array = np.array(NB)

with pm.Model() as Model1:
    # Data
    YEAR = pm.Data('YEAR', lookup_year)
    AGP = pm.Data('AGP', lookup_agp)
    COUNT = pm.Data('COUNT', Data.COUNT)
    POP = pm.Data('POP', Data.POP)

    N_s = len(Data.MOI.unique())
    N_a = len(Data.AGP.unique())
    N_t = len(Data.YEAR.unique())
    
    # Priors
    alpha_0 = pm.Normal('alpha_0', mu = 0, sigma = 100000)
    beta_0 = pm.Normal('beta_0', mu = 0, sigma = 100000)

    sigma_alpha_s = pm.HalfNormal('sigma_alpha_s', sigma = 1)
    alpha_s = pm.ICAR("alpha_s", W=NB_array, sigma=sigma_alpha_s)
    
    sigma_beta_s = pm.HalfNormal('sigma_beta_s', sigma = 1)
    beta_s = pm.ICAR("beta_s", W=NB_array, sigma=sigma_beta_s)
    
    # Agp Terms
    sigma_alpha_a = pm.HalfNormal('sigma_alpha_a', sigma = 1)
    alpha_a_raw = pm.GaussianRandomWalk('alpha_a_raw', mu = 0, sigma = sigma_alpha_a, init_dist = pm.Normal.dist(0, 1), shape = N_a)
    alpha_a = pm.Deterministic('alpha_a', alpha_a_raw - alpha_a_raw[0])
    
    # Time Terms
    sigma_alpha_t = pm.HalfNormal('sigma_alpha_t', sigma = 1)
    alpha_t_raw = pm.GaussianRandomWalk('alpha_t_raw', mu = 0, sigma = sigma_alpha_t, init_dist = pm.Normal.dist(0, 1), shape = N_t)
    alpha_t = pm.Deterministic('alpha_t', alpha_t_raw - alpha_t_raw[0])

    # Interaction

    sigma_gamma = pm.HalfNormal('sigma_gamma', sigma = 1)
    gamma_raw = pm.GaussianRandomWalk('gamma_raw', mu = 0, sigma = sigma_gamma, 
                               init_dist = pm.Normal.dist(0, 1), shape = (N_a, N_t))
    gamma = pm.Deterministic('gamma', gamma_raw - gamma_raw[:, 0][:, None])

    sigma_xi = pm.HalfNormal('sigma_xi', sigma = 1)
    xi = pm.Normal('xi', mu = 0, sigma = sigma_xi, shape = (N_a, N_s))
    
    # Log Mean Model
    log_m = alpha_0 + beta_0 * YEAR
    log_m += alpha_s[lookup_moi] + beta_s[lookup_moi] * YEAR
    log_m += alpha_a[lookup_agp] + alpha_t[lookup_year]
    log_m += gamma[lookup_agp, lookup_year] + xi[lookup_agp, lookup_moi]

    # Likelihood
    mu = pm.math.exp(log_m) * POP
    d = pm.Poisson('d', mu=mu, observed=COUNT)
    
    #Compiled_model = nutpie.compile_pymc_model(Model1, backend='jax', gradient_backend='jax')
Compiled_model = nutpie.compile_pymc_model(Model1)
Trace_Model1 = nutpie.sample(
        Compiled_model,
        chains=4,
        draws=5000,
        tune=20000,
        target_accept=0.9,
        init_mean=0,
        seed=1011
    )

with Model1:
  PosteriorPredict_Model1 = pm.sample_posterior_predictive(Trace_Model1)
PP = PosteriorPredict_Model1['posterior_predictive']
PP.to_netcdf(path = "icar_poisson2/PPmale.nc")