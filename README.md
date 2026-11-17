# MCMC from First Principles: Theory, Proofs and a MATLAB Implementation

An applied and statistical mathematics project written in **MATLAB** (core
functions only, no toolboxes). It implements random-walk Metropolis, a Gibbs
sampler, parallel tempering and convergence diagnostics from scratch, proves the
theory behind them, and checks every result numerically against an exact answer.
A short R script cross-checks the diagnostics with the `coda` package.

The aim is not just to run MCMC but to understand why it works, how accurate it
is, how to tune it, and how to tell when it has failed.

## Questions

1. Why does a Metropolis-Hastings or Gibbs sampler converge to the right distribution?
2. How accurate is an MCMC average? How many independent draws is a correlated chain worth?
3. How should the random-walk step size be chosen, and where does the "0.234 acceptance rate" come from?
4. How do you detect a sampler that is stuck, and what fixes it?
5. Does the sampler reproduce a posterior that can be computed exactly?

## Key findings

*(Fill in after running: how closely the estimated autocorrelation times match
the exact formula; the optimal step size and acceptance rate found in 50
dimensions compared with 2.38 and 0.234; what R-hat and effective sample size
showed on the bimodal target with and without parallel tempering; how closely the
MCMC regression posterior matched the exact one, and the ESS gain from a
covariance-tuned proposal.)*

## Mathematics covered

Everything marked "proved" is proved in `reports/report.md`; items marked
"cited" are standard theorems stated with a reference.

| Topic | Result | Status |
|---|---|---|
| Validity | Detailed balance implies stationarity; the Metropolis-Hastings kernel satisfies detailed balance; Gibbs is Metropolis-Hastings with acceptance 1 | proved |
| Convergence and CLT | Convergence of the chain and the Markov-chain central limit theorem | cited |
| Accuracy | N Var(mean) tends to sigma^2 tau, where tau = 1 + 2 sum of autocorrelations | proved |
| AR(1) chain | tau = (1 + phi) / (1 - phi) | proved |
| Gibbs, bivariate normal | The x-coordinates form an AR(1) chain with phi = rho^2, so tau = (1 + rho^2) / (1 - rho^2) | proved |
| Step size | E[min(1, e^W)] = 2 Phi(-sigma/2) for W ~ N(-sigma^2/2, sigma^2); high-dimensional acceptance rate 2 Phi(-l/2) | proved |
| Optimal scaling | Efficiency proportional to l^2 * 2 Phi(-l/2) peaks at l = 2.38 with acceptance 0.234 | optimum computed numerically; diffusion-limit argument cited |
| Diagnostics | Split R-hat, effective sample size (Geyer initial positive sequence) | implemented; properties cited |
| Parallel tempering | Swap acceptance exp((b_k - b_k+1)(log pi(x_k+1) - log pi(x_k))) preserves the product target | proved |
| Bayesian regression | Exact normal-inverse-gamma posterior, including marginal means and variances | proved |

## Experiments

| Script | What it does |
|---|---|
| `run_01_verify_theory.m` | (a) builds the exact Metropolis-Hastings matrix on a 5-state space and checks stationarity, detailed balance and the decay rate; (b) compares estimated autocorrelation times of the Gibbs sampler with (1 + rho^2) / (1 - rho^2) for rho = 0.5, 0.9, 0.99, and the variance of chain averages with sigma^2 tau / N over 1,000 independent chains; (c) checks the acceptance lemma by Monte Carlo |
| `run_02_scaling_acceptance.m` | Random-walk Metropolis on a 50-dimensional standard normal over a grid of step sizes: acceptance rate versus 2 Phi(-l/2), mean squared jump distance and ESS, and the empirical versus theoretical optimum |
| `run_03_multimodal.m` | A 10-dimensional two-mode target with the modes about 16 standard deviations apart. Four chains started in different modes: random walks of several step sizes versus parallel tempering with 10 temperatures. Reports R-hat, ESS, per-chain mode occupancy and mode switches |
| `run_04_bayes_regression.m` | Bayesian linear regression with unknown noise variance and two strongly correlated predictors. Isotropic versus covariance-tuned proposals, four dispersed chains, every posterior mean and standard deviation compared with the exact values using Monte Carlo standard errors |

One lesson the third experiment is designed to show: **per-chain effective
sample size can look excellent while the chains are stuck in different modes**.
Only the between-chain comparison (R-hat) exposes the problem.

## Tools

- MATLAB, core functions only (`randn`, `fft`, `chol`, `eig`, `erfc`, `fminbnd`)
- R with the `coda` package, for the cross-check of ESS and R-hat

## Reproduce

1. In MATLAB, set the current folder to the project root and run `scripts/run_all.m`.
   It runs `tests/run_tests.m` first and then the four experiments; allow a few minutes.
2. Results appear in `results/` (CSV files) and `results/figures/` (PNG files).
3. Optional, in R from the project root:
   `install.packages("coda")`, then `source("R/cross_check_coda.R")`.

To run a single experiment, for example `run('scripts/run_02_scaling_acceptance.m')`.
Every script sets a random seed.

## Project structure

```
src/
  mh_sampler.m            random-walk Metropolis (scalar or matrix proposal factor)
  gibbs_bivariate.m       Gibbs sampler for a bivariate normal
  parallel_tempering.m    replica-exchange sampler
  autocorr_fft.m          sample autocorrelation via the FFT
  iact.m                  integrated autocorrelation time (Geyer initial positive sequence)
  ess.m                   effective sample size
  rhat.m                  split R-hat
  nig_posterior.m         exact conjugate regression posterior
  log_post_regression.m   log posterior of (beta, log sigma^2) with the Jacobian term
  log_mix_d.m             log density of the bimodal target
  summarise_modes.m       mode-occupancy and switching diagnostics
  normal_cdf.m, write_csv.m   small utilities
scripts/                  run_01 to run_04 and run_all.m
tests/run_tests.m         unit tests with known answers
R/cross_check_coda.R      cross-check of ESS and R-hat with coda
reports/report.md         theorems, proofs, experiments, limitations, references
results/                  generated tables and figures
```

## Tests

`tests/run_tests.m` checks, against known answers: the autocorrelation time of
iid noise (1) and of an AR(1) chain (9 for phi = 0.8); R-hat for identical and for
separated chains; Metropolis moments and acceptance rate for a standard normal;
Gibbs correlation and variance; the acceptance lemma; the closed-form regression
posterior; and agreement between Metropolis output and the exact posterior mean.

## Limitations

- The central limit theorem for Markov chains, the diffusion limit behind the
  0.234 rule and the properties of R-hat are cited, not proved
- The 0.234 rule is asymptotic and can be a poor guide for skewed, heavy-tailed or
  strongly dependent targets
- ESS estimation is noisy for short chains or very slow mixing, and is computed
  per coordinate; the ESS of the whole distribution can be smaller
- The proposal covariance is estimated from a pilot run and then frozen, which keeps
  the chain Markov but is not fully adaptive
- The tempering ladder was chosen by hand for one target; real problems need tuning
  or adaptive ladders
- All targets are low-dimensional with known answers so the code can be graded;
  real posteriors are rarely so well behaved

## References

- Tierney (1994), Markov chains for exploring posterior distributions, *Annals of Statistics*
- Roberts, Gelman and Gilks (1997), Weak convergence and optimal scaling of random walk Metropolis algorithms, *Annals of Applied Probability*
- Roberts and Rosenthal (2004), General state space Markov chains and MCMC algorithms, *Probability Surveys*
- Gelman and Rubin (1992), Inference from iterative simulation using multiple sequences, *Statistical Science*
- Geyer (1992), Practical Markov chain Monte Carlo, *Statistical Science*
- Vehtari et al. (2021), Rank-normalization, folding, and localization: an improved R-hat, *Bayesian Analysis*
