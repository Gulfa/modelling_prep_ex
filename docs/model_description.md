# Model description

`seirvax` contains one model in two versions: a deterministic one (ODEs) and a stochastic one (discrete-time binomial draws). Both versions share the same parameters and outputs.

## Structure

The population is divided into **groups** $g = 1,\dots,n$ and two **vaccination strata** $v \in \{0, 1\}$ (unvaccinated, vaccinated). In the R interface, groups are the combinations of age group $a$ and region $r$. The model code only sees $n = n_\text{age} \cdot n_\text{region}$ groups, numbered with age varying fastest: $g = (r-1)\,n_\text{age} + a$.

Every group and stratum has the compartments

```
            λ_g·sus_gv             σ                γ·p_hosp_gv                  (1 - p_death_hosp_gv) / dur_hosp_gv
  S_gv ───────────────► E_gv ──────────► I_gv ───────────────────► H_gv ────────────────────────────────► R_gv
                                           │                          │
                                           │                          └─ p_death_hosp_gv / dur_hosp_gv ──► D_hosp_gv
                                           ├─ γ·(1-p_hosp_gv)·p_death_comm_gv ──► D_comm_gv
                                           └─ γ·(1-p_hosp_gv)·(1-p_death_comm_gv) ──► R_gv
```

- **S:** susceptible.
- **E:** exposed (latent, not infectious). Mean duration is the latent period, $1/\sigma$.
- **I:** infectious. Mean duration is the infectious period, $1/\gamma$. When people leave I, they are hospitalised with probability `p_hosp`. Otherwise they die in the community with probability `p_death_comm` or recover.
- **H:** in hospital, not infectious. Mean length of stay is `dur_hosp`. On leaving, patients die with probability `p_death_hosp` or recover.
- **R:** recovered and immune. There is no waning.
- **D_hosp, D_comm:** cumulative deaths in hospital and in the community.

Vaccination is an **initial coverage**, `vacc_coverage`: that fraction of each group starts in $S_{g1}$. Optionally, **`vacc_schedule`** gives doses during the run as well (see [Vaccination over time](#vaccination-over-time)). Nobody loses vaccinated status.

## Force of infection

$$
\lambda_g(t) = \beta \sum_{h} C_{gh}\, \min\big(m_g(t), m_h(t)\big)\; \frac{\text{trans}_{h0} I_{h0} + \text{trans}_{h1} I_{h1}}{N_h}
$$

- $N_h$ is the living population of group $h$ ($S+E+I+H+R$, both strata).
- $C_{gh}$ is the number of contacts per day a person in group $g$ has with people in group $h$.
- $\text{sus}_{gv}$ and $\text{trans}_{hv}$ are the relative susceptibility and infectiousness.
- $m_g(t)$ is the contact multiplier from `beta_schedule`: the fraction of group $g$'s contacts that still happen. Contacts between $g$ and $h$ are scaled by $\min(m_g, m_h)$, so the more restricted group decides:
  - A multiplier applied to everyone scales all transmission by exactly that factor.
  - "65+ at 0.5" halves every contact involving someone aged 65+, including their contacts with other age groups.
  - "Region A at 0.25" cuts all contacts involving region A residents to 25%, including their contacts with other regions.
  - Contact rates stay reciprocal, because the factor is the same from both sides.
  - Where two restricted groups meet, the stricter multiplier applies; the two are not multiplied together.

New infections in group $g$, stratum $v$, occur at rate $\lambda_g\,\text{sus}_{gv}\,S_{gv}$.

### Combining age and region mixing

$$
C_{(a,r),(b,s)} = C^\text{age}_{ab} \cdot M^\text{region}_{rs}
\qquad\Longleftrightarrow\qquad
C = M^\text{region} \otimes C^\text{age}
$$

- `C_age[a, b]` is contacts per day an age-$a$ person has with age-$b$ people.
- `M_region[r, s]` is the fraction of a region-$r$ resident's contacts made with residents of region $s$. Rows sum to 1.

This assumes that everyone has the same age-mixing pattern and that the region matrix only redistributes where contacts happen. A 1×1 matrix removes the corresponding dimension.

### Mixing that changes over time

`mixing_schedule` switches the matrices at given times, so the *pattern* of contacts can change, not just its level:

```r
build_model(..., C_age = C_term, M_region = M_normal,
  mixing_schedule = list(
    list(time = 40, C_age = C_holiday),                  # schools close
    list(time = 40, M_region = M_restricted),            # ... and travel is restricted
    list(time = 70, C_age = C_term, M_region = M_normal) # both back to normal
  ))
```

Each entry replaces `C_age`, `M_region` or both from its time onwards; whichever is not given carries over. `C_age` and `M_region` themselves apply from time 0. Internally this is one `n × n × n_time` array, interpolated piecewise-constant, so the cost is a little memory and nothing at run time.

`beta_schedule` and `mixing_schedule` are independent and compose: the multipliers scale whichever matrix is in force. Use `mixing_schedule` when the *shape* changes (school holidays, closed borders) and `beta_schedule` when the *level* does (distancing, a group reducing its contacts).

**R0 is calibrated from the mixing at time 0.** `R_eff(model, time)` uses the matrix in force at `time`, so it reflects later changes.

## Vaccination over time

`vacc_schedule` gives doses per day from given times, for everyone or for chosen age groups and regions:

```r
build_model(...,
  vacc_schedule = data.frame(time  = c(40, 70),
                             doses = c(5000, 0),     # 5,000 a day from day 40 to day 70
                             age   = "65+"))
```

- Rows follow the same rules as `beta_schedule`: from `time` onwards, the groups a row matches (`NA`, or a missing `age`/`region` column, matches all) receive `doses` per day between them, shared in proportion to their population, until a later row matching them overrides it. Before the first row, no doses are given. `doses = 0` ends a campaign.
- Doses go to **unvaccinated people who are not ill**: $S_{g0}$ and $R_{g0}$, in proportion to their size. They move to $S_{g1}$ and $R_{g1}$ respectively, so a dose given to someone who has already recovered changes nothing, as in reality. People in E, I or H are not vaccinated.
- **Protection is immediate**: a vaccinated person has the vaccinated parameters from the moment of the dose. To approximate a delay until protection, start the schedule correspondingly later.
- In group $g$, the doses $d_g(t)$ are given at per-capita rate $\min\big(d_g(t) / (S_{g0} + R_{g0}),\ 10\big)$ per day among the eligible. This delivers exactly $d_g$ a day while enough eligible people remain; when they run out, the cap makes vaccination slow and stop instead of driving the compartments negative.
- Doses given are recorded as `cum_doses` and `inc_doses`, under `vacc_status = "vaccinated"`. The initial coverage is not counted.
- `R_eff(model, time)` uses the initial coverage only, not the doses given since.

Because doses change the dynamics from their start time onwards, unlike `vacc_coverage`, a vaccination schedule *can* be applied to a projection that continues from a fit (see [Uncertainty in projections](#uncertainty-in-projections)).

## Vaccine effects

The vaccinated stratum uses the unvaccinated parameters multiplied by one minus the vaccine effectiveness (VE). Each VE can vary by age and region.

| Vaccinated parameter | = unvaccinated × |
|---|---|
| `sus` | `1 - ve_infection` |
| `trans` | `1 - ve_transmission` |
| `p_hosp` | `1 - ve_hosp` |
| `p_death_hosp`, `p_death_comm` | `1 - ve_death` |
| `dur_hosp` | replaced by `dur_hosp_vacc` if given |

`ve_hosp` and `ve_death` are protection **given infection**. The total protection against hospitalisation is $1-(1-\text{ve\_infection})(1-\text{ve\_hosp})$.

## Calibrating β from R0

When `R0` is given, β is chosen so that the spectral radius of the next-generation matrix is `R0`. That radius is computed for a fully susceptible, unvaccinated population with all contact multipliers set to 1:

$$
K_{gh} = \beta\, \text{sus}_{g0}\, C_{gh}\, \text{trans}_{h0}\, \frac{N_g}{N_h} \cdot \frac{1}{\gamma}
$$

`R_eff(model, time)` uses the same construction over both strata, with the initial susceptible populations and the multipliers at `time`. `print(model)` shows both values.

**Known limitation:** the next-generation R0 describes generation-to-generation growth once infections follow the dominant eigenvector's age and region distribution. If the epidemic is seeded differently (for example in one region or one age group), early growth can be faster or slower than R0 suggests until that distribution is reached. The same applies to the calibrated β.

## Stochastic version

The stochastic model uses tau-leaping with a fixed time step `dt`. Each step:

| Transition | Draw |
|---|---|
| leaving S | $n_S = \text{Bin}(S, 1 - e^{-(\lambda\,\text{sus} + v)\,dt})$, with $v$ the vaccination rate (0 in the vaccinated stratum) |
| S → E | $\text{Bin}(n_S, \lambda\,\text{sus} / (\lambda\,\text{sus} + v))$; the rest of $n_S$ is vaccinated |
| R vaccinated | $\text{Bin}(R, 1 - e^{-v\,dt})$ |
| E → I | $\text{Bin}(E, 1 - e^{-\sigma dt})$ |
| leaving I | $n_I = \text{Bin}(I, 1 - e^{-\gamma dt})$ |
| I → H | $\text{Bin}(n_I, p_\text{hosp})$ |
| I → D_comm | $\text{Bin}(n_I - n_{IH}, p_\text{death,comm})$ |
| leaving H | $n_H = \text{Bin}(H, 1 - e^{-dt/\text{dur\_hosp}})$ |
| H → D_hosp | $\text{Bin}(n_H, p_\text{death,hosp})$ |

The chained binomials guarantee that competing outcomes never exceed the number of people leaving a compartment, and counts stay whole numbers.

**Time-step bias.** In discrete time, a person stays in each stage for $dt/(1-e^{-\text{rate}\cdot dt}) \approx 1/\text{rate} + dt/2$ on average. The infectious period is therefore slightly too long, and epidemics are slightly larger than in the ODE model. For R0 = 2 the mean final size is:

| `dt` | final size vs. ODE |
|---|---|
| 1 | +6.1% |
| 0.5 | +3.2% |
| 0.25 | +1.7% |
| 0.1 (default) | +0.7% |
| 0.05 | +0.4% |

## Plot styling

Plots use ColorBrewer **Dark2** for groups (age, region, scenario, MCMC chain). With more than eight groups the palette is interpolated rather than recycled, so no two groups share a colour. Spaghetti lines fade as they get more numerous but never below 30% opacity; `alpha` overrides that on `plot()`. Maps keep viridis, which suits a continuous scale.

## Outputs

All outputs are recorded daily and stratified by age, region and vaccination status. [`get_output()`](../R/output.R) sums over any strata you don't ask for.

| Variable | Meaning |
|---|---|
| `S`, `E`, `I`, `H`, `R` | number in each compartment (`H` = hospital occupancy) |
| `D_hosp`, `D_comm`, `D` | cumulative deaths (hospital, community, both) |
| `cum_infections`, `cum_hosp`, `cum_doses` | cumulative infections, hospital admissions, doses given by `vacc_schedule` |
| `inc_infections`, `inc_hosp`, `inc_deaths_hosp`, `inc_deaths_comm`, `inc_deaths`, `inc_doses` | new infections / admissions / deaths / doses during day $(t-1, t]$ |

Initial infections (`I0`, `E0`) are not counted in `cum_infections` or `inc_infections`.

## Fitting to hospital admissions

`fit_model()` estimates `R0`, and optionally `hosp_scale` and `nb_size` (see `fit_parameters()`), from observed daily hospital admissions. All other parameters stay fixed at the values in the model you pass in.

**Observation model.** Data can be one national series, or several series by region and/or age. For series $k$ on day $t$, the modelled admissions are

$$
\mu_k(t) = \sum_g \text{obs\_map}_{kg} \sum_v \text{inc\_hosp}_{gv}(t),
$$

where `obs_map[k, g]` is 1 if group $g$ belongs to series $k$. Observed admissions are assumed to be $\text{NegBin}(\text{size} = \text{nb\_size},\ \text{mean} = \mu_k(t))$. A smaller `nb_size` means noisier data; the variance is $\mu + \mu^2/\text{nb\_size}$. Missing values (`NA`) are skipped. The comparison is at the end of both odin files, and when the model is only simulated it's ignored.

**Fitted parameters.**
- **`R0`:** β is recalibrated from R0 with the next-generation matrix, exactly as in `build_model()`.
- **`hosp_scale`:** multiplies `p_hosp` in every group and in both vaccination strata, capped at 1.
- **`nb_size`:** the dispersion above.

With data from the early exponential phase only, the growth rate (and so R0) is well identified. The *level* of admissions is not: a larger `hosp_scale` can be traded against fewer or later initial infections. Fit `hosp_scale` only when the data covers more of the epidemic, or when the seeding is well known.

**Methods.**

| `method` | Likelihood | Algorithm | Uncertainty |
|---|---|---|---|
| `"optim"` | deterministic model | `stats::optim` (L-BFGS-B) to the maximum a posteriori | normal approximation from the curvature at the optimum |
| `"mcmc"` | deterministic model | monty's adaptive random-walk sampler, started at the `optim` estimate and proposal | posterior samples |
| `"pmcmc"` | stochastic model: particle filter with `n_particles` particles | monty random walk, proposal from the `optim` curvature | posterior samples, including the epidemic's own randomness |

For what particle MCMC is doing — the filter, why the likelihood it returns is noisy, and how to choose the number of particles — see Endo, van Leeuwen and Baguelin, [*Introduction to particle Markov-chain Monte Carlo for disease dynamics modellers*](https://doi.org/10.1016/j.epidem.2019.100363), Epidemics 29 (2019), which works through the method with R code.

With default uniform priors, the `optim` estimate is the maximum likelihood estimate within the prior range.

**Time alignment.** Day 0 of the data is time 0 of the model, and `admissions` on day `t` are admissions during day $(t-1, t]$. If the epidemic started before the data, choose the model's seeding (`I0`) to match. Fitting the start date isn't supported yet.

## Uncertainty in projections

`project_fit()` and `plot(fit)` propagate parameter uncertainty by Monte Carlo: draw parameters, rebuild the model for each draw, run it, and take quantiles across draws at each time point.

- **Where the draws come from.** For `"mcmc"` and `"pmcmc"`, a random subsample (without replacement) of the pooled post-burn-in samples. For `"optim"`, draws from a multivariate normal at the estimate with the covariance implied by the curvature there, rejected until they fall inside the prior's range. Draws are whole parameter vectors, so correlations between parameters are kept.
- **Each draw rebuilds the model**, so a drawn `R0` is recalibrated to its own β.
- **The projection starts where the data ends.** `project_fit()` reports from `from`, which defaults to the last day of the data, so a projection is about the future rather than a re-prediction of the fitted period. `from = 0` gives the whole epidemic; `plot(fit)` uses that to show data and model together.
- **For `pmcmc` fits the start is the filtered state.** Each draw runs the particle filter over the data up to `from` and continues from a state sampled from it, so the projection carries on the epidemic that was actually observed — including where it happened to be above or below the model's average. Conditioning this way narrows the projection: in one check, cumulative admissions after 90 further days spanned 910–11,298 conditioned against 365–13,383 unconditioned. The deterministic methods have no such state and simply run from time 0.
- The result behaves like a model run: each draw is one of its simulations, so `get_output()`, `summary()` and `plot()` work on it directly, and `$draws` holds the parameter values behind each draw. `seed` makes the choice reproducible.

What the band does **not** include:

| | |
|---|---|
| Process noise | only with `project_fit(type = "stochastic", n_sim = ...)`; the default is deterministic, one trajectory per draw |
| Observation noise | only in `plot(fit)`, as the lighter outer band |
| Parameters that were not fitted | fixed at the model's values — severity, seeding, contact matrices, durations, population |
| Model structure | not represented at all |

`project_summary()` applies the same machinery to end-of-projection outcomes, broken down by age group and/or region and summarised over draws. The measures are named after the output variables they come from — `cum_infections`, `cum_hosp` (admissions) and `D` (deaths) at the end of the projection, plus `peak_hosp_occupancy` and `peak_hosp_day` from `H` — so `summary()`, `project_summary()` and `get_output()` all use one vocabulary. (The fitting *data* column is `admissions`, since that is a column of observations, not a model output.) It also reports each stratum's `share` of the national total, computed within each draw, which is the honest way to ask "what fraction of admissions are in 65+?" Shares are given for the totals only — strata peak at different times, so their peaks do not add up to the national peak.

A caveat specific to breakdowns: if only `R0` is fitted, the *shape* of the age and region distribution is fixed by the contact matrix and the severity parameters, so its interval will be very narrow however uncertain the totals are. That narrowness reflects the assumptions, not evidence from the data.

**Scenarios on top of a fit.** `project_fit()` and `project_summary()` take a `model` argument, so the fitted parameters can be applied to a modified model while keeping the posterior draws:

```r
project_summary(fit, model = update_model(fit$base_model,
                  beta_schedule = data.frame(time = 100, mult = 0.6)),
                t_max = 365, measures = "cum_hosp")
```

The alternative, `update_model(fit$model, ...)`, gives a single trajectory at the point estimate and throws the uncertainty away. One caveat: a drawn `R0` is recalibrated for whatever mixing the scenario has, so changing a contact matrix holds R0 fixed rather than β.

One more thing to keep in mind: the bands are **pointwise**. The median curve is the median at each day, not any single draw, and with spread-out peak timing it is flatter than every individual trajectory.

Note that a scenario starting from a filtered state cannot change the past: `vacc_coverage`, `I0` and anything else that only sets initial conditions has no effect, because the state already encodes them. Only parameters governing the dynamics after `from` matter — `beta_schedule`, `mixing_schedule` and `vacc_schedule` entries from `from` onwards. So vaccinating from the end of the data is a `vacc_schedule` starting at `from`, not a change to `vacc_coverage`.
