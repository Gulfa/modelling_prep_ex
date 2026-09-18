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

Vaccination is an **initial coverage** only: a fraction `vacc_coverage` of each group starts in $S_{g1}$. Nobody changes vaccination status during a run.

## Force of infection

$$
\lambda_g(t) = \beta \sqrt{m_g(t)} \sum_{h} C_{gh}\, \sqrt{m_h(t)}\; \frac{\text{trans}_{h0} I_{h0} + \text{trans}_{h1} I_{h1}}{N_h}
$$

- $N_h$ is the living population of group $h$ ($S+E+I+H+R$, both strata).
- $C_{gh}$ is the number of contacts per day a person in group $g$ has with people in group $h$.
- $\text{sus}_{gv}$ and $\text{trans}_{hv}$ are the relative susceptibility and infectiousness.
- $m_g(t)$ is the contact multiplier from `beta_schedule`. Contacts between $g$ and $h$ are scaled by $\sqrt{m_g m_h}$:
  - A multiplier applied to everyone scales all transmission by exactly that factor.
  - A multiplier applied to one region reduces its internal contacts by $m$ and its contacts with other regions by $\sqrt{m}$.

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
| S → E | $\text{Bin}(S, 1 - e^{-\lambda\,\text{sus}\,dt})$ |
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

## Outputs

All outputs are recorded daily and stratified by age, region and vaccination status. [`get_output()`](../R/output.R) sums over any strata you don't ask for.

| Variable | Meaning |
|---|---|
| `S`, `E`, `I`, `H`, `R` | number in each compartment (`H` = hospital occupancy) |
| `D_hosp`, `D_comm`, `D` | cumulative deaths (hospital, community, both) |
| `cum_infections`, `cum_hosp` | cumulative infections, hospital admissions |
| `inc_infections`, `inc_hosp`, `inc_deaths_hosp`, `inc_deaths_comm`, `inc_deaths` | new infections / admissions / deaths during day $(t-1, t]$ |

Initial infections (`I0`, `E0`) are not counted in `cum_infections` or `inc_infections`.
