# Model equations

## Variable definitions

### Equation subscripts

The following letters are used for subscripts to identify the dimension
and indexing for model variables.

| Subscript | Definition            | Number |
|:----------|:----------------------|:-------|
| $`y`$     | Year                  | 1.1    |
| $`m`$     | Season                | 1.2    |
| $`a`$     | Age                   | 1.3    |
| $`r`$     | Region (spatial area) | 1.4    |
| $`f`$     | Fleet                 | 1.5    |
| $`s`$     | Stock                 | 1.6    |
| $`i`$     | Index of abundance    | 1.7    |
| $`\ell`$  | Length bin            | 1.8    |

### Fixed parameters

The parameters here are set up as fixed inputs prior to fitting the
model.

| Description | Symbol | Number |
|:---|:---|:---|
| Stock weight at age | $`w_{y,m,a,s}`$ | 2.1 |
| Fishery weight at age | $`w^F_{y,m,a,f,s}`$ | 2.2 |
| Maturity ogive | $`m_{y,a,s}`$ | 2.3 |
| Fecundity | $`f_{y,a,s}`$ | 2.4 |
| Natural mortality (instantaneous per year) | $`M_{y,a,s}`$ | 2.5 |
| Length-at-age probability | $`\textrm{Pr}(\ell \mid a )_{y,m,s}`$ | 2.6 |
| Fractional parameter (between 0 - 1) | $`\Delta`$ | 2.7 |
| Season length (relative to year) | $`\Delta_m`$ | 2.8 |
| Spawn timing (within season) | $`\Delta_{sp}`$ | 2.9 |
| Season of spawning (subscript) | $`m_\textrm{sp}`$ | 2.10 |
| Season of recruitment (subscript) | $`m_\textrm{rec}`$ | 2.11 |
| Natal spawning (proportion of mature individuals that spawn) | $`n_{r,s}`$ | 2.12 |
| Stock scaling parameter | $`r_s`$ | 2.13 |
| Survey timing | $`\Delta_i`$ | 2.14 |
| Survey sampling coverage | $`\delta_{i,r,s}`$ | 2.15 |

### Estimated parameters

The parameters here are set up to be either estimated or fixed in the
model. Such parameters are identified as $`x`$ which is estimated over
all real numbers and transformed to the appropriate model parameter
described below.

- Unfished recruitment is scaled by an additional user parameter $`r_s`$
  which is intended to aid convergence. For multi-stock models, $`r_s`$
  should be proportional to the expected stock size, i.e., large values
  for large stocks.
- Maturity and natural mortality can either be estimated or defined by
  the user as above.
- Fishing mortality is fixed to zero when the corresponding catch is
  less than $`10^{-8}`$.

| Description | Symbol | Number |
|:---|:---|:---|
| Estimated parameter (all real numbers) | $`x`$ | 3.1 |
| Unfished recruitment | $`R_{0,s} = r_s \exp(x^{R0}_s)`$ | 3.2 |
| Beverton-Holt stock-recruit steepness | $`h_s = 0.8 \times \textrm{logit}^{-1}(x^h_s) + 0.2`$ | 3.3 |
| Ricker stock-recruit steepness | $`h_s = \exp(x^h_s) + 0.2`$ | 3.4 |
| Age of 50% maturity ($`A`$ is the maximum age) | $`a^{50}_s=A\times\textrm{logit}^{-1}(x^{a50}_s)`$ | 3.5 |
| Age of 95% maturity | $`a^{95}_s = a^{50}_s + \exp(x^{a95}_s)`$ | 3.6 |
| Natural mortality | $`M_s = \exp(x^M_s)`$ | 3.7 |
| Recruitment deviates | $`x^R_{y,s}`$ | 3.8 |
| Standard deviation of recruitment deviates | $`\sigma^R_s=\exp(x^{\sigma_R}_s)`$ | 3.9 |
| Fleet catchability by stock | $`q^F_{f,s}= \exp(x^{qF}_{f,s})`$ | 3.10 |
| Fishing mortality | $`F_{y,m,f,r} = \begin{cases}  \exp(x^{\textrm{Fmult}}_f) \quad & y = y_{\textrm{ref}}, m = m_{\textrm{ref}}, r = r_{\textrm{ref}}\\ \exp(x^{\textrm{Fmult}}_f + x^{\textrm{Fdev}}_{y,m,r}) \quad & \textrm{otherwise} \end{cases}`$ | 3.11 |
| Selectivity - length of full selectivity ($`f`$ and $`i`$ subscripts are interchangeable) | $`\mu_f = L_{\textrm{max}} \times \textrm{logit}^{-1}(x^{\mu}_f)`$ | 3.12 |
| Selectivity - width of ascending limb | $`\sigma^{\textrm{asc}}_f = \exp(x^{\textrm{asc}}_f)`$ | 3.13 |
| Selectivity - width of descending limb | $`\sigma^{\textrm{dsc}}_f = \exp(x^{\textrm{dsc}}_f)`$ | 3.14 |
| Base movement parameters from origin $`r`$ to destination $`r'`$ | $`x^b_{m,a,r,r',s}`$ | 3.15 |
| Attractivity term - movement | $`x^g_{y,m,a,r',s}`$ | 3.16 |
| Viscosity term - movement | $`x^v_{y,m,a,s}`$ | 3.17 |
| Initial equilibrium fishing mortality | $`F^{\textrm{eq}}_{m,f,r}=\exp(x^{\textrm{Feq}}_{m,f,r})`$ | 3.18 |
| Deviations from the equilibrium age structure | $`x^{\textrm{Req}}_{a,s}`$ | 3.19 |

### Derived variables

This section defines additional variables derived from data or estimated
parameters described in the previous sections.

- Selectivity is reported here in terms of length. The corresponding
  age-based selectivity by stock is obtained from the length-at-age
  probability key and is seasonally-varying based on the growth pattern.
- Movement is parameterized with three arrays and several configurations
  are possible.
- Stock-recruit functions use the steepness parameterization, along with
  the unfished recruitment and the unfished spawning output per recruit
  ($`\phi_0`$). In seasonal and multi-region models, the population
  dynamics model is used to numerically obtain $`\phi_0`$ by setting
  $`F_{y,m,f,r} = 0`$, recruitment to 1, and all other parameters to
  constant seasonal values. $`\phi_0`$ is the equilibrium spawning
  output at the end of this numerical spool-up.

| Description | Equation | Number |
|:---|:---|:---|
| Selectivity at length | $`s_{\ell,f} = \begin{cases} \exp\left(-0.5\left[\dfrac{L_\ell - \mu_f}{\sigma^{\textrm{asc}}_f}\right]^2\right) &, L_\ell < \mu_f \\ \exp\left(-0.5\left[\dfrac{L_\ell - \mu_f}{\sigma^{\textrm{dsc}}_f}\right]^2\right) &, L_\ell \ge \mu_f \end{cases}`$ | 4.1 |
| Selectivity at age | $`s_{y,m,a,f,s} = \sum_\ell \textrm{Pr}(\ell \mid a )_{y,m,s} \times s_{\ell,f}`$ | 4.2 |
| Maturity (if estimated) | $`m_{a,s} = \left[1 + \exp\left(-\log(19)\frac{a-a^{50}_s}{a^{95}_s - a^{50}_s}\right)\right]^{-1}`$ | 4.3 |
| Movement from $`r`$ to $`r'`$ | $`\textrm{mov}_{y,m,a,r,r',s} = \dfrac{\exp(x^b_{m,a,r,r',s} + x^g_{y,m,a,r,s} + x^v_{y,m,a,s})}{\sum_{r'} \exp(x^b_{m,a,r,r',s} + x^g_{y,m,a,r,s} + x^v_{y,m,a,s})}`$ | 4.4 |
| Beverton-Holt stock recruit parameters | $`\begin{align} \alpha_s &= \dfrac{4h_s}{(1-h_s)\phi_{0,s}}\\ \beta_s &= \dfrac{5h_s-1}{(1-h_s)R_{0,s}\phi_{0,s}}\end{align}`$ | 4.5 |
| Ricker stock recruit parameters | $`\begin{align} \alpha_s &= \dfrac{5h^{1.25}}{\phi_{0,s}}\\ \beta_s &= \dfrac{\log([5h_s]^{1.25})}{R_{0,s}\phi_{0,s}}\end{align}`$ | 4.6 |

## Population dynamics

The following equations project the population forward in time.

- To obtain the initial abundance $`N_{y=1,m=1,a,r,s}`$ array in
  seasonal and multi-region models, a numerical spool-up is performed
  with the seasonal fishing mortality equal to
  $`F^{\textrm{eq}}_{m,f,r}`$, recruitment to 1, and all other
  parameters set to constant seasonal values from the first year of the
  model. From this initialization, the equilibrium spawners per recruit
  $`\phi_{eq}`$ is the final spawning output, and the seasonal numbers
  per recruit $`\textrm{NPR}^{\textrm{eq}}_{m,a,r,s}`$ is obtained from
  the abundance array. The initial abundance is the product of the
  equilibrium recruitment and numbers per recruit.
- It is possible that some proportion of the mature population do not
  contribute to the annual spawning based on the natal spawning
  parameter specifying the spatial spawning pattern. Thus, there is a
  distinction between potential spawners and realized spawners. The
  unfished replacement line of the stock-recruit relationship
  ($`R_s = 1/\phi_0`$) is based on the realized spawning in equilibrium.

| Description | Equation | Number |
|:---|:---|:---|
| Index of fishing effort (instantaneous per season) | $`F_{y,m,f,r}`$ | 5.1 |
| Stock abundance | $`N_{y,m,a,r,s}`$ | 5.2 |
| Equilibrium recruitment | $`R_{\textrm{eq},s} = \begin{cases} \dfrac{\alpha_s \phi_{eq,s} - 1}{\beta_s\phi_{eq,s}} & \textrm{Beverton-Holt}\\ \dfrac{\log(\alpha_s \phi_{eq,s})}{\beta_s\phi_{eq,s}} & \textrm{Ricker} \end{cases}`$ | 5.3 |
| Initial abundance | $`N_{y=1,m=1,a,r,s} = R_{\textrm{eq},s}\exp(x^{\textrm{Req}}_{a,s}) \times \textrm{NPR}^{\textrm{eq}}_{m=1,a,r,s}`$ | 5.4 |
| Fishing mortality by time, age, fleet, region, stock | $`F_{y,m,a,f,r,s} = q^F_{f,s}s_{y,m,a,f,s}F_{y,m,f,r}`$ | 5.5 |
| Total mortality (instantaneous per season) | $`Z_{y,m,a,r,s} = \Delta_m M_{y,a,s} + \sum_f F_{y,m,a,f,r,s}`$ | 5.6 |
| Seasonal abundance without incoming recruitment (after survival and movement) | $`N_{y,m,a,r',s} = \sum_r N_{y,m-1,a,r,s}\exp(-Z_{y,m-1,a,r,s}) \textrm{mov}_{y,m-1,a,r,r',s}`$ | 5.7 |
| Potential spawners (PS) | $`N^{\textrm{PS}}_{y,a,r,s} = N_{y,m = m_\textrm{sp},a,r,s}\exp(-\Delta_\textrm{sp}Z_{y,m=m_\textrm{sp},a,r,s})m_{y,a,s}`$ | 5.8 |
| Active spawners | $`N^{\textrm{S}}_{y,a,r,s} = N^{\textrm{PS}}_{y,m = m_{\textrm{sp}},a,r,s}\times n_{r,s}`$ | 5.9 |
| Spawning output | $`S_{y,r,s} = \sum_a N^{\textrm{S}}_{y,a,r,s} f_{y,a,s}`$ | 5.10 |
| Recruitment: Beverton-Holt | $`R_{y,s} = \dfrac{\alpha_s \sum_r S_{y,r,s}}{1 + \beta_s \sum_r S_{y,r,s}}\exp(x^R_{y,s} - 0.5(\sigma^R_s)^2)`$ | 5.11 |
| Recruitment: Ricker | $`R_{y,s} = \alpha_s \sum_r S_{y,r,s} \exp(-\beta_s \times \sum_r S_{y,r,s})\exp(x^R_{y,s} - 0.5(\sigma^R_s)^2)`$ | 5.12 |
| Seasonal abundance (incoming recruitment and advancing age class when $`m = m_{\textrm{rec}}`$) | $`N_{y,m,a,r',s} = \begin{cases} R_{y,s}& a = 0\\ \sum_r N_{y,m-1,a-1,r,s}\exp(-Z_{y,m-1,a-1,r,s}) \textrm{mov}_{y,m-1,a-1,r,r',s} & a = 1, \ldots, A-1\\ \sum_{a'=A-1}^A\sum_r N_{y,m-1,a',r,s}\exp(-Z_{y,m-1,a',r,s})\textrm{mov}_{y,m-1,a',r,r',s} & a = A\end{cases}`$ | 5.13 |

## Report variables

Here, we calculate additional variables that are not needed for the
population dynamics model, but are of interest for fitting the model or
for reporting.

- In a multi-region and/or seasonal model, we may want a summary fishing
  mortality (per year) across all regions and fleets ($`F_{y,a,s}`$)
  which calculated from the Baranov equation with natural mortality
  $`M_{y,a,s}`$, total stock abundance at the beginning of the year
  $`N_{y,a,s}`$, and total catch $`C^N_{y,a,s}`$. The summary total
  mortality (per year) is then $`Z_{y,a,s} = F_{y,a,s} + M_{y,a,s}`$.
- Vulnerable biomass is the availability of the stock to individual
  fleets.
- When fitting to close-kin genetic data, we can calculate the
  probability of parent-offspring pairs (POP) with the cohort year of
  the offspring is $`y`$, the parental age at capture is $`a'`$, and the
  capture year of the parent $`t`$.
- The half-offspring pair probability is calculated from the parental
  probability in years $`i`$ and $`j`$, the cohort year of the older and
  younger sibling, respectively, and the parental survival from year
  $`i`$ to year $`j`$. The parental age is not observed, so we calculate
  the probability across all potential ages and follow each cohort from
  $`i`$ to $`j`$.

| Description | Equation | Number |
|:---|:---|:---|
| Equilibrium catch (abundance, age) | $`C^{Neq}_{m,a,f,r,s} = \dfrac{F^{\textrm{eq}}_{m,a,f,r,s}}{Z^{\textrm{eq}}_{m,a,r,s}} (1 - \exp(-Z^{\textrm{eq}}_{y,m,a,r,s})) N^{\textrm{eq}}_{m,a,r,s}`$ | 6.1 |
| Equilibrium catch (biomass) | $`C^{Beq}_{m,f,r,s} = \sum_a w^F_{y=1,m,a,f,s} C^{Neq}_{m,a,f,r,s}`$ | 6.2 |
| Catch (abundance, age) | $`C^N_{y,m,a,f,r,s} = \dfrac{F_{y,m,a,f,r,s}}{Z_{y,m,a,r,s}} (1 - \exp(-Z_{y,m,a,r,s})) N_{y,m,a,r,s}`$ | 6.3 |
| Catch (abundance, length) | $`C^N_{y,m,l,f,r,s} = \sum_a C^N_{y,m,a,f,r,s} \times s_{\ell,f} \textrm{Pr}(\ell\mid a)_{y,m,s}`$ | 6.4 |
| Catch (biomass) | $`C^B_{y,m,f,r,s} = \sum_a w^F_{y,m,a,f,s} C^N_{y,m,a,f,r,s}`$ | 6.5 |
| Total biomass | $`B_{y,m,r,s} = \sum_a w_{y,m,a,s} N_{y,m,a,r,s}`$ | 6.6 |
| Vulnerable biomass | $`V_{y,m,r,s} = \sum_a s_{y,m,a,f,s} w^F_{y,m,a,f,s} N_{y,m,a,r,s}`$ | 6.7 |
| Index age composition | $`I^A_{y,m,a,i,s} = q_i \sum_r s_{y,m,a,i,s} N_{y,m,a,r,s} \exp(-\Delta_i \times Z_{y,m,a,r,s})\delta_{i,r,s}`$ | 6.8 |
| Index length composition | $`I^L_{y,m,\ell,i,s} = \sum_a I^A_{y,m,a,i,s} \textrm{Pr}(\ell\mid a)_{y,m,s}`$ | 6.9 |
| Biomass index | $`I_{y,m,i} = q_i \sum_a \sum_s w_{y,m,a,s} I^A_{y,m,a,i,s}`$ | 6.10 |
| Annual catch at age (all seasons, fleets, and regions) | $`C^N_{y,a,s} = \sum_m\sum_f\sum_rC^N_{y,m,a,f,r,s}`$ | 6.11 |
| Abundance at year start across regions | $`N_{y,a,s} = \sum_r N_{y,m=1,a,r,s}`$ | 6.12 |
| Summary fishing mortality per year $`F_{y,a,s}`$ | $`C^N_{y,a,s} = \dfrac{F_{y,a,s}}{F_{y,a,s} + M_{y,a,s}}(1 - \exp(-[F_{y,a,s} + M_{y,a,s}])) N_{y,a,s}`$ | 6.13 |
| Parent-offspring probability | $`\textrm{Pr}(\textrm{POP} \mid y, t, a',s) = 2 \times \dfrac{f_{y,a = y - (t - a'),s}}{\sum_a f_{y,a,s}N^{\textrm{S}}_{y,a,s}}`$ | 6.14 |
| Half-sibling probability | $`\textrm{Pr}(\textrm{HSP} \mid i,j,s) = 4 \times \sum_a c_{i,a,s} \times \textrm{surv}_{ijs} \times c_{j,a,s}`$ | 6.15 |
| Parental probability for older sibling $`i`$ | $`c_{i,a,s} = \dfrac{N_{y=i,a,s}f_{y=i,a,s}}{\sum_{a'}N_{y=i,a',s}f_{y=i,a',s}}`$ | 6.16 |
| Parental survival from year $`i`$ to $`j`$ | $`\textrm{surv}_{ijs} = \exp\left(-\sum\limits_{t=0}^{j-i-1}Z_{i+t,a+t,s}\right)`$ | 6.17 |
| Parental probability for younger sibling $`j`$ | $`c_{j,a,s} = \dfrac{f_{y=j,a,s}}{\sum_{a'}N_{y=j,a',s}f_{y=j,a',s}}`$ | 6.18 |

## Objective function

The objective function is the sum of the negative log-likelihoods,
negative log-priors, and penalty function.

### Likelihoods

The statistical distributions used for the likelihoods of the data are
described. The dimensions of the data are given below as well as the
corresponding model variable for the predicted value, which is typically
summed across stocks, (except for stock composition). Composition data
are presented as proportions $`p`$ and a separate table provides the
mean and variance of the various likelihood options.

The close-kin likelihood uses the ratio of matches (for either
parent-offspring or sibling matches) and the number of pairwise
comparisons ($`N`$).

| Data | Symbol | Predicted | Distribution | Variance | Number |
|:---|:---|:---|:---|:---|:---|
| Equilibrium catch (biomass) | $`C^{Beq}_{m,f,r}`$ | $`\sum_s \hat{C}^{Beq}_{m,f,r,s}`$ | Lognormal | $`\sigma = 0.01`$ | 7.1 |
| Catch (biomass) | $`C^B_{y,m,f,r}`$ | $`\sum_s \hat{C}^B_{y,m,f,r,s}`$ | Lognormal | $`\sigma^C_{y,m,f,r}`$ | 7.2 |
| Catch at age | $`p^{CN}_{y,m,a,f,r}`$ | $`\sum_s \hat{C}^N_{y,m,a,f,r,s}/\sum_a\sum_s \hat{C}^N_{y,m,a,f,r,s}`$ | Composition | See next table | 7.3 |
| Catch at length | $`p^{CN}_{y,m,l,f,r}`$ | $`\sum_s \hat{C}^N_{y,m,\ell,f,r,s}/\sum_\ell\sum_s \hat{C}^N_{y,m,\ell,f,r,s}`$ | Composition | See next table | 7.4 |
| Total indices | $`I_{y,m,i}`$ | $`\hat{I}_{y,m,i}`$ | Lognormal | $`\sigma^I_{y,m,i}`$ | 7.5 |
| Index at age | $`p^{IN}_{y,m,a,i}`$ | $`\sum_s \hat{I}^N_{y,m,a,i,s}/\sum_a\sum_s \hat{I}^N_{y,m,a,i,s}`$ | Composition | See next table | 7.6 |
| Index at length | $`p^{IN}_{y,m,\ell,i}`$ | $`\sum_s \hat{I}^N_{y,m,\ell,i,s}/\sum_\ell\sum_s \hat{I}^N_{y,m,\ell,i,s}`$ | Composition | See next table | 7.7 |
| Stock composition | $`p^{SC}_{y,m,a,f,r,s}`$ | $`\hat{C}^N_{y,m,a,f,r,s}/\sum_s \hat{C}^N_{y,m,a,f,r,s}`$ | Composition | See next table | 7.8 |
| Parent-offspring pairs | $`p^{POP}_{y,t,a,s}`$ | $`\hat{p}^{POP}_{y,t,a,s}`$ | Binomial | $`N\hat{p}^{POP}(1 - \hat{p}^{POP})`$ | 7.9 |
| Half-sibling pairs | $`p^{HSP}_{i,j,s}`$ | $`\hat{p}^{HSP}_{i,j,s}`$ | Binomial | $`N\hat{p}^{HSP}(1 - \hat{p}^{HSP})`$ | 7.10 |
| Tag | $`p^{mov}_{y,m,a,r,r',s}`$ | $`\hat{mov}_{y,m,a,r,r',s}`$ | Multinomial | See next table | 7.11 |

Potential distributions for the likelihoods of composition data, which
are presented as proportions, and the predicted mean and variance. $`N`$
is the sample size for each composition vector and $`\theta`$ is a
tuning parameter for the Dirichlet-multinomial distribution, both
provided as user inputs. For lognormal and logitnormal likelihoods,
$`\sigma`$ is a user-provided variance term. $`N`$ is unique to each
vector observation, e.g., age composition by season, fleet, and region
while $`\theta`$ is unique to fleet or survey.

| Distribution | Mean | Variance | Number |
|:---|:---|:---|:---|
| Multinomial | $`N\hat{p}`$ | $`N\hat{p}(1-\hat{p})`$ | 8.1 |
| Dirichlet-multinomial (Type 1) | $`N\hat{p}`$ | $`N\hat{p}(1-\hat{p})\dfrac{N+\theta N}{1+\theta N}`$ | 8.2 |
| Dirichlet-multinomial (Type 2) | $`N\hat{p}`$ | $`N\hat{p}(1-\hat{p})\dfrac{N+\theta}{1+\theta}`$ | 8.3 |
| Lognormal | $`\log(\hat{p})`$ | $`\hat{p}^{-1}`$ or $`\sigma`$ | 8.4 |
| Logitnormal | $`\log\left(\frac{\hat{p}}{\hat{1-p}}\right)`$ | $`\sigma`$ | 8.5 |

### Priors

Prior distributions for various parameters are described here.

| Description | Distribution | Equation | Number |
|:---|:---|:---|:---|
| Deviations from the equilibrium age structure | Lognormal | $`x^{\textrm{Req}}_{a,s} \sim N(0, \sigma^R_s)`$ | 9.1 |
| Recruitment deviates | Lognormal | $`x^R_{y,s} \sim N(0, \sigma^R_s)`$ | 9.2 |

### Penalty function

A quadratic penalty to the objective function when any $`F_{y,m,f,r}`$
exceeds the specified maximum.

``` math
\textrm{Penalty} = \sum_y\sum_m\sum_f\sum_r
\begin{cases}
0.1 (F_{max} - F_{y,m,f,r})^2 & F_{y,m,f,r} \ge F_{max}\\
0 & \textrm{otherwise}
\end{cases}
```
