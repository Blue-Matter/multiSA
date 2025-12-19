# Model equations

## Variable definitions

### Equation subscripts

The following letters are used for subscripts to identify the dimension
and indexing for model variables.

| Subscript | Definition            | Number |
|:----------|:----------------------|:-------|
| $y$       | Year                  | 1.1    |
| $m$       | Season                | 1.2    |
| $a$       | Age                   | 1.3    |
| $r$       | Region (spatial area) | 1.4    |
| $f$       | Fleet                 | 1.5    |
| $s$       | Stock                 | 1.6    |
| $i$       | Index of abundance    | 1.7    |
| $\ell$    | Length bin            | 1.8    |

### Fixed parameters

The parameters here are set up as fixed inputs prior to fitting the
model.

| Description                                                  | Symbol                                        | Number |
|:-------------------------------------------------------------|:----------------------------------------------|:-------|
| Stock weight at age                                          | $w_{y,m,a,s}$                                 | 2.1    |
| Fishery weight at age                                        | $w_{y,m,a,f,s}^{F}$                           | 2.2    |
| Maturity ogive                                               | $m_{y,a,s}$                                   | 2.3    |
| Fecundity                                                    | $f_{y,a,s}$                                   | 2.4    |
| Natural mortality (instantaneous per year)                   | $M_{y,a,s}$                                   | 2.5    |
| Length-at-age probability                                    | $\text{Pr}\left( \ell \mid a \right)_{y,m,s}$ | 2.6    |
| Fractional parameter (between 0 - 1)                         | $\Delta$                                      | 2.7    |
| Season length (relative to year)                             | $\Delta_{m}$                                  | 2.8    |
| Spawn timing (within season)                                 | $\Delta_{sp}$                                 | 2.9    |
| Season of spawning (subscript)                               | $m_{\text{sp}}$                               | 2.10   |
| Season of recruitment (subscript)                            | $m_{\text{rec}}$                              | 2.11   |
| Natal spawning (proportion of mature individuals that spawn) | $n_{r,s}$                                     | 2.12   |
| Stock scaling parameter                                      | $r_{s}$                                       | 2.13   |
| Survey timing                                                | $\Delta_{i}$                                  | 2.14   |
| Survey sampling coverage                                     | $\delta_{i,r,s}$                              | 2.15   |

### Estimated parameters

The parameters here are set up to be either estimated or fixed in the
model. Such parameters are identified as $x$ which is estimated over all
real numbers and transformed to the appropriate model parameter
described below.

- Unfished recruitment is scaled by an additional user parameter $r_{s}$
  which is intended to aid convergence. For multi-stock models, $r_{s}$
  should be proportional to the expected stock size, i.e., large values
  for large stocks.
- Maturity and natural mortality can either be estimated or defined by
  the user as above.
- Fishing mortality is fixed to zero when the corresponding catch is
  less than $10^{- 8}$.

| Description                                                                           | Symbol                                                                                                         | Number |
|:--------------------------------------------------------------------------------------|:---------------------------------------------------------------------------------------------------------------|:-------|
| Estimated parameter (all real numbers)                                                | $x$                                                                                                            | 3.1    |
| Unfished recruitment                                                                  | $R_{0,s} = r_{s}\exp\left( x_{s}^{R0} \right)$                                                                 | 3.2    |
| Beverton-Holt stock-recruit steepness                                                 | $h_{s} = 0.8 \times \text{logit}^{- 1}\left( x_{s}^{h} \right) + 0.2$                                          | 3.3    |
| Ricker stock-recruit steepness                                                        | $h_{s} = \exp\left( x_{s}^{h} \right) + 0.2$                                                                   | 3.4    |
| Age of 50% maturity ($A$ is the maximum age)                                          | $a_{s}^{50} = A \times \text{logit}^{- 1}\left( x_{s}^{a50} \right)$                                           | 3.5    |
| Age of 95% maturity                                                                   | $a_{s}^{95} = a_{s}^{50} + \exp\left( x_{s}^{a95} \right)$                                                     | 3.6    |
| Natural mortality                                                                     | $M_{s} = \exp\left( x_{s}^{M} \right)$                                                                         | 3.7    |
| Recruitment deviates                                                                  | $x_{y,s}^{R}$                                                                                                  | 3.8    |
| Standard deviation of recruitment deviates                                            | $\sigma_{s}^{R} = \exp\left( x_{s}^{\sigma_{R}} \right)$                                                       | 3.9    |
| Fleet catchability by stock                                                           | $q_{f,s}^{F} = \exp\left( x_{f,s}^{qF} \right)$                                                                | 3.10   |
| Fishing mortality                                                                     | $F_{y,m,f,r} = \begin{cases}                                                                                   
                                                                                         {\exp\left( x_{f}^{\text{Fmult}} \right)\quad} & {y = y_{\text{ref}},m = m_{\text{ref}},r = r_{\text{ref}}} \\  
                                                                                         {\exp\left( x_{f}^{\text{Fmult}} + x_{y,m,r}^{\text{Fdev}} \right)\quad} & \text{otherwise}                     
                                                                                         \end{cases}$                                                                                                    | 3.11   |
| Selectivity - length of full selectivity ($f$ and $i$ subscripts are interchangeable) | $\mu_{f} = L_{\text{max}} \times \text{logit}^{- 1}\left( x_{f}^{\mu} \right)$                                 | 3.12   |
| Selectivity - width of ascending limb                                                 | $\sigma_{f}^{\text{asc}} = \exp\left( x_{f}^{\text{asc}} \right)$                                              | 3.13   |
| Selectivity - width of descending limb                                                | $\sigma_{f}^{\text{dsc}} = \exp\left( x_{f}^{\text{dsc}} \right)$                                              | 3.14   |
| Base movement parameters from origin $r$ to destination $r\prime$                     | $x_{m,a,r,r\prime,s}^{b}$                                                                                      | 3.15   |
| Attractivity term - movement                                                          | $x_{y,m,a,r\prime,s}^{g}$                                                                                      | 3.16   |
| Viscosity term - movement                                                             | $x_{y,m,a,s}^{v}$                                                                                              | 3.17   |
| Initial equilibrium fishing mortality                                                 | $F_{m,f,r}^{\text{eq}} = \exp\left( x_{m,f,r}^{\text{Feq}} \right)$                                            | 3.18   |
| Deviations from the equilibrium age structure                                         | $x_{a,s}^{\text{Req}}$                                                                                         | 3.19   |

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
  ($\phi_{0}$). In seasonal and multi-region models, the population
  dynamics model is used to numerically obtain $\phi_{0}$ by setting
  $F_{y,m,f,r} = 0$, recruitment to 1, and all other parameters to
  constant seasonal values. $\phi_{0}$ is the equilibrium spawning
  output at the end of this numerical spool-up.

| Description                            | Equation                                                                                                                                                                                                                   | Number |
|:---------------------------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:-------|
| Selectivity at length                  | $s_{\ell,f} = \begin{cases}                                                                                                                                                                                                
                                          {\exp\left( - 0.5\left\lbrack \frac{L_{\ell} - \mu_{f}}{\sigma_{f}^{\text{asc}}} \right\rbrack^{2} \right)} & {,L_{\ell} < \mu_{f}} \\                                                                                      
                                          {\exp\left( - 0.5\left\lbrack \frac{L_{\ell} - \mu_{f}}{\sigma_{f}^{\text{dsc}}} \right\rbrack^{2} \right)} & {,L_{\ell} \geq \mu_{f}}                                                                                      
                                          \end{cases}$                                                                                                                                                                                                                | 4.1    |
| Selectivity at age                     | $s_{y,m,a,f,s} = \sum_{\ell}\text{Pr}\left( \ell \mid a \right)_{y,m,s} \times s_{\ell,f}$                                                                                                                                 | 4.2    |
| Maturity (if estimated)                | $m_{a,s} = \left\lbrack 1 + \exp\left( - \log(19)\frac{a - a_{s}^{50}}{a_{s}^{95} - a_{s}^{50}} \right) \right\rbrack^{- 1}$                                                                                               | 4.3    |
| Movement from $r$ to $r\prime$         | $\text{mov}_{y,m,a,r,r\prime,s} = \frac{\exp\left( x_{m,a,r,r\prime,s}^{b} + x_{y,m,a,r,s}^{g} + x_{y,m,a,s}^{v} \right)}{\sum_{r\prime}\exp\left( x_{m,a,r,r\prime,s}^{b} + x_{y,m,a,r,s}^{g} + x_{y,m,a,s}^{v} \right)}$ | 4.4    |
| Beverton-Holt stock recruit parameters | $\begin{aligned}                                                                                                                                                                                                           
                                          \alpha_{s} & {= \frac{4h_{s}}{\left( 1 - h_{s} \right)\phi_{0,s}}} \\                                                                                                                                                       
                                          \beta_{s} & {= \frac{5h_{s} - 1}{\left( 1 - h_{s} \right)R_{0,s}\phi_{0,s}}}                                                                                                                                                
                                          \end{aligned}$                                                                                                                                                                                                              | 4.5    |
| Ricker stock recruit parameters        | $\begin{aligned}                                                                                                                                                                                                           
                                          \alpha_{s} & {= \frac{5h^{1.25}}{\phi_{0,s}}} \\                                                                                                                                                                            
                                          \beta_{s} & {= \frac{\log\left( \left\lbrack 5h_{s} \right\rbrack^{1.25} \right)}{R_{0,s}\phi_{0,s}}}                                                                                                                       
                                          \end{aligned}$                                                                                                                                                                                                              | 4.6    |

## Population dynamics

The following equations project the population forward in time.

- To obtain the initial abundance $N_{y = 1,m = 1,a,r,s}$ array in
  seasonal and multi-region models, a numerical spool-up is performed
  with the seasonal fishing mortality equal to $F_{m,f,r}^{\text{eq}}$,
  recruitment to 1, and all other parameters set to constant seasonal
  values from the first year of the model. From this initialization, the
  equilibrium spawners per recruit $\phi_{eq}$ is the final spawning
  output, and the seasonal numbers per recruit
  $\text{NPR}_{m,a,r,s}^{\text{eq}}$ is obtained from the abundance
  array. The initial abundance is the product of the equilibrium
  recruitment and numbers per recruit.
- It is possible that some proportion of the mature population do not
  contribute to the annual spawning based on the natal spawning
  parameter specifying the spatial spawning pattern. Thus, there is a
  distinction between potential spawners and realized spawners. The
  unfished replacement line of the stock-recruit relationship
  ($R_{s} = 1/\phi_{0}$) is based on the realized spawning in
  equilibrium.

| Description                                                                                 | Equation                                                                                                                                        | Number |
|:--------------------------------------------------------------------------------------------|:------------------------------------------------------------------------------------------------------------------------------------------------|:-------|
| Index of fishing effort (instantaneous per season)                                          | $F_{y,m,f,r}$                                                                                                                                   | 5.1    |
| Stock abundance                                                                             | $N_{y,m,a,r,s}$                                                                                                                                 | 5.2    |
| Equilibrium recruitment                                                                     | $R_{\text{eq},s} = \begin{cases}                                                                                                                
                                                                                               \frac{\alpha_{s}\phi_{eq,s} - 1}{\beta_{s}\phi_{eq,s}} & \text{Beverton-Holt} \\                                                                 
                                                                                               \frac{\log\left( \alpha_{s}\phi_{eq,s} \right)}{\beta_{s}\phi_{eq,s}} & \text{Ricker}                                                            
                                                                                               \end{cases}$                                                                                                                                     | 5.3    |
| Initial abundance                                                                           | $N_{y = 1,m = 1,a,r,s} = R_{\text{eq},s}\exp\left( x_{a,s}^{\text{Req}} \right) \times \text{NPR}_{m = 1,a,r,s}^{\text{eq}}$                    | 5.4    |
| Fishing mortality by time, age, fleet, region, stock                                        | $F_{y,m,a,f,r,s} = q_{f,s}^{F}s_{y,m,a,f,s}F_{y,m,f,r}$                                                                                         | 5.5    |
| Total mortality (instantaneous per season)                                                  | $Z_{y,m,a,r,s} = \Delta_{m}M_{y,a,s} + \sum_{f}F_{y,m,a,f,r,s}$                                                                                 | 5.6    |
| Seasonal abundance without incoming recruitment (after survival and movement)               | $N_{y,m,a,r\prime,s} = \sum_{r}N_{y,m - 1,a,r,s}\exp\left( - Z_{y,m - 1,a,r,s} \right)\text{mov}_{y,m,a,r,r\prime,s}$                           | 5.7    |
| Potential spawners (PS)                                                                     | $N_{y,a,r,s}^{\text{PS}} = N_{y,m = m_{\text{sp}},a,r,s}\exp\left( - \Delta_{\text{sp}}Z_{y,m = m_{\text{sp}},a,r,s} \right)m_{y,a,s}$          | 5.8    |
| Active spawners                                                                             | $N_{y,a,r,s}^{\text{S}} = N_{y,m = m_{\text{sp}},a,r,s}^{\text{PS}} \times n_{r,s}$                                                             | 5.9    |
| Spawning output                                                                             | $S_{y,r,s} = \sum_{a}N_{y,a,r,s}^{\text{S}}f_{y,a,s}$                                                                                           | 5.10   |
| Recruitment: Beverton-Holt                                                                  | $R_{y,s} = \frac{\alpha_{s}\sum_{r}S_{y,r,s}}{1 + \beta_{s}\sum_{r}S_{y,r,s}}\exp\left( x_{y,s}^{R} \right)$                                    | 5.11   |
| Recruitment: Ricker                                                                         | $R_{y,s} = \alpha_{s}\sum_{r}S_{y,r,s}\exp\left( - \beta_{s} \times \sum_{r}S_{y,r,s} \right)\exp\left( x_{y,s}^{R} \right)$                    | 5.12   |
| Seasonal abundance (incoming recruitment and advancing age class when $m = m_{\text{rec}}$) | $N_{y,m,a,r\prime,s} = \begin{cases}                                                                                                            
                                                                                               {R_{y,s}\text{mov}_{y,m,1,r,r\prime,s}} & {a = 1} \\                                                                                             
                                                                                               {\sum_{r}N_{y,m - 1,a - 1,r,s}\exp\left( - Z_{y,m - 1,a - 1,r,s} \right)\text{mov}_{y,m,a,r,r\prime,s}} & {a = 2,\ldots,A - 1} \\                
                                                                                               {\sum_{a\prime = A - 1}^{A}\sum_{r}N_{y,m - 1,a\prime,r,s}\exp\left( - Z_{y,m - 1,a\prime,r,s} \right)\text{mov}_{y,m,a,r,r\prime,s}} & {a = A}  
                                                                                               \end{cases}$                                                                                                                                     | 5.13   |

## Report variables

Here, we calculate additional variables that are not needed for the
population dynamics model, but are of interest for fitting the model or
for reporting.

- In a multi-region and/or seasonal model, we may want a summary fishing
  mortality (per year) across all regions and fleets ($F_{y,a,s}$) which
  calculated from the Baranov equation with natural mortality
  $M_{y,a,s}$, total stock abundance at the beginning of the year
  $N_{y,a,s}$, and total catch $C_{y,a,s}^{N}$. The summary total
  mortality (per year) is then $Z_{y,a,s} = F_{y,a,s} + M_{y,a,s}$.
- Vulnerable biomass is the availability of the stock to individual
  fleets.
- When fitting to close-kin genetic data, we can calculate the
  probability of parent-offspring pairs (POP) with the cohort year of
  the offspring is $y$, the parental age at capture is $a\prime$, and
  the capture year of the parent $t$.
- The half-offspring pair probability is calculated from the parental
  probability in years $i$ and $j$, the cohort year of the older and
  younger sibling, respectively, and the parental survival from year $i$
  to year $j$. The parental age is not observed, so we calculate the
  probability across all potential ages and follow each cohort from $i$
  to $j$.

| Description                                            | Equation                                                                                                                                                                  | Number |
|:-------------------------------------------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:-------|
| Equilibrium catch (abundance, age)                     | $C_{m,a,f,r,s}^{Neq} = \frac{F_{m,a,f,r,s}^{\text{eq}}}{Z_{m,a,r,s}^{\text{eq}}}\left( 1 - \exp\left( - Z_{y,m,a,r,s}^{\text{eq}} \right) \right)N_{m,a,r,s}^{\text{eq}}$ | 6.1    |
| Equilibrium catch (biomass)                            | $C_{m,f,r,s}^{Beq} = \sum_{a}w_{y = 1,m,a,f,s}^{F}C_{m,a,f,r,s}^{Neq}$                                                                                                    | 6.2    |
| Catch (abundance, age)                                 | $C_{y,m,a,f,r,s}^{N} = \frac{F_{y,m,a,f,r,s}}{Z_{y,m,a,r,s}}\left( 1 - \exp\left( - Z_{y,m,a,r,s} \right) \right)N_{y,m,a,r,s}$                                           | 6.3    |
| Catch (abundance, length)                              | $C_{y,m,l,f,r,s}^{N} = \sum_{a}C_{y,m,a,f,r,s}^{N}\text{Pr}\left( \ell \mid a \right)_{y,m,s}$                                                                            | 6.4    |
| Catch (biomass)                                        | $C_{y,m,f,r,s}^{B} = \sum_{a}w_{y,m,a,f,s}^{F}C_{y,m,a,f,r,s}^{N}$                                                                                                        | 6.5    |
| Total biomass                                          | $B_{y,m,r,s} = \sum_{a}w_{y,m,a,s}N_{y,m,a,r,s}$                                                                                                                          | 6.6    |
| Vulnerable biomass                                     | $V_{y,m,r,s} = \sum_{a}s_{y,m,a,f,s}w_{y,m,a,f,s}^{F}N_{y,m,a,r,s}$                                                                                                       | 6.7    |
| Index age composition                                  | $I_{y,m,a,i,s}^{A} = q_{i}\sum_{r}s_{y,m,a,i,s}N_{y,m,a,r,s}\exp\left( - \Delta_{i} \times Z_{y,m,a,r,s} \right)\delta_{i,r,s}$                                           | 6.8    |
| Index length composition                               | $I_{y,m,\ell,i,s}^{L} = \sum_{a}I_{y,m,a,i,s}^{A}\text{Pr}\left( \ell \mid a \right)_{y,m,s}$                                                                             | 6.9    |
| Biomass index                                          | $I_{y,m,i} = q_{i}\sum_{a}\sum_{s}w_{y,m,a,s}I_{y,m,a,i,s}^{A}$                                                                                                           | 6.10   |
| Annual catch at age (all seasons, fleets, and regions) | $C_{y,a,s}^{N} = \sum_{m}\sum_{f}\sum_{r}C_{y,m,a,f,r,s}^{N}$                                                                                                             | 6.11   |
| Abundance at year start across regions                 | $N_{y,a,s} = \sum_{r}N_{y,m = 1,a,r,s}$                                                                                                                                   | 6.12   |
| Summary fishing mortality per year $F_{y,a,s}$         | $C_{y,a,s}^{N} = \frac{F_{y,a,s}}{F_{y,a,s} + M_{y,a,s}}\left( 1 - \exp\left( - \left\lbrack F_{y,a,s} + M_{y,a,s} \right\rbrack \right) \right)N_{y,a,s}$                | 6.13   |
| Parent-offspring probability                           | $\text{Pr}\left( \text{POP} \mid y,t,a\prime,s \right) = 2 \times \frac{f_{y,a = y - {(t - a\prime)},s}}{\sum_{a}f_{y,a,s}N_{y,a,s}^{\text{S}}}$                          | 6.14   |
| Half-sibling probability                               | $\text{Pr}\left( \text{HSP} \mid i,j,s \right) = 4 \times \sum_{a}c_{i,a,s} \times \text{surv}_{ijs} \times c_{j,a,s}$                                                    | 6.15   |
| Parental probability for older sibling $i$             | $c_{i,a,s} = \frac{N_{y = i,a,s}f_{y = i,a,s}}{\sum_{a\prime}N_{y = i,a\prime,s}f_{y = i,a\prime,s}}$                                                                     | 6.16   |
| Parental survival from year $i$ to $j$                 | $\text{surv}_{ijs} = \exp\left( - \sum\limits_{t = 0}^{j - i - 1}Z_{i + t,a + t,s} \right)$                                                                               | 6.17   |
| Parental probability for younger sibling $j$           | $c_{j,a,s} = \frac{f_{y = j,a,s}}{\sum_{a\prime}N_{y = j,a\prime,s}f_{y = j,a\prime,s}}$                                                                                  | 6.18   |

## Objective function

The objective function is the sum of the negative log-likelihoods,
negative log-priors, and penalty function.

### Likelihoods

The statistical distributions used for the likelihoods of the data are
described. The dimensions of the data are given below as well as the
corresponding model variable for the predicted value, which is typically
summed across stocks, (except for stock composition). Composition data
are presented as proportions $p$ and a separate table provides the mean
and variance of the various likelihood options.

The close-kin likelihood uses the ratio of matches (for either
parent-offspring or sibling matches) and the number of pairwise
comparisons ($N$).

| Data                        | Symbol                 | Predicted                                                                                          | Distribution | Variance                                                     | Number |
|:----------------------------|:-----------------------|:---------------------------------------------------------------------------------------------------|:-------------|:-------------------------------------------------------------|:-------|
| Equilibrium catch (biomass) | $C_{m,f,r}^{Beq}$      | $\sum_{s}{\widehat{C}}_{m,f,r,s}^{Beq}$                                                            | Lognormal    | $\sigma = 0.01$                                              | 7.1    |
| Catch (biomass)             | $C_{y,m,f,r}^{B}$      | $\sum_{s}{\widehat{C}}_{y,m,f,r,s}^{B}$                                                            | Lognormal    | $\sigma_{y,m,f,r}^{C}$                                       | 7.2    |
| Catch at age                | $p_{y,m,a,f,r}^{CN}$   | $\sum_{s}{\widehat{C}}_{y,m,a,f,r,s}^{N}/\sum_{a}\sum_{s}{\widehat{C}}_{y,m,a,f,r,s}^{N}$          | Composition  | See next table                                               | 7.3    |
| Catch at length             | $p_{y,m,l,f,r}^{CN}$   | $\sum_{s}{\widehat{C}}_{y,m,\ell,f,r,s}^{N}/\sum_{\ell}\sum_{s}{\widehat{C}}_{y,m,\ell,f,r,s}^{N}$ | Composition  | See next table                                               | 7.4    |
| Total indices               | $I_{y,m,i}$            | ${\widehat{I}}_{y,m,i}$                                                                            | Lognormal    | $\sigma_{y,m,i}^{I}$                                         | 7.5    |
| Index at age                | $p_{y,m,a,i}^{IN}$     | $\sum_{s}{\widehat{I}}_{y,m,a,i,s}^{N}/\sum_{a}\sum_{s}{\widehat{I}}_{y,m,a,i,s}^{N}$              | Composition  | See next table                                               | 7.6    |
| Index at length             | $p_{y,m,\ell,i}^{IN}$  | $\sum_{s}{\widehat{I}}_{y,m,\ell,i,s}^{N}/\sum_{\ell}\sum_{s}{\widehat{I}}_{y,m,\ell,i,s}^{N}$     | Composition  | See next table                                               | 7.7    |
| Stock composition           | $p_{y,m,a,f,r,s}^{SC}$ | ${\widehat{C}}_{y,m,a,f,r,s}^{N}/\sum_{s}{\widehat{C}}_{y,m,a,f,r,s}^{N}$                          | Composition  | See next table                                               | 7.8    |
| Parent-offspring pairs      | $p_{y,t,a,s}^{POP}$    | ${\widehat{p}}_{y,t,a,s}^{POP}$                                                                    | Binomial     | $N{\widehat{p}}^{POP}\left( 1 - {\widehat{p}}^{POP} \right)$ | 7.9    |
| Half-sibling pairs          | $p_{i,j,s}^{HSP}$      | ${\widehat{p}}_{i,j,s}^{HSP}$                                                                      | Binomial     | $N{\widehat{p}}^{HSP}\left( 1 - {\widehat{p}}^{HSP} \right)$ | 7.10   |
| Tag                         | TBD                    |                                                                                                    |              |                                                              | 7.11   |

Potential distributions for the likelihoods of composition data, which
are presented as proportions, and the predicted mean and variance. $N$
is the sample size for each composition vector and $\theta$ is a tuning
parameter for the Dirichlet-multinomial distribution, both provided as
user inputs. $N$ is unique to each vector observation, e.g., age
composition by season, fleet, and region while $\theta$ is unique to
fleet or survey.

| Distribution                            | Mean                             | Variance                                                                      | Number |
|:----------------------------------------|:---------------------------------|:------------------------------------------------------------------------------|:-------|
| Multinomial                             | $N\widehat{p}$                   | $N\widehat{p}\left( 1 - \widehat{p} \right)$                                  | 8.1    |
| Dirichlet-multinomial (Type 1)          | $N\widehat{p}$                   | $N\widehat{p}\left( 1 - \widehat{p} \right)\frac{N + \theta N}{1 + \theta N}$ | 8.2    |
| Dirichlet-multinomial (Type 2)          | $N\widehat{p}$                   | $N\widehat{p}\left( 1 - \widehat{p} \right)\frac{N + \theta}{1 + \theta}$     | 8.3    |
| Lognormal (summed across positive bins) | $\log\left( \widehat{p} \right)$ | ${\widehat{p}}^{- 1}$                                                         | 8.4    |

### Priors

Prior distributions for various parameters are described here.

| Description                                   | Distribution | Equation                                                                                          | Number |
|:----------------------------------------------|:-------------|:--------------------------------------------------------------------------------------------------|:-------|
| Deviations from the equilibrium age structure | Lognormal    | $x_{a,s}^{\text{Req}} \sim N\left( - 0.5\left( \sigma_{s}^{R} \right)^{2},\sigma_{s}^{R} \right)$ | 9.1    |
| Recruitment deviates                          | Lognormal    | $x_{y,s}^{R} \sim N\left( - 0.5\left( \sigma_{s}^{R} \right)^{2},\sigma_{s}^{R} \right)$          | 9.2    |

### Penalty function

A quadratic penalty to the objective function when any $F_{y,m,f,r}$
exceeds the specified maximum.

$$\text{Penalty} = \sum\limits_{y}\sum\limits_{m}\sum\limits_{f}\sum\limits_{r}\begin{cases}
{0.1\left( F_{max} - F_{y,m,f,r} \right)^{2}} & {F_{y,m,f,r} \geq F_{max}} \\
0 & \text{otherwise}
\end{cases}$$
