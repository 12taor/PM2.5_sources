APCA
================
Rachel Tao
1/2/2021

Load in data

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   .default = col_double(),
    ##   date = col_character()
    ## )
    ## ℹ Use `spec()` for the full column specifications.

Means, correlations and variances

``` r
conc_corr <- 
  conc %>% 
  select(-date)

corr <- 
  summarize(conc_corr, corr = cor(conc_corr, use = "pairwise.complete.obs"))

conc_summarize <- 
  conc_corr %>% 
  pivot_longer(everything(),
               names_to = "element",
               values_to = "concentration") %>% 
  group_by(element)

means <- conc_summarize %>% 
  summarize(mean = mean(concentration, na.rm = TRUE)) %>% 
  arrange(element)
```

    ## `summarise()` ungrouping output (override with `.groups` argument)

``` r
vars <- conc_summarize %>% 
  summarize(variance = var(concentration, na.rm = TRUE)) %>% 
  arrange(element)
```

    ## `summarise()` ungrouping output (override with `.groups` argument)

Scaled factor analysis - must provide number of factors here

``` r
# set mean = 0, SD = 1
# remove potassium so that 4th of July fireworks do not alter overall results
zspec_pm25 <- 
  conc %>% 
  select(-date, -pm25, -K) %>% 
  scale(center = TRUE, scale = TRUE)

# factor analysis - must provide number of factors!
factors <- 
  fa(zspec_pm25,
     nfactors = 6, max.iter = 500,
     rotate = "varimax",
     scores = "regression",
     SMC  = TRUE,
     fm = "pa")
```

    ## Warning in fa.stats(r = r, f = f, phi = phi, n.obs = n.obs, np.obs = np.obs, :
    ## The estimated weights for the factor scores are probably incorrect. Try a
    ## different factor score estimation method.

    ## Warning in fac(r = r, nfactors = nfactors, n.obs = n.obs, rotate = rotate, : An
    ## ultra-Heywood case was detected. Examine the results carefully

``` r
# calculate factor scores, loadings, and weights (aka standardized scoring coefficients)
fact <- as_tibble(factors$scores)
loadings <- factors$loadings
factor_stat <- factors$weights %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("element") %>% 
  arrange(element)
```

calculate absolute factor scores for each day (Equations 5-7 in
statistical methods of Kavouras 2001)

``` r
# join standardized scoring coefficients with means and variances
coeff <- left_join(means, vars, by = "element") %>% 
  left_join(factor_stat, by = "element")

# calculate minimum values to find absolute factor scores (equations 5 and 6)
min <- coeff %>% 
  pivot_longer(contains("PA"),
               names_to = "factor_num",
               names_prefix = "PA",
               values_to = "factor_coeff") %>% 
  mutate(minimum = factor_coeff*(-mean/(variance^(1/2)))) %>% 
  select(element, factor_num, minimum) %>% 
  group_by(factor_num) %>% 
  summarize(sum = sum(minimum, na.rm = TRUE)) %>% 
  pivot_wider(factor_num:sum,
              names_from = factor_num,
              names_prefix = "minimum_",
              values_from = sum)
```

    ## `summarise()` ungrouping output (override with `.groups` argument)

``` r
# reattach factor scores to dates and calculate absolute factor scores for each day
# (equation 7)
a <- conc %>% 
  select(date) %>% 
  cbind(fact) %>% 
  pivot_longer(contains("PA"),
               names_to = "factor_num",
               names_prefix = "PA",
               values_to = "factor_score") %>% 
  mutate(factor_num = as.numeric(factor_num)) %>% 
  cbind(min) %>% 
  pivot_longer(contains("minimum"),
               names_to = "min_num",
               names_prefix = "minimum_",
               values_to = "minimum") %>% 
  mutate(min_num = as.numeric(min_num)) %>% 
  filter(factor_num == min_num) %>% 
  mutate(factor_score = factor_score - minimum) %>% 
  select(-min_num, -minimum)
```

factor means

``` r
# this is a check - not necessary for next steps
factor_means <- a %>% 
  select(factor_num, factor_score) %>% 
  group_by(factor_num) %>% 
  summarize(mean_score = mean(factor_score, na.rm = TRUE))
```

    ## `summarise()` ungrouping output (override with `.groups` argument)

Regress each element on factors to see source contribution

``` r
# define 'regress' as regression function with y as dependent variable and x as independent variable(s)
# return tidy output for beta values and r-squared value
regress = function(y) {
  
  z = lm(y ~ x)
  
  r.squared = broom::glance(z) %>% select(r.squared)

  betas = broom::tidy(z) %>% pivot_wider(term:estimate,
                                         names_from = "term",
                                         values_from = "estimate")
  
  return(cbind(r.squared, betas))
  
}

# independent variables for regression will be factors from factor analysis. define these as 'x'
x <- a %>% 
  pivot_wider(everything(),
              names_from = factor_num,
              names_prefix = "factor_",
              values_from = factor_score) %>% 
  select(contains("factor")) %>% as.matrix()

# dependent variable for regression will be actual element concentration. define this as y.
y <- conc %>% select(-date)

# repeat 'regress' for each element using 'x' as factors and 'y' as concentrations)
regress <- map(y, regress) %>% 
  enframe(name = "element") %>% 
  unnest(c(value)) %>% 
  janitor::clean_names() %>% 
  rename_at(vars(contains("factor")), funs(str_replace(., "xfactor", "beta")))
```

    ## Warning: `funs()` is deprecated as of dplyr 0.8.0.
    ## Please use a list of either functions or lambdas: 
    ## 
    ##   # Simple named list: 
    ##   list(mean = mean, median = median)
    ## 
    ##   # Auto named with `tibble::lst()`: 
    ##   tibble::lst(mean, median)
    ## 
    ##   # Using lambdas
    ##   list(~ mean(., trim = .2), ~ median(., na.rm = TRUE))
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_warnings()` to see where this warning was generated.

``` r
# dataset storing r-squred values
r <- regress %>% 
  select(element, r_squared)

# dataset storing beta values
betas <- regress %>% 
  select(element, starts_with("beta")) %>% 
  nest(betas = starts_with("beta"))
```

calculate source contributions

``` r
# combine betas with factor scores in a list 
sources = lapply(betas$betas, cbind, a)

# in the above step, element names are automatically removed. rename list elements
names(sources) <- betas$element

# calculate source contributions by multiplying absolute factor scores by beta values
sources <- lapply(sources, as_tibble) %>% 
  bind_rows(.id = "element") %>%
  pivot_longer(starts_with("beta"),
               names_to = "beta_num",
               names_prefix = "beta_",
               values_to = "beta_value") %>% 
  mutate(
    beta_num = as.numeric(beta_num)
  ) %>% 
  filter(factor_num == beta_num) %>% 
  mutate(
    srce = factor_score*beta_value
  ) %>% 
  mutate(srce_num = as.character(beta_num)) %>% 
  select(element, date, srce_num, srce) 

conc_long <- conc %>% 
  pivot_longer(!date,
               names_to = "element",
               values_to = "concentration") %>% 
  arrange(element)

source_conc <- left_join(sources, conc_long, by = c('element', 'date'))

# calculate mean source contributions for each element (is this equivalent to laodings?)
mean_source <- source_conc %>% 
  group_by(element, srce_num) %>% 
  summarize(source_mean = mean(srce, na.rm = TRUE),
            MeanConc = mean(concentration, na.rm = TRUE)) %>% 
  pivot_wider(names_from = "srce_num",
              names_prefix = "source_",
              values_from = "source_mean")
```

    ## `summarise()` regrouping output by 'element' (override with `.groups` argument)

### Final datasets:

1)  original loadings from factor analysis
2)  mean source contributions for each element, mean predicted
    concentration, percent error
3)  source contributions to total pm2.5

<!-- end list -->

``` r
# results of factor analysis, including loadings 
factors
```

    ## Factor Analysis using method =  pa
    ## Call: fa(r = zspec_pm25, nfactors = 6, rotate = "varimax", scores = "regression", 
    ##     SMC = TRUE, max.iter = 500, fm = "pa")
    ## Standardized loadings (pattern matrix) based upon correlation matrix
    ##      PA1   PA2   PA4   PA6   PA5   PA3    h2     u2 com
    ## bc  0.11  0.29  0.51  0.27  0.14 -0.02 0.450  0.550 2.6
    ## Na  0.21  0.12  0.46 -0.05 -0.02  0.76 0.848  0.152 1.9
    ## Al  0.93  0.01  0.19  0.06  0.21  0.04 0.945  0.055 1.2
    ## Si  0.94  0.02 -0.01  0.18  0.08 -0.01 0.923  0.077 1.1
    ## S   0.25  0.11  0.84 -0.02  0.05  0.07 0.792  0.208 1.2
    ## Cl -0.03 -0.03 -0.15  0.03  0.05  0.74 0.582  0.418 1.1
    ## Ca  0.61  0.20  0.08  0.40  0.12  0.17 0.624  0.376 2.3
    ## Ti  0.56  0.04  0.10  0.18  0.54  0.01 0.647  0.353 2.3
    ## V   0.08  0.88  0.19  0.11  0.10  0.03 0.838  0.162 1.2
    ## Cr  0.07  0.02  0.08  0.16  0.18 -0.03 0.071  0.929 2.8
    ## Mn  0.08  0.09  0.03  0.50  0.03  0.02 0.268  0.732 1.1
    ## Fe  0.49  0.19  0.20  0.83  0.19 -0.02 1.035 -0.035 2.0
    ## Ni  0.02  0.89  0.07  0.14  0.06  0.03 0.816  0.184 1.1
    ## Cu  0.04  0.08  0.11  0.09  0.64  0.03 0.435  0.565 1.2
    ## Zn  0.02  0.33  0.22  0.32  0.20 -0.02 0.294  0.706 3.5
    ## Se -0.02  0.00  0.15  0.02  0.04 -0.02 0.024  0.976 1.3
    ## Br  0.06  0.16  0.30  0.12  0.17  0.12 0.174  0.826 3.1
    ## Ba  0.09  0.02  0.05 -0.01  0.44  0.02 0.207  0.793 1.1
    ## Pb  0.11  0.21  0.23  0.13  0.31  0.01 0.221  0.779 3.4
    ## 
    ##                        PA1  PA2  PA4  PA6  PA5  PA3
    ## SS loadings           2.83 1.94 1.57 1.44 1.24 1.18
    ## Proportion Var        0.15 0.10 0.08 0.08 0.07 0.06
    ## Cumulative Var        0.15 0.25 0.33 0.41 0.47 0.54
    ## Proportion Explained  0.28 0.19 0.15 0.14 0.12 0.12
    ## Cumulative Proportion 0.28 0.47 0.62 0.76 0.88 1.00
    ## 
    ## Mean item complexity =  1.9
    ## Test of the hypothesis that 6 factors are sufficient.
    ## 
    ## The degrees of freedom for the null model are  171  and the objective function was  8.82 with Chi Square of  38597.72
    ## The degrees of freedom for the model are 72  and the objective function was  0.42 
    ## 
    ## The root mean square of the residuals (RMSR) is  0.02 
    ## The df corrected root mean square of the residuals is  0.03 
    ## 
    ## The harmonic number of observations is  3886 with the empirical chi square  585.45  with prob <  1.7e-81 
    ## The total number of observations was  4383  with Likelihood Chi Square =  1847.57  with prob <  0 
    ## 
    ## Tucker Lewis Index of factoring reliability =  0.89
    ## RMSEA index =  0.075  and the 90 % confidence intervals are  0.072 0.078
    ## BIC =  1243.81
    ## Fit based upon off diagonal values = 0.99

``` r
# calculate mean predicted concentration and percent error
SA_Varimax <- mean_source %>% 
  ungroup(element) %>% 
  left_join(r, by = "element") %>% 
  mutate(
    PredConc = reduce(select(., contains("source")), `+`),
    Pct_error = (PredConc - MeanConc)/MeanConc
  )

# source contributions to total pm2.5
pm25 <- source_conc %>% 
  filter(element == 'pm25') %>% 
  select(-element) %>% 
  pivot_wider(names_from = srce_num,
              names_prefix = "source_",
              values_from = srce)
```

### Next steps (likely to copy into new doc):

  - Do general time series of PM2.5 and MI
  - Total PM2.5 and MI
  - APCA, PCP, PCP + APCA?
  - when running APCA for NYC, extract loadings and for each source
    contribution look at seasonal patterns (break calendar into 4
    seasons of 3 months each), look at concentrations/contributions by
    season (traffic higher in winter, regional in winter), day of week
    (traffic higher on weekdays)
  - scores: correlations between sources and PM2.5 + patterns
