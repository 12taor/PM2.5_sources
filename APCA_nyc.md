APCA NYC
================
Rachel Tao
1/23/2021

Load in data

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   .default = col_double(),
    ##   date_local = col_date(format = "")
    ## )
    ## ℹ Use `spec()` for the full column specifications.

    ## Warning: Missing column names filled in: 'X1' [1]

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   .default = col_double(),
    ##   date_local = col_date(format = "")
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
  select(-date, -pm25) %>% 
  scale(center = TRUE, scale = TRUE)

# factor analysis - must provide number of factors!
factors <- 
  fa(zspec_pm25,
     nfactors = 6, max.iter = 500,
     rotate = "varimax",
     scores = "regression",
     SMC  = TRUE,
     fm = "pa")

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
    ##             PA1   PA2   PA3   PA4   PA5   PA6     h2    u2 com
    ## aluminum   0.05  0.68  0.05 -0.13  0.13  0.23 0.5523 0.448 1.4
    ## arsenic    0.39  0.08  0.21  0.29  0.05  0.09 0.3009 0.699 2.8
    ## barium     0.20  0.17  0.09  0.47 -0.08  0.01 0.3099 0.690 1.8
    ## bromine    0.61  0.09  0.24  0.14  0.02  0.30 0.5523 0.448 2.0
    ## cadmium    0.01  0.02 -0.02 -0.08 -0.02  0.00 0.0072 0.993 1.5
    ## calcium    0.60  0.50  0.06 -0.19 -0.04 -0.20 0.6960 0.304 2.4
    ## chromium  -0.02 -0.01 -0.04  0.01  0.43 -0.03 0.1859 0.814 1.0
    ## copper     0.50  0.24  0.14 -0.06  0.39  0.14 0.5096 0.490 2.8
    ## iron       0.47  0.53  0.18  0.10  0.58  0.01 0.8838 0.116 3.2
    ## lead       0.65  0.17  0.16  0.32 -0.05  0.04 0.5870 0.413 1.8
    ## magnesium  0.00  0.09  0.04 -0.01 -0.01  0.20 0.0483 0.952 1.5
    ## manganese  0.66  0.38  0.14 -0.06  0.19 -0.20 0.6817 0.318 2.1
    ## nickel     0.63  0.08  0.03  0.48 -0.05 -0.29 0.7254 0.275 2.4
    ## selenium   0.50  0.14  0.41  0.41 -0.05  0.04 0.6095 0.390 3.1
    ## silicon    0.13  0.83  0.19  0.09 -0.04  0.20 0.7938 0.206 1.3
    ## sulfate    0.24  0.20  0.91  0.20 -0.02  0.09 0.9699 0.030 1.4
    ## sulfur     0.24  0.20  0.91  0.20 -0.01  0.10 0.9790 0.021 1.4
    ## titanium   0.22  0.68  0.20  0.29  0.03  0.00 0.6299 0.370 1.8
    ## nitrates   0.67  0.01  0.23  0.13  0.07  0.01 0.5235 0.477 1.3
    ## vanadium   0.61  0.26  0.34  0.31  0.02 -0.16 0.6762 0.324 2.7
    ## zinc       0.84  0.05 -0.04 -0.09 -0.01  0.05 0.7166 0.283 1.0
    ## 
    ##                        PA1  PA2  PA3  PA4  PA5  PA6
    ## SS loadings           4.62 2.61 2.30 1.20 0.75 0.46
    ## Proportion Var        0.22 0.12 0.11 0.06 0.04 0.02
    ## Cumulative Var        0.22 0.34 0.45 0.51 0.55 0.57
    ## Proportion Explained  0.39 0.22 0.19 0.10 0.06 0.04
    ## Cumulative Proportion 0.39 0.61 0.80 0.90 0.96 1.00
    ## 
    ## Mean item complexity =  1.9
    ## Test of the hypothesis that 6 factors are sufficient.
    ## 
    ## The degrees of freedom for the null model are  210  and the objective function was  13.21 with Chi Square of  24049.42
    ## The degrees of freedom for the model are 99  and the objective function was  0.91 
    ## 
    ## The root mean square of the residuals (RMSR) is  0.02 
    ## The df corrected root mean square of the residuals is  0.04 
    ## 
    ## The harmonic number of observations is  1829 with the empirical chi square  476.87  with prob <  8.2e-51 
    ## The total number of observations was  1829  with Likelihood Chi Square =  1650.38  with prob <  1.5e-279 
    ## 
    ## Tucker Lewis Index of factoring reliability =  0.862
    ## RMSEA index =  0.093  and the 90 % confidence intervals are  0.089 0.097
    ## BIC =  906.74
    ## Fit based upon off diagonal values = 0.99
    ## Measures of factor score adequacy             
    ##                                                    PA1  PA2  PA3  PA4  PA5  PA6
    ## Correlation of (regression) scores with factors   0.94 0.92 0.98 0.82 0.86 0.72
    ## Multiple R square of scores with factors          0.89 0.85 0.95 0.67 0.75 0.52
    ## Minimum correlation of possible factor scores     0.77 0.70 0.91 0.34 0.49 0.03

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

### Next steps:

  - APCA, PCP, PCP + APCA?
  - when running APCA for NYC, extract loadings and for each source
    contribution look at seasonal patterns (break calendar into 4
    seasons of 3 months each), look at concentrations/contributions by
    season (traffic higher in winter, regional in winter), day of week
    (traffic higher on weekdays)
  - scores: correlations between sources and PM2.5 + patterns
