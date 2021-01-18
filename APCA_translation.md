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

Scaled factor analysis

``` r
# set mean =0, SD = 1

zspec_pm25 <- 
  conc %>% 
  select(-date, -pm25, -K) %>% 
  scale(center = TRUE, scale = TRUE)

# factor analysis

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
loadings <- print(factors$loadings, cutoff = 0)
```

    ## 
    ## Loadings:
    ##    PA1    PA2    PA4    PA6    PA5    PA3   
    ## bc  0.114  0.295  0.505  0.275  0.139 -0.020
    ## Na  0.207  0.123  0.463 -0.049 -0.017  0.757
    ## Al  0.927  0.009  0.190  0.056  0.213  0.039
    ## Si  0.940  0.022 -0.007  0.182  0.079 -0.010
    ## S   0.249  0.106  0.844 -0.019  0.046  0.070
    ## Cl -0.031 -0.030 -0.152  0.031  0.046  0.745
    ## Ca  0.610  0.203  0.081  0.397  0.125  0.173
    ## Ti  0.556  0.040  0.097  0.179  0.543  0.007
    ## V   0.082  0.880  0.185  0.112  0.097  0.028
    ## Cr  0.074  0.024  0.076  0.160  0.180 -0.031
    ## Mn  0.079  0.087  0.029  0.502  0.034  0.023
    ## Fe  0.488  0.188  0.202  0.827  0.189 -0.025
    ## Ni  0.023  0.887  0.074  0.138  0.055  0.025
    ## Cu  0.044  0.084  0.107  0.093  0.637  0.029
    ## Zn  0.023  0.325  0.219  0.316  0.198 -0.023
    ## Se -0.021  0.004  0.146  0.021  0.040 -0.019
    ## Br  0.059  0.157  0.300  0.117  0.170  0.116
    ## Ba  0.089  0.019  0.054 -0.006  0.443  0.023
    ## Pb  0.112  0.207  0.227  0.128  0.313  0.010
    ## 
    ##                  PA1   PA2   PA4   PA6   PA5   PA3
    ## SS loadings    2.826 1.944 1.568 1.438 1.236 1.184
    ## Proportion Var 0.149 0.102 0.083 0.076 0.065 0.062
    ## Cumulative Var 0.149 0.251 0.334 0.409 0.474 0.537

``` r
factor_stat <- factors$weights %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("element") %>% 
  arrange(element)
```

fix intercept (Kavouras 2001)

``` r
coeff <- left_join(means, vars, by = "element") %>% 
  left_join(factor_stat, by = "element")

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

reattach factor scores to dates

``` r
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

factor_means <- a %>% 
  select(factor_num, factor_score) %>% 
  group_by(factor_num) %>% 
  summarize(mean_score = mean(factor_score, na.rm = TRUE))
```

    ## `summarise()` ungrouping output (override with `.groups` argument)

Regress each element on factors to see source contribution

``` r
x <- a %>% 
  pivot_wider(everything(),
              names_from = factor_num,
              names_prefix = "factor_",
              values_from = factor_score) %>% 
  select(contains("factor")) %>% as.matrix()

regress = function(y) {
  
  z = lm(y ~ x)
  
  r.squared = broom::glance(z) %>% select(r.squared)

  betas = broom::tidy(z) %>% pivot_wider(term:estimate,
                                         names_from = "term",
                                         values_from = "estimate")
  
  return(cbind(r.squared, betas))
  
}

y <- conc %>% select(-date)


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
r <- regress %>% 
  select(element, r_squared)

betas <- regress %>% 
  select(element, starts_with("beta")) %>% 
  nest(betas = starts_with("beta"))
```

combine betas with factor scores for final datasets

``` r
sources = lapply(betas$betas, cbind, a)

names(sources) <- betas$element

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

mean_source <- source_conc %>% 
  group_by(element, srce_num) %>% 
  summarize(source_mean = mean(srce, na.rm = TRUE),
            MeanConc = mean(concentration, na.rm = TRUE)) %>% 
  pivot_wider(names_from = "srce_num",
              names_prefix = "source_",
              values_from = "source_mean")
```

    ## `summarise()` regrouping output by 'element' (override with `.groups` argument)

``` r
SA_Varimax <- mean_source %>% 
  ungroup(element) %>% 
  left_join(r, by = "element") %>% 
  mutate(
    PredConc = reduce(select(., contains("source")), `+`),
    Pct_error = (PredConc - MeanConc)/MeanConc
  )

pm25 <- source_conc %>% 
  filter(element == 'pm25') %>% 
  select(-element) %>% 
  pivot_wider(names_from = srce_num,
              names_prefix = "source_",
              values_from = srce)
```

also extract loadings

PM2.5 regression gives contribution of each source to total pm2.5

Do general time series of PM2.5 and MI Total PM2.5 and MI APCA, PCP, PCP
+ APCA? when running APCA for NYC, extract loadings and for each source
contribution look at seasonal patterns (break calendar into 4 seasons of
3 months each), look at concentrations/contributions by season (traffic
higher in winter, regional in winter), day of week (traffic higher on
weekdays) scores: correlations between sources and PM2.5 + patterns
