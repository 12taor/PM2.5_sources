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
zspec_pm25 <- 
  conc %>% 
  select(-date, -pm25, -K) %>% 
  scale(center = TRUE, scale = TRUE)

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

Regress each element (starting with aluminum) on factors to see source
contribution

What I want is to go from pm25 to Pb (8:28)

``` r
#this part is for iterating over each chemical

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

# regress(conc$Al)

# i <- c(8:28)

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
betas <- regress %>% 
  select(element, starts_with("beta")) %>% 
  nest(betas = starts_with("beta"))

# what happens to the r.squared value?
```

combine betas with factor scores for final datasets

``` r
# make a separate dataset for each chemical (use map?)


Al <- betas %>% filter(element == 'Al')


# 
# i <- c(1:21)
# betas[ i, ] <- apply(betas[ i, ], 1,
#                     function(x)
#                       as.character(x))
# 
# betas <- map(a[ , i], regress)


sources = lapply(betas$betas, cbind, a)

names(sources) <- betas$element

lapply(sources, as_tibble) 
```

    ## $pm25
    ## # A tibble: 26,298 x 9
    ##    beta_1 beta_2 beta_4 beta_6 beta_5 beta_3 date       factor_num factor_score
    ##     <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <date>          <dbl>        <dbl>
    ##  1  1542.  1343.  5728.   638.   953.   9.10 1999-01-01          1       NA    
    ##  2  1542.  1343.  5728.   638.   953.   9.10 1999-01-01          2       NA    
    ##  3  1542.  1343.  5728.   638.   953.   9.10 1999-01-01          4       NA    
    ##  4  1542.  1343.  5728.   638.   953.   9.10 1999-01-01          6       NA    
    ##  5  1542.  1343.  5728.   638.   953.   9.10 1999-01-01          5       NA    
    ##  6  1542.  1343.  5728.   638.   953.   9.10 1999-01-01          3       NA    
    ##  7  1542.  1343.  5728.   638.   953.   9.10 1999-01-02          1        0.413
    ##  8  1542.  1343.  5728.   638.   953.   9.10 1999-01-02          2        2.07 
    ##  9  1542.  1343.  5728.   638.   953.   9.10 1999-01-02          4        0.360
    ## 10  1542.  1343.  5728.   638.   953.   9.10 1999-01-02          6       -0.130
    ## # … with 26,288 more rows
    ## 
    ## $bc
    ## # A tibble: 26,298 x 9
    ##    beta_1 beta_2 beta_4 beta_6 beta_5 beta_3 date       factor_num factor_score
    ##     <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <date>          <dbl>        <dbl>
    ##  1   40.8   127.   249.   98.9   67.7  -25.9 1999-01-01          1       NA    
    ##  2   40.8   127.   249.   98.9   67.7  -25.9 1999-01-01          2       NA    
    ##  3   40.8   127.   249.   98.9   67.7  -25.9 1999-01-01          4       NA    
    ##  4   40.8   127.   249.   98.9   67.7  -25.9 1999-01-01          6       NA    
    ##  5   40.8   127.   249.   98.9   67.7  -25.9 1999-01-01          5       NA    
    ##  6   40.8   127.   249.   98.9   67.7  -25.9 1999-01-01          3       NA    
    ##  7   40.8   127.   249.   98.9   67.7  -25.9 1999-01-02          1        0.413
    ##  8   40.8   127.   249.   98.9   67.7  -25.9 1999-01-02          2        2.07 
    ##  9   40.8   127.   249.   98.9   67.7  -25.9 1999-01-02          4        0.360
    ## 10   40.8   127.   249.   98.9   67.7  -25.9 1999-01-02          6       -0.130
    ## # … with 26,288 more rows
    ## 
    ## $Na
    ## # A tibble: 26,298 x 9
    ##    beta_1 beta_2 beta_4 beta_6 beta_5 beta_3 date       factor_num factor_score
    ##     <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <date>          <dbl>        <dbl>
    ##  1   27.9   16.3   68.2  -5.02  -6.91   121. 1999-01-01          1       NA    
    ##  2   27.9   16.3   68.2  -5.02  -6.91   121. 1999-01-01          2       NA    
    ##  3   27.9   16.3   68.2  -5.02  -6.91   121. 1999-01-01          4       NA    
    ##  4   27.9   16.3   68.2  -5.02  -6.91   121. 1999-01-01          6       NA    
    ##  5   27.9   16.3   68.2  -5.02  -6.91   121. 1999-01-01          5       NA    
    ##  6   27.9   16.3   68.2  -5.02  -6.91   121. 1999-01-01          3       NA    
    ##  7   27.9   16.3   68.2  -5.02  -6.91   121. 1999-01-02          1        0.413
    ##  8   27.9   16.3   68.2  -5.02  -6.91   121. 1999-01-02          2        2.07 
    ##  9   27.9   16.3   68.2  -5.02  -6.91   121. 1999-01-02          4        0.360
    ## 10   27.9   16.3   68.2  -5.02  -6.91   121. 1999-01-02          6       -0.130
    ## # … with 26,288 more rows
    ## 
    ## $Al
    ## # A tibble: 26,298 x 9
    ##    beta_1 beta_2 beta_4 beta_6 beta_5 beta_3 date       factor_num factor_score
    ##     <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <date>          <dbl>        <dbl>
    ##  1   35.0  0.136   7.32   1.57   8.41   1.07 1999-01-01          1       NA    
    ##  2   35.0  0.136   7.32   1.57   8.41   1.07 1999-01-01          2       NA    
    ##  3   35.0  0.136   7.32   1.57   8.41   1.07 1999-01-01          4       NA    
    ##  4   35.0  0.136   7.32   1.57   8.41   1.07 1999-01-01          6       NA    
    ##  5   35.0  0.136   7.32   1.57   8.41   1.07 1999-01-01          5       NA    
    ##  6   35.0  0.136   7.32   1.57   8.41   1.07 1999-01-01          3       NA    
    ##  7   35.0  0.136   7.32   1.57   8.41   1.07 1999-01-02          1        0.413
    ##  8   35.0  0.136   7.32   1.57   8.41   1.07 1999-01-02          2        2.07 
    ##  9   35.0  0.136   7.32   1.57   8.41   1.07 1999-01-02          4        0.360
    ## 10   35.0  0.136   7.32   1.57   8.41   1.07 1999-01-02          6       -0.130
    ## # … with 26,288 more rows
    ## 
    ## $Si
    ## # A tibble: 26,298 x 9
    ##    beta_1 beta_2 beta_4 beta_6 beta_5 beta_3 date       factor_num factor_score
    ##     <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <date>          <dbl>        <dbl>
    ##  1   63.9   2.08  -2.56   10.8   1.53 -0.867 1999-01-01          1       NA    
    ##  2   63.9   2.08  -2.56   10.8   1.53 -0.867 1999-01-01          2       NA    
    ##  3   63.9   2.08  -2.56   10.8   1.53 -0.867 1999-01-01          4       NA    
    ##  4   63.9   2.08  -2.56   10.8   1.53 -0.867 1999-01-01          6       NA    
    ##  5   63.9   2.08  -2.56   10.8   1.53 -0.867 1999-01-01          5       NA    
    ##  6   63.9   2.08  -2.56   10.8   1.53 -0.867 1999-01-01          3       NA    
    ##  7   63.9   2.08  -2.56   10.8   1.53 -0.867 1999-01-02          1        0.413
    ##  8   63.9   2.08  -2.56   10.8   1.53 -0.867 1999-01-02          2        2.07 
    ##  9   63.9   2.08  -2.56   10.8   1.53 -0.867 1999-01-02          4        0.360
    ## 10   63.9   2.08  -2.56   10.8   1.53 -0.867 1999-01-02          6       -0.130
    ## # … with 26,288 more rows
    ## 
    ## $S
    ## # A tibble: 26,298 x 9
    ##    beta_1 beta_2 beta_4 beta_6 beta_5 beta_3 date       factor_num factor_score
    ##     <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <date>          <dbl>        <dbl>
    ##  1   196.   70.5   844.  -32.5   5.30   4.91 1999-01-01          1       NA    
    ##  2   196.   70.5   844.  -32.5   5.30   4.91 1999-01-01          2       NA    
    ##  3   196.   70.5   844.  -32.5   5.30   4.91 1999-01-01          4       NA    
    ##  4   196.   70.5   844.  -32.5   5.30   4.91 1999-01-01          6       NA    
    ##  5   196.   70.5   844.  -32.5   5.30   4.91 1999-01-01          5       NA    
    ##  6   196.   70.5   844.  -32.5   5.30   4.91 1999-01-01          3       NA    
    ##  7   196.   70.5   844.  -32.5   5.30   4.91 1999-01-02          1        0.413
    ##  8   196.   70.5   844.  -32.5   5.30   4.91 1999-01-02          2        2.07 
    ##  9   196.   70.5   844.  -32.5   5.30   4.91 1999-01-02          4        0.360
    ## 10   196.   70.5   844.  -32.5   5.30   4.91 1999-01-02          6       -0.130
    ## # … with 26,288 more rows
    ## 
    ## $Cl
    ## # A tibble: 26,298 x 9
    ##    beta_1 beta_2 beta_4 beta_6 beta_5 beta_3 date       factor_num factor_score
    ##     <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <date>          <dbl>        <dbl>
    ##  1  -3.07  -2.44  -19.3   3.86   7.78   69.5 1999-01-01          1       NA    
    ##  2  -3.07  -2.44  -19.3   3.86   7.78   69.5 1999-01-01          2       NA    
    ##  3  -3.07  -2.44  -19.3   3.86   7.78   69.5 1999-01-01          4       NA    
    ##  4  -3.07  -2.44  -19.3   3.86   7.78   69.5 1999-01-01          6       NA    
    ##  5  -3.07  -2.44  -19.3   3.86   7.78   69.5 1999-01-01          5       NA    
    ##  6  -3.07  -2.44  -19.3   3.86   7.78   69.5 1999-01-01          3       NA    
    ##  7  -3.07  -2.44  -19.3   3.86   7.78   69.5 1999-01-02          1        0.413
    ##  8  -3.07  -2.44  -19.3   3.86   7.78   69.5 1999-01-02          2        2.07 
    ##  9  -3.07  -2.44  -19.3   3.86   7.78   69.5 1999-01-02          4        0.360
    ## 10  -3.07  -2.44  -19.3   3.86   7.78   69.5 1999-01-02          6       -0.130
    ## # … with 26,288 more rows
    ## 
    ## $K
    ## # A tibble: 26,298 x 9
    ##    beta_1 beta_2 beta_4 beta_6 beta_5 beta_3 date       factor_num factor_score
    ##     <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <date>          <dbl>        <dbl>
    ##  1   8.82   4.06   7.74   3.39   7.15   2.35 1999-01-01          1       NA    
    ##  2   8.82   4.06   7.74   3.39   7.15   2.35 1999-01-01          2       NA    
    ##  3   8.82   4.06   7.74   3.39   7.15   2.35 1999-01-01          4       NA    
    ##  4   8.82   4.06   7.74   3.39   7.15   2.35 1999-01-01          6       NA    
    ##  5   8.82   4.06   7.74   3.39   7.15   2.35 1999-01-01          5       NA    
    ##  6   8.82   4.06   7.74   3.39   7.15   2.35 1999-01-01          3       NA    
    ##  7   8.82   4.06   7.74   3.39   7.15   2.35 1999-01-02          1        0.413
    ##  8   8.82   4.06   7.74   3.39   7.15   2.35 1999-01-02          2        2.07 
    ##  9   8.82   4.06   7.74   3.39   7.15   2.35 1999-01-02          4        0.360
    ## 10   8.82   4.06   7.74   3.39   7.15   2.35 1999-01-02          6       -0.130
    ## # … with 26,288 more rows
    ## 
    ## $Ca
    ## # A tibble: 26,298 x 9
    ##    beta_1 beta_2 beta_4 beta_6 beta_5 beta_3 date       factor_num factor_score
    ##     <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <date>          <dbl>        <dbl>
    ##  1   10.8   3.87  0.784   6.31   2.01   3.65 1999-01-01          1       NA    
    ##  2   10.8   3.87  0.784   6.31   2.01   3.65 1999-01-01          2       NA    
    ##  3   10.8   3.87  0.784   6.31   2.01   3.65 1999-01-01          4       NA    
    ##  4   10.8   3.87  0.784   6.31   2.01   3.65 1999-01-01          6       NA    
    ##  5   10.8   3.87  0.784   6.31   2.01   3.65 1999-01-01          5       NA    
    ##  6   10.8   3.87  0.784   6.31   2.01   3.65 1999-01-01          3       NA    
    ##  7   10.8   3.87  0.784   6.31   2.01   3.65 1999-01-02          1        0.413
    ##  8   10.8   3.87  0.784   6.31   2.01   3.65 1999-01-02          2        2.07 
    ##  9   10.8   3.87  0.784   6.31   2.01   3.65 1999-01-02          4        0.360
    ## 10   10.8   3.87  0.784   6.31   2.01   3.65 1999-01-02          6       -0.130
    ## # … with 26,288 more rows
    ## 
    ## $Ti
    ## # A tibble: 26,298 x 9
    ##    beta_1 beta_2 beta_4 beta_6 beta_5 beta_3 date       factor_num factor_score
    ##     <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <date>          <dbl>        <dbl>
    ##  1   1.56 0.0780  0.221  0.444   2.28 0.0428 1999-01-01          1       NA    
    ##  2   1.56 0.0780  0.221  0.444   2.28 0.0428 1999-01-01          2       NA    
    ##  3   1.56 0.0780  0.221  0.444   2.28 0.0428 1999-01-01          4       NA    
    ##  4   1.56 0.0780  0.221  0.444   2.28 0.0428 1999-01-01          6       NA    
    ##  5   1.56 0.0780  0.221  0.444   2.28 0.0428 1999-01-01          5       NA    
    ##  6   1.56 0.0780  0.221  0.444   2.28 0.0428 1999-01-01          3       NA    
    ##  7   1.56 0.0780  0.221  0.444   2.28 0.0428 1999-01-02          1        0.413
    ##  8   1.56 0.0780  0.221  0.444   2.28 0.0428 1999-01-02          2        2.07 
    ##  9   1.56 0.0780  0.221  0.444   2.28 0.0428 1999-01-02          4        0.360
    ## 10   1.56 0.0780  0.221  0.444   2.28 0.0428 1999-01-02          6       -0.130
    ## # … with 26,288 more rows
    ## 
    ## $V
    ## # A tibble: 26,298 x 9
    ##    beta_1 beta_2 beta_4 beta_6 beta_5 beta_3 date       factor_num factor_score
    ##     <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <date>          <dbl>        <dbl>
    ##  1  0.289   3.42  0.640  0.345  0.349 0.0679 1999-01-01          1       NA    
    ##  2  0.289   3.42  0.640  0.345  0.349 0.0679 1999-01-01          2       NA    
    ##  3  0.289   3.42  0.640  0.345  0.349 0.0679 1999-01-01          4       NA    
    ##  4  0.289   3.42  0.640  0.345  0.349 0.0679 1999-01-01          6       NA    
    ##  5  0.289   3.42  0.640  0.345  0.349 0.0679 1999-01-01          5       NA    
    ##  6  0.289   3.42  0.640  0.345  0.349 0.0679 1999-01-01          3       NA    
    ##  7  0.289   3.42  0.640  0.345  0.349 0.0679 1999-01-02          1        0.413
    ##  8  0.289   3.42  0.640  0.345  0.349 0.0679 1999-01-02          2        2.07 
    ##  9  0.289   3.42  0.640  0.345  0.349 0.0679 1999-01-02          4        0.360
    ## 10  0.289   3.42  0.640  0.345  0.349 0.0679 1999-01-02          6       -0.130
    ## # … with 26,288 more rows
    ## 
    ## $Cr
    ## # A tibble: 26,298 x 9
    ##    beta_1  beta_2 beta_4 beta_6 beta_5  beta_3 date       factor_num
    ##     <dbl>   <dbl>  <dbl>  <dbl>  <dbl>   <dbl> <date>          <dbl>
    ##  1 0.0280 0.00791 0.0387 0.0667  0.124 -0.0171 1999-01-01          1
    ##  2 0.0280 0.00791 0.0387 0.0667  0.124 -0.0171 1999-01-01          2
    ##  3 0.0280 0.00791 0.0387 0.0667  0.124 -0.0171 1999-01-01          4
    ##  4 0.0280 0.00791 0.0387 0.0667  0.124 -0.0171 1999-01-01          6
    ##  5 0.0280 0.00791 0.0387 0.0667  0.124 -0.0171 1999-01-01          5
    ##  6 0.0280 0.00791 0.0387 0.0667  0.124 -0.0171 1999-01-01          3
    ##  7 0.0280 0.00791 0.0387 0.0667  0.124 -0.0171 1999-01-02          1
    ##  8 0.0280 0.00791 0.0387 0.0667  0.124 -0.0171 1999-01-02          2
    ##  9 0.0280 0.00791 0.0387 0.0667  0.124 -0.0171 1999-01-02          4
    ## 10 0.0280 0.00791 0.0387 0.0667  0.124 -0.0171 1999-01-02          6
    ## # … with 26,288 more rows, and 1 more variable: factor_score <dbl>
    ## 
    ## $Mn
    ## # A tibble: 26,298 x 9
    ##    beta_1 beta_2 beta_4 beta_6 beta_5 beta_3 date       factor_num factor_score
    ##     <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <date>          <dbl>        <dbl>
    ##  1  0.172  0.207 0.0281   1.00 0.0646 0.0842 1999-01-01          1       NA    
    ##  2  0.172  0.207 0.0281   1.00 0.0646 0.0842 1999-01-01          2       NA    
    ##  3  0.172  0.207 0.0281   1.00 0.0646 0.0842 1999-01-01          4       NA    
    ##  4  0.172  0.207 0.0281   1.00 0.0646 0.0842 1999-01-01          6       NA    
    ##  5  0.172  0.207 0.0281   1.00 0.0646 0.0842 1999-01-01          5       NA    
    ##  6  0.172  0.207 0.0281   1.00 0.0646 0.0842 1999-01-01          3       NA    
    ##  7  0.172  0.207 0.0281   1.00 0.0646 0.0842 1999-01-02          1        0.413
    ##  8  0.172  0.207 0.0281   1.00 0.0646 0.0842 1999-01-02          2        2.07 
    ##  9  0.172  0.207 0.0281   1.00 0.0646 0.0842 1999-01-02          4        0.360
    ## 10  0.172  0.207 0.0281   1.00 0.0646 0.0842 1999-01-02          6       -0.130
    ## # … with 26,288 more rows
    ## 
    ## $Fe
    ## # A tibble: 26,298 x 9
    ##    beta_1 beta_2 beta_4 beta_6 beta_5 beta_3 date       factor_num factor_score
    ##     <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <date>          <dbl>        <dbl>
    ##  1   18.3   7.46   7.47   28.0   7.71 -0.992 1999-01-01          1       NA    
    ##  2   18.3   7.46   7.47   28.0   7.71 -0.992 1999-01-01          2       NA    
    ##  3   18.3   7.46   7.47   28.0   7.71 -0.992 1999-01-01          4       NA    
    ##  4   18.3   7.46   7.47   28.0   7.71 -0.992 1999-01-01          6       NA    
    ##  5   18.3   7.46   7.47   28.0   7.71 -0.992 1999-01-01          5       NA    
    ##  6   18.3   7.46   7.47   28.0   7.71 -0.992 1999-01-01          3       NA    
    ##  7   18.3   7.46   7.47   28.0   7.71 -0.992 1999-01-02          1        0.413
    ##  8   18.3   7.46   7.47   28.0   7.71 -0.992 1999-01-02          2        2.07 
    ##  9   18.3   7.46   7.47   28.0   7.71 -0.992 1999-01-02          4        0.360
    ## 10   18.3   7.46   7.47   28.0   7.71 -0.992 1999-01-02          6       -0.130
    ## # … with 26,288 more rows
    ## 
    ## $Ni
    ## # A tibble: 26,298 x 9
    ##    beta_1 beta_2 beta_4 beta_6 beta_5 beta_3 date       factor_num factor_score
    ##     <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <date>          <dbl>        <dbl>
    ##  1 0.0928   3.44  0.162  0.439  0.152 0.0904 1999-01-01          1       NA    
    ##  2 0.0928   3.44  0.162  0.439  0.152 0.0904 1999-01-01          2       NA    
    ##  3 0.0928   3.44  0.162  0.439  0.152 0.0904 1999-01-01          4       NA    
    ##  4 0.0928   3.44  0.162  0.439  0.152 0.0904 1999-01-01          6       NA    
    ##  5 0.0928   3.44  0.162  0.439  0.152 0.0904 1999-01-01          5       NA    
    ##  6 0.0928   3.44  0.162  0.439  0.152 0.0904 1999-01-01          3       NA    
    ##  7 0.0928   3.44  0.162  0.439  0.152 0.0904 1999-01-02          1        0.413
    ##  8 0.0928   3.44  0.162  0.439  0.152 0.0904 1999-01-02          2        2.07 
    ##  9 0.0928   3.44  0.162  0.439  0.152 0.0904 1999-01-02          4        0.360
    ## 10 0.0928   3.44  0.162  0.439  0.152 0.0904 1999-01-02          6       -0.130
    ## # … with 26,288 more rows
    ## 
    ## $Cu
    ## # A tibble: 26,298 x 9
    ##     beta_1 beta_2 beta_4 beta_6 beta_5 beta_3 date       factor_num factor_score
    ##      <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <date>          <dbl>        <dbl>
    ##  1 -0.0572  0.202  0.283  0.209   2.95  0.140 1999-01-01          1       NA    
    ##  2 -0.0572  0.202  0.283  0.209   2.95  0.140 1999-01-01          2       NA    
    ##  3 -0.0572  0.202  0.283  0.209   2.95  0.140 1999-01-01          4       NA    
    ##  4 -0.0572  0.202  0.283  0.209   2.95  0.140 1999-01-01          6       NA    
    ##  5 -0.0572  0.202  0.283  0.209   2.95  0.140 1999-01-01          5       NA    
    ##  6 -0.0572  0.202  0.283  0.209   2.95  0.140 1999-01-01          3       NA    
    ##  7 -0.0572  0.202  0.283  0.209   2.95  0.140 1999-01-02          1        0.413
    ##  8 -0.0572  0.202  0.283  0.209   2.95  0.140 1999-01-02          2        2.07 
    ##  9 -0.0572  0.202  0.283  0.209   2.95  0.140 1999-01-02          4        0.360
    ## 10 -0.0572  0.202  0.283  0.209   2.95  0.140 1999-01-02          6       -0.130
    ## # … with 26,288 more rows
    ## 
    ## $Zn
    ## # A tibble: 26,298 x 9
    ##    beta_1 beta_2 beta_4 beta_6 beta_5 beta_3 date       factor_num factor_score
    ##     <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <date>          <dbl>        <dbl>
    ##  1 0.0479   4.06   2.83   3.30   3.21 -0.421 1999-01-01          1       NA    
    ##  2 0.0479   4.06   2.83   3.30   3.21 -0.421 1999-01-01          2       NA    
    ##  3 0.0479   4.06   2.83   3.30   3.21 -0.421 1999-01-01          4       NA    
    ##  4 0.0479   4.06   2.83   3.30   3.21 -0.421 1999-01-01          6       NA    
    ##  5 0.0479   4.06   2.83   3.30   3.21 -0.421 1999-01-01          5       NA    
    ##  6 0.0479   4.06   2.83   3.30   3.21 -0.421 1999-01-01          3       NA    
    ##  7 0.0479   4.06   2.83   3.30   3.21 -0.421 1999-01-02          1        0.413
    ##  8 0.0479   4.06   2.83   3.30   3.21 -0.421 1999-01-02          2        2.07 
    ##  9 0.0479   4.06   2.83   3.30   3.21 -0.421 1999-01-02          4        0.360
    ## 10 0.0479   4.06   2.83   3.30   3.21 -0.421 1999-01-02          6       -0.130
    ## # … with 26,288 more rows
    ## 
    ## $Se
    ## # A tibble: 26,298 x 9
    ##     beta_1   beta_2 beta_4  beta_6 beta_5  beta_3 date       factor_num
    ##      <dbl>    <dbl>  <dbl>   <dbl>  <dbl>   <dbl> <date>          <dbl>
    ##  1 -0.0165 -0.00102  0.102 0.00811 0.0321 -0.0198 1999-01-01          1
    ##  2 -0.0165 -0.00102  0.102 0.00811 0.0321 -0.0198 1999-01-01          2
    ##  3 -0.0165 -0.00102  0.102 0.00811 0.0321 -0.0198 1999-01-01          4
    ##  4 -0.0165 -0.00102  0.102 0.00811 0.0321 -0.0198 1999-01-01          6
    ##  5 -0.0165 -0.00102  0.102 0.00811 0.0321 -0.0198 1999-01-01          5
    ##  6 -0.0165 -0.00102  0.102 0.00811 0.0321 -0.0198 1999-01-01          3
    ##  7 -0.0165 -0.00102  0.102 0.00811 0.0321 -0.0198 1999-01-02          1
    ##  8 -0.0165 -0.00102  0.102 0.00811 0.0321 -0.0198 1999-01-02          2
    ##  9 -0.0165 -0.00102  0.102 0.00811 0.0321 -0.0198 1999-01-02          4
    ## 10 -0.0165 -0.00102  0.102 0.00811 0.0321 -0.0198 1999-01-02          6
    ## # … with 26,288 more rows, and 1 more variable: factor_score <dbl>
    ## 
    ## $Br
    ## # A tibble: 26,298 x 9
    ##    beta_1 beta_2 beta_4 beta_6 beta_5 beta_3 date       factor_num factor_score
    ##     <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <date>          <dbl>        <dbl>
    ##  1 0.0519  0.206  0.450  0.131  0.313  0.160 1999-01-01          1       NA    
    ##  2 0.0519  0.206  0.450  0.131  0.313  0.160 1999-01-01          2       NA    
    ##  3 0.0519  0.206  0.450  0.131  0.313  0.160 1999-01-01          4       NA    
    ##  4 0.0519  0.206  0.450  0.131  0.313  0.160 1999-01-01          6       NA    
    ##  5 0.0519  0.206  0.450  0.131  0.313  0.160 1999-01-01          5       NA    
    ##  6 0.0519  0.206  0.450  0.131  0.313  0.160 1999-01-01          3       NA    
    ##  7 0.0519  0.206  0.450  0.131  0.313  0.160 1999-01-02          1        0.413
    ##  8 0.0519  0.206  0.450  0.131  0.313  0.160 1999-01-02          2        2.07 
    ##  9 0.0519  0.206  0.450  0.131  0.313  0.160 1999-01-02          4        0.360
    ## 10 0.0519  0.206  0.450  0.131  0.313  0.160 1999-01-02          6       -0.130
    ## # … with 26,288 more rows
    ## 
    ## $Ba
    ## # A tibble: 26,298 x 9
    ##    beta_1 beta_2 beta_4 beta_6 beta_5 beta_3 date       factor_num factor_score
    ##     <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <date>          <dbl>        <dbl>
    ##  1  0.442 0.0318  0.377 -0.149   6.17  0.321 1999-01-01          1       NA    
    ##  2  0.442 0.0318  0.377 -0.149   6.17  0.321 1999-01-01          2       NA    
    ##  3  0.442 0.0318  0.377 -0.149   6.17  0.321 1999-01-01          4       NA    
    ##  4  0.442 0.0318  0.377 -0.149   6.17  0.321 1999-01-01          6       NA    
    ##  5  0.442 0.0318  0.377 -0.149   6.17  0.321 1999-01-01          5       NA    
    ##  6  0.442 0.0318  0.377 -0.149   6.17  0.321 1999-01-01          3       NA    
    ##  7  0.442 0.0318  0.377 -0.149   6.17  0.321 1999-01-02          1        0.413
    ##  8  0.442 0.0318  0.377 -0.149   6.17  0.321 1999-01-02          2        2.07 
    ##  9  0.442 0.0318  0.377 -0.149   6.17  0.321 1999-01-02          4        0.360
    ## 10  0.442 0.0318  0.377 -0.149   6.17  0.321 1999-01-02          6       -0.130
    ## # … with 26,288 more rows
    ## 
    ## $Pb
    ## # A tibble: 26,298 x 9
    ##    beta_1 beta_2 beta_4 beta_6 beta_5  beta_3 date       factor_num factor_score
    ##     <dbl>  <dbl>  <dbl>  <dbl>  <dbl>   <dbl> <date>          <dbl>        <dbl>
    ##  1  0.313  0.791  0.935  0.394   1.69 0.00971 1999-01-01          1       NA    
    ##  2  0.313  0.791  0.935  0.394   1.69 0.00971 1999-01-01          2       NA    
    ##  3  0.313  0.791  0.935  0.394   1.69 0.00971 1999-01-01          4       NA    
    ##  4  0.313  0.791  0.935  0.394   1.69 0.00971 1999-01-01          6       NA    
    ##  5  0.313  0.791  0.935  0.394   1.69 0.00971 1999-01-01          5       NA    
    ##  6  0.313  0.791  0.935  0.394   1.69 0.00971 1999-01-01          3       NA    
    ##  7  0.313  0.791  0.935  0.394   1.69 0.00971 1999-01-02          1        0.413
    ##  8  0.313  0.791  0.935  0.394   1.69 0.00971 1999-01-02          2        2.07 
    ##  9  0.313  0.791  0.935  0.394   1.69 0.00971 1999-01-02          4        0.360
    ## 10  0.313  0.791  0.935  0.394   1.69 0.00971 1999-01-02          6       -0.130
    ## # … with 26,288 more rows

``` r
sources1 <- bind_rows(sources, .id = "element") %>%
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

# # need to figure out how to iterate over the following:
# sources1 <- sources[["Al"]] %>% 
#   pivot_longer(starts_with("beta"),
#                names_to = "beta_num",
#                names_prefix = "beta_",
#                values_to = "beta_value") %>% 
#   mutate(
#     beta_num = as.numeric(beta_num)
#   ) %>% 
#   filter(factor_num == beta_num) %>% 
#   mutate(
#     srce = factor_score*beta_value
#   ) %>% 
#   mutate(srce_num = as.character(beta_num)) %>% 
#   select(date, srce_num, srce) 


# need to add back r.squared value here

conc2 <- conc %>% 
  pivot_longer(pm25:Pb,
               names_to = "element",
               values_to = "concentration") %>% 
  arrange(element)

mean_source <- left_join(sources1, conc2, by = c('element', 'date')) %>% 
  group_by(element, srce_num) %>% 
  summarize(source_mean = mean(srce, na.rm = TRUE),
            MeanConc = mean(concentration, na.rm = TRUE)) %>% 
  pivot_wider(names_from = "srce_num",
              names_prefix = "source_",
              values_from = "source_mean")
```

    ## `summarise()` regrouping output by 'element' (override with `.groups` argument)

``` r
# mean_source <- full_join(sources1, conc, by = "date") %>% 
#   group_by(element, srce_num) %>%
#   summarize(source_mean = mean(srce, na.rm = TRUE),
#             MeanConc = mean(Al, na.rm = TRUE)) %>% #also here!
#   pivot_wider(srce_num:MeanConc,
#               names_from = "srce_num",
#               names_prefix = "source_",
#               values_from = "source_mean")

SA_Varimax <- mean_source %>% 
  ungroup(element) %>% 
  mutate(
    PredConc = reduce(select(., contains("source")), `+`),
    Pct_error = (PredConc - MeanConc)/MeanConc
  )
```

also extract loadings

PM2.5 regression gives contribution of each source to total pm2.5

Do general time series of PM2.5 and MI Total PM2.5 and MI APCA, PCP, PCP
+ APCA? when running APCA for NYC, extract loadings and for each source
contribution look at seasonal patterns (break calendar into 4 seasons of
3 months each), look at concentrations/contributions by season (traffic
higher in winter, regional in winter), day of week (traffic higher on
weekdays) scores: correlations between sources and PM2.5 + patterns
