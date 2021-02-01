APCA NYC
================
Rachel Tao
1/23/2021

Load in data

``` r
rename_pollutants = function(df){
  data = df %>% 
  rename(
    date = date_local,
    aluminum = aluminum_pm2_5_lc_88104, 
    arsenic = arsenic_pm2_5_lc_88103,
    barium = barium_pm2_5_lc_88107,
    bromine = bromine_pm2_5_lc_88109,
    cadmium = cadmium_pm2_5_lc_88110,
    calcium = calcium_pm2_5_lc_88111,
    chromium = chromium_pm2_5_lc_88112,
    copper = copper_pm2_5_lc_88114,
    iron = iron_pm2_5_lc_88126,
    lead = lead_pm2_5_lc_88128,
    magnesium = magnesium_pm2_5_lc_88140,
    manganese = manganese_pm2_5_lc_88132,
    potassium = potassium_pm2_5_lc_88180,
    nickel = nickel_pm2_5_lc_88136,
    selenium = selenium_pm2_5_lc_88154,
    silicon = silicon_pm2_5_lc_88165,
    sulfate = sulfate_pm2_5_lc_88403,
    sulfur = sulfur_pm2_5_lc_88169,
    titanium = titanium_pm2_5_lc_88161,
    nitrates = total_nitrate_pm2_5_lc_88306,
    vanadium = vanadium_pm2_5_lc_88164,
    zinc = zinc_pm2_5_lc_88167
    )
}

pm25_tot <- read_csv("./data/pm2.5_total_components_daily_nyc_avg.csv") %>% 
  rename(
    date = date_local,
    pm25 = pm2.5_total
  ) %>% 
  select(date, pm25)

conc <- read_csv("./data/NYC_pm2.5_components_daily.csv") %>% 
  rename_pollutants() %>% 
  select(-X1) %>% 
  left_join(pm25_tot, by = "date") %>% 
  relocate(date, pm25)

# multiply all concentrations by 1000 for units
i <- c(2:24)
conc[ , i] <- apply(conc[ , i], 2,
                    function(x)
                      x*1000)
```

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

vars <- conc_summarize %>% 
  summarize(variance = var(concentration, na.rm = TRUE)) %>% 
  arrange(element)
```

Scaled factor analysis - must provide number of factors here

``` r
# set mean = 0, SD = 1
zspec_pm25 <- 
  conc %>% 
  select(-date, -pm25, -sulfur, -potassium) %>% 
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
y <- conc %>% select(-date, -sulfur)

# repeat 'regress' for each element using 'x' as factors and 'y' as concentrations)
regress <- map(y, regress) %>% 
  enframe(name = "element") %>% 
  unnest(c(value)) %>% 
  janitor::clean_names() %>% 
  rename_at(vars(contains("factor")), funs(str_replace(., "xfactor", "beta")))

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

# calculate mean source contributions for each element
mean_source <- source_conc %>% 
  group_by(element, srce_num) %>% 
  summarize(source_mean = mean(srce, na.rm = TRUE),
            MeanConc = mean(concentration, na.rm = TRUE)) %>% 
  pivot_wider(names_from = "srce_num",
              names_prefix = "source_",
              values_from = "source_mean")
```

### Source contributions to each component of PM2.5

``` r
# calculate mean predicted concentration and percent error
SA_Varimax <- mean_source %>% 
  ungroup(element) %>% 
  left_join(r, by = "element") %>% 
  mutate(
    PredConc = reduce(select(., contains("source")), `+`),
    Pct_error = round(((PredConc - MeanConc)/MeanConc)*100, digits = 2)
  )

i <- c(2:11)
SA_Varimax[ , i] <- apply(SA_Varimax[ , i], 2,
                    function(x)
                      round(x, digits = 3))

SA_Varimax %>% kable()
```

| element   |  MeanConc | source\_1 | source\_2 | source\_3 | source\_4 | source\_5 | source\_6 | r\_squared | PredConc | Pct\_error |
| :-------- | --------: | --------: | --------: | --------: | --------: | --------: | --------: | ---------: | -------: | ---------: |
| aluminum  |    20.906 |   \-0.191 |    20.211 |   \-1.471 |     3.126 |     1.879 |     1.235 |      0.619 |   24.789 |      18.58 |
| arsenic   |     0.657 |     0.278 |     0.033 |     0.041 |     0.150 |     0.032 |     0.012 |      0.330 |    0.545 |    \-17.03 |
| barium    |     6.227 |     0.562 |     0.885 |     3.222 |     1.814 |   \-1.115 |     0.960 |      0.917 |    6.327 |       1.61 |
| bromine   |     3.312 |     1.274 |     0.114 |   \-0.021 |     1.257 |   \-0.099 |     0.196 |      0.589 |    2.721 |    \-17.84 |
| cadmium   |     1.510 |   \-0.092 |     0.081 |   \-0.069 |     0.149 |   \-0.080 |   \-0.009 |      0.011 |  \-0.020 |   \-101.34 |
| calcium   |    56.074 |     1.747 |    19.033 |   \-0.663 |    32.400 |   \-2.785 |   \-1.174 |      0.711 |   48.560 |    \-13.40 |
| chromium  |     1.628 |   \-0.196 |   \-0.111 |   \-0.036 |   \-0.421 |     3.338 |   \-0.234 |      0.294 |    2.339 |      43.70 |
| copper    |     4.483 |     0.580 |     0.455 |     0.058 |     2.109 |     1.393 |     0.447 |      0.699 |    5.043 |      12.48 |
| iron      |   111.841 |    13.293 |    28.919 |     1.949 |    25.810 |    36.427 |     1.433 |      0.960 |  107.831 |     \-3.59 |
| lead      |     3.293 |     1.344 |     0.396 |     0.359 |     2.260 |   \-0.367 |     0.067 |      0.629 |    4.059 |      23.24 |
| magnesium |     7.682 |     0.061 |     0.680 |     0.452 |     0.371 |   \-0.730 |     1.481 |      0.125 |    2.316 |    \-69.86 |
| manganese |     2.373 |     0.407 |     0.688 |   \-0.019 |     1.513 |     0.363 |   \-0.080 |      0.685 |    2.871 |      21.00 |
| nickel    |     8.590 |     2.605 |     0.533 |     1.188 |     4.985 |     0.001 |   \-1.676 |      0.934 |    7.637 |    \-11.10 |
| nitrates  |  1811.742 |   889.636 |  \-63.063 |  \-49.501 |   962.873 |   122.526 |  \-51.797 |      0.619 | 1810.675 |     \-0.06 |
| pm25      | 11855.453 |  4860.086 |  1824.951 |     7.394 |  1718.873 |   484.892 |   366.934 |      0.766 | 9263.131 |    \-21.87 |
| potassium |    47.913 |     8.508 |     9.661 |     2.338 |    12.887 |   \-0.146 |     6.067 |      0.214 |   39.315 |    \-17.95 |
| selenium  |     0.755 |     0.679 |     0.128 |     0.048 |     0.139 |   \-0.050 |   \-0.018 |      0.773 |    0.926 |      22.59 |
| silicon   |    72.073 |     9.324 |    59.260 |     0.364 |     5.147 |   \-7.317 |     1.260 |      0.913 |   68.038 |     \-5.60 |
| sulfate   |  3008.629 |  1860.004 |   706.139 |    33.329 | \-313.679 |  \-17.717 |   148.965 |      0.756 | 2417.041 |    \-19.66 |
| titanium  |     3.445 |     0.520 |     2.146 |     0.322 |     0.508 |     0.004 |     0.023 |      0.691 |    3.522 |       2.26 |
| vanadium  |     4.449 |     2.122 |     1.088 |     0.277 |     1.975 |     0.142 |   \-0.253 |      0.708 |    5.350 |      20.25 |
| zinc      |    26.876 |     3.272 |   \-1.006 |   \-0.329 |    25.477 |   \-2.426 |     0.701 |      0.900 |   25.690 |     \-4.41 |

``` r
# calculate proportion of each element coming from each source 
# (negative values from SA_Varimax converted to 0%)
SA_proportion_long <- source_conc %>% 
  group_by(element, srce_num) %>% 
  summarize(source_mean = mean(srce, na.rm = TRUE),
            MeanConc = mean(concentration, na.rm = TRUE)) %>% 
  mutate(
    source_mean = if_else(source_mean < 0, 0, source_mean)) %>% 
  ungroup() %>% 
  group_by(element) %>% 
  mutate(
    total = sum(source_mean),
    source_percent = round((source_mean/total)*100, digits = 2))

# make the table nicer to look at
SA_proportion <- SA_proportion_long %>% 
  relocate(element, srce_num, source_percent) %>% 
  pivot_wider(element: source_percent,
              names_from = srce_num,
              names_prefix = "source_",
              values_from = source_percent) %>% 
  kable()

SA_proportion
```

| element   | source\_1 | source\_2 | source\_3 | source\_4 | source\_5 | source\_6 |
| :-------- | --------: | --------: | --------: | --------: | --------: | --------: |
| aluminum  |      0.00 |     76.41 |      0.00 |     11.82 |      7.10 |      4.67 |
| arsenic   |     50.97 |      6.08 |      7.44 |     27.50 |      5.81 |      2.20 |
| barium    |      7.56 |     11.89 |     43.29 |     24.37 |      0.00 |     12.90 |
| bromine   |     44.84 |      4.00 |      0.00 |     44.25 |      0.00 |      6.91 |
| cadmium   |      0.00 |     35.11 |      0.00 |     64.89 |      0.00 |      0.00 |
| calcium   |      3.29 |     35.79 |      0.00 |     60.92 |      0.00 |      0.00 |
| chromium  |      0.00 |      0.00 |      0.00 |      0.00 |    100.00 |      0.00 |
| copper    |     11.51 |      9.03 |      1.15 |     41.82 |     27.63 |      8.86 |
| iron      |     12.33 |     26.82 |      1.81 |     23.94 |     33.78 |      1.33 |
| lead      |     30.36 |      8.95 |      8.11 |     51.07 |      0.00 |      1.52 |
| magnesium |      2.01 |     22.32 |     14.85 |     12.19 |      0.00 |     48.63 |
| manganese |     13.69 |     23.16 |      0.00 |     50.92 |     12.23 |      0.00 |
| nickel    |     27.97 |      5.73 |     12.76 |     53.53 |      0.01 |      0.00 |
| nitrates  |     45.04 |      0.00 |      0.00 |     48.75 |      6.20 |      0.00 |
| pm25      |     52.47 |     19.70 |      0.08 |     18.56 |      5.23 |      3.96 |
| potassium |     21.56 |     24.48 |      5.93 |     32.66 |      0.00 |     15.37 |
| selenium  |     68.30 |     12.92 |      4.80 |     13.98 |      0.00 |      0.00 |
| silicon   |     12.37 |     78.64 |      0.48 |      6.83 |      0.00 |      1.67 |
| sulfate   |     67.67 |     25.69 |      1.21 |      0.00 |      0.00 |      5.42 |
| titanium  |     14.77 |     60.91 |      9.14 |     14.42 |      0.10 |      0.66 |
| vanadium  |     37.87 |     19.42 |      4.94 |     35.24 |      2.53 |      0.00 |
| zinc      |     11.11 |      0.00 |      0.00 |     86.51 |      0.00 |      2.38 |

``` r
# bar graph showing proportion of each element coming from each source
SA_bar <- SA_proportion_long %>% 
  ungroup() %>% 
  mutate(element = fct_reorder(element, source_percent)) %>% 
  ggplot(aes(x = element, y = source_mean, fill = srce_num)) +
  geom_col(position = "fill") +
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("") +
  xlab("")

SA_bar
```

<img src="APCA_nyc_files/figure-gfm/unnamed-chunk-8-1.png" width="90%" />

``` r
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
