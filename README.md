# Manually calculate Durbin impacts
Van Kinh Nguyen

``` r
library(TMB)
library(spdep)
library(Matrix)
library(dplyr)
library(sf)
set.seed(123)
```

Data Simulation for a model

$$
y = \rho W y + \beta_1 x_1 + \beta_2 x_2 + \beta_3 x_3 + \theta_1 W x_1 + \theta_2 W x_2
$$

Note that there is no $W x_3$.

``` r
# True parameter values
true_rho <- 0.6
true_beta <- c(5, 1.5, -1, 0.8) # Intercept, b1, b2, b3
true_theta <- c(0.75, -0.5)   # t1, t2
true_sigma <- 1.0

# simulate a spatial grid
n_dim <- 20
n_obs <- n_dim * n_dim 
bbox <- st_bbox(c(xmin = 0, ymin = 0, xmax = n_dim, ymax = n_dim))
grid <- st_make_grid(bbox, n = n_dim, square = TRUE)
nb <- poly2nb(grid, queen = TRUE)
W_listw <- nb2listw(nb, style = "W", zero.policy = TRUE)
W_dense <- spdep::listw2mat(W_listw)
W <- as(W_dense, "sparseMatrix")

# Generate covariates
x1 <- rnorm(n_obs)
x2 <- rnorm(n_obs)
x3 <- rnorm(n_obs)

# Create model matrices
X <- model.matrix(~ x1 + x2 + x3)
WX_sub <- as.matrix(W %*% X[, c("x1", "x2")]) # Lag only x1 and x2
# Generate the response variable y
I <- Diagonal(n_obs)
A <- I - true_rho * W
mu_non_spatial <- X %*% true_beta + WX_sub %*% true_theta
epsilon <- rnorm(n_obs, mean = 0, sd = true_sigma)
y <- as.vector(solve(A, mu_non_spatial + epsilon))

# simulate data frame 
sim_data <- st_sf(data.frame(y, x1, x2, x3), geometry = st_geometry(grid))
```

## TMB Implementation

``` r
compile("SDLM.cpp", flags = '-Wno-Wunused-variable')
```

``` r
dyn.load(dynlib("SDLM"))

# Prepare data and parameters for TMB
data_list <- list(
  y = y,
  X = X,
  WX = WX_sub,
  W = W,
  W_eigs = eigen(W)$values
)

params_list <- list(
  beta = rep(0, ncol(X)),
  theta = rep(0, ncol(WX_sub)),
  log_sigma = 0,
  atanh_rho = 0
)

obj <- MakeADFun(data_list, params_list, DLL = "SDLM", silent = TRUE)

fit_tmb <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr)
sdr_tmb <- sdreport(obj)
```

## lagsarlm with Durbin argument

``` r
library(spatialreg)
fit_sarlm <- lagsarlm(
  formula = y ~ x1 + x2 + x3,
  data = sim_data,
  listw = W_listw,
  Durbin = ~ x1 + x2, 
  method = "eigen"
)
```

## lagsarlm without Durbin argument

Manually create lag with `create_WX` (making `impacts` not correct if
apply directly).

``` r
lag_data <- sim_data %>% 
  st_drop_geometry() %>% 
  select(x1, x2) %>% 
  create_WX(listw = W_listw, prefix = 'lag') %>% 
  cbind(sim_data)

fit_manual <- lagsarlm(
  formula = y ~ x1 + x2 + x3 + lag.x1 + lag.x2,
  data = lag_data,
  listw = W_listw,
  Durbin = FALSE,
  method = "eigen"
)
```

## Comparision of coefficients

TMB’s

``` r
sdr_tmb
```

    sdreport(.) result
                 Estimate Std. Error
    beta       5.37400064 0.57346664
    beta       1.46935194 0.05101161
    beta      -0.98547741 0.04768317
    beta       0.94039157 0.04597507
    theta      0.88525331 0.18727521
    theta     -0.21943856 0.15709764
    log_sigma -0.06490111 0.03574329
    atanh_rho  0.64895303 0.06705832

with Durbin args

``` r
summary(fit_sarlm)
```
                 Estimate Std. Error  z value  Pr(>|z|)
    (Intercept)  5.374016   0.583322   9.2128 < 2.2e-16
    x1           1.469352   0.050855  28.8930 < 2.2e-16
    x2          -0.985477   0.047741 -20.6421 < 2.2e-16
    x3           0.940391   0.046047  20.4224 < 2.2e-16
    lag.x1       0.885249   0.182901   4.8401 1.298e-06
    lag.x2      -0.219443   0.158969  -1.3804    0.1675

    Rho: 0.57096, LR test value: 120.33, p-value: < 2.22e-16

without Durbin args

``` r
summary(fit_manual)
```
                 Estimate Std. Error  z value  Pr(>|z|)
    (Intercept)  5.374016   0.583322   9.2128 < 2.2e-16
    x1           1.469352   0.050855  28.8930 < 2.2e-16
    x2          -0.985477   0.047741 -20.6421 < 2.2e-16
    x3           0.940391   0.046047  20.4224 < 2.2e-16
    lag.x1       0.885249   0.182901   4.8401 1.298e-06
    lag.x2      -0.219443   0.158969  -1.3804    0.1675

    Rho: 0.57096, LR test value: 120.33, p-value: < 2.22e-16
## Compare impacts estimates

First, load the function for manual impact calculation

``` r
impact_manual <- function(lagmodel, fix_names, lag_names, W) {
  betas <- coef(fit_manual)[fix_names]
  thetas <- coef(fit_manual)[lag_names]
  names(thetas) <- gsub("lag\\.", "", names(thetas))
  rho <- lagmodel$rho
  S_rho <- solve(I - rho * W)
  impacts_list <- list()
  for (var in names(betas)) {
    beta_k <- betas[var]
    theta_k <- thetas[var]
    if (is.null(theta_k) | is.na(theta_k)) theta_k = 0
    M <- S_rho %*% (I * beta_k + W * theta_k)
    direct_impact <- mean(diag(M))
    indirect_impact <- total_impact - direct_impact
    impacts_list[[var]] <- c(
      direct = direct_impact,
      indirect = indirect_impact,
      total = total_impact)
  }
  do.call(rbind, impacts_list)
}
```

Impact with Durbin args

``` r
impacts(fit_sarlm, listw = W_listw) 
```

    Impact measures (mixed, exact):
          Direct  Indirect     Total
    x1  1.669649  3.818466  5.488115
    x2 -1.076334 -1.732098 -2.808432
    x3  1.002765  1.189102  2.191867

Impact without Durbin args, manual lag-covariates

``` r
fix_names = c('x1', 'x2', 'x3')
lag_names = c('lag.x1', 'lag.x2')
impact_manual(fit_manual, fix_names, lag_names, W)
```

          direct  indirect     total
    x1  1.669649  3.818466  5.488115
    x2 -1.076334 -1.732098 -2.808432
    x3  1.002765  1.189102  2.191867
