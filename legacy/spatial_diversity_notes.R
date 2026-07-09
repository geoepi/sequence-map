# choose subset of PCs for "diversity"
K <- 5
pc_cols <- paste0("PC_", 1:K)

prov_df <- sea_data_in %>%
  group_by(province, X, Y) %>%  # X,Y are province centroids
  summarise(
    n = n(),
    PC1_mean = mean(PC_1, na.rm = TRUE), # individual PCs
    PC2_mean = mean(PC_2, na.rm = TRUE),
    PC3_mean = mean(PC_3, na.rm = TRUE),
    PC4_mean = mean(PC_4, na.rm = TRUE),
    PC5_mean = mean(PC_5, na.rm = TRUE),
    # within-province dispersion across PCs 1..K
    # div_trace = province diversity
    div_trace = {
      Z <- as.matrix(dplyr::pick(dplyr::all_of(pc_cols)))
      Z <- Z[complete.cases(Z), , drop = FALSE]
      if (nrow(Z) >= 2) sum(diag(stats::cov(Z))) else NA_real_
    },
    year_med = stats::median(year, na.rm = TRUE), # median year/dominant host
    host_mode = {
      h <- host[!is.na(host)]
      if (length(h) == 0) NA_character_ else names(which.max(table(h)))
    },
    .groups = "drop"
  ) %>%
  mutate(
    host_mode = factor(host_mode),
    log_div_trace = log(div_trace)  # log scale
  )


locs <- as.matrix(prov_df[, c("X", "Y")])

spde <- inla.spde2.pcmatern(
  mesh_dom, alpha = 2,
  prior.range = c(diff(range(prov_df$Y))/3, 0.5),  # P(range < r0)=p
  prior.sigma = c(2, 0.01),                        # relax vs your earlier sigma=1
  constr = TRUE
)

srf_idx <- inla.spde.make.index("srf_idx", spde$n.spde)
A.mat  <- inla.spde.make.A(mesh_dom, loc = locs)



fmdv.lst = list(c(srf_idx, # index for mesh
                  list(intercept1 = 1)), # custom intercept (see other example)
                list(year = prov_df[,"year_med"],
                     host = prov_df[,"host_mode"])) 

stk_pc1 <- inla.stack(data = list(Y = prov_df$PC1_mean), # response/dependent variable
                       A = list(A.mat, 1), # matrix to align with mesh
                       effects = fmdv.lst, # variables of interest  
                       tag = "base.0") # arbitrary name/label

# weakly-informative priors for iid effects
pcprec <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))

form_pc1 <- Y ~ -1 + intercept1 +
  f(srf_idx, model = spde) +
  f(year_med, model = "iid", constr = TRUE, hyper = pcprec) +
  f(host_mode, model = "iid", constr = TRUE, hyper = pcprec)

m_pc1 <- inla(
  form_pc1,
  data = inla.stack.data(stk_pc1),
  family = "gaussian",
  verbose = TRUE,
  control.predictor = list(A = inla.stack.A(stk_pc1), compute = TRUE),
  control.compute = list(waic = TRUE, config=TRUE)
)


#############  Diversity
prov_div <- prov_df %>% filter(is.finite(log_div_trace))

locs_d <- as.matrix(prov_div[, c("X", "Y")])
A_d    <- inla.spde.make.A(mesh_dom, loc = locs_d)
srf_idx <- inla.spde.make.index("srf_idx", spde$n.spde)



stk_div <- inla.stack(
  data = list(y = prov_div$log_div_trace),
  A = list(A_d, 1),
  effects = list(
    idx,
    data.frame(intercept = 1,
               year_med = prov_div$year_med,
               host_mode = prov_div$host_mode)
  ),
  tag = "div"
)

div.lst = list(c(srf_idx, # index for mesh
                  list(intercept1 = 1)), # custom intercept (see other example)
                list(year = prov_div[,"year_med"],
                     host = prov_div[,"host_mode"])) 

stk_div <- inla.stack(data = list(Y = prov_div$log_div_trace), # response/dependent variable
                      A = list(A_d , 1), # matrix to align with mesh
                      effects = div.lst, # variables of interest  
                      tag = "div.0") # arbitrary name/label


form_div <- Y ~ -1 + intercept1 +
  f(srf_idx, model = spde) +
  f(year_med, model = "iid", constr = TRUE, hyper = pcprec) +
  f(host_mode, model = "iid", constr = TRUE, hyper = pcprec)

m_div <- inla(
  form_div,
  data = inla.stack.data(stk_div),
  family = "gaussian",
  verbose=TRUE,
  control.predictor = list(A = inla.stack.A(stk_div), compute = TRUE),
  control.compute = list(waic = TRUE, config=TRUE)
)


##########  Rasters
# vn_rast defines the raster template you already have
grid_xy <- xyFromCell(vn_rast, 1:ncell(vn_rast))
Ap_grid <- inla.spde.make.A(mesh_dom, loc = grid_xy)

# Extract intercept and spatial field mean at mesh nodes
beta0_pc1 <- m_pc1$summary.fixed["intercept", "mean"]
u_pc1     <- m_pc1$summary.random$srf$mean

pc1_grid_mean <- beta0_pc1 + drop(Ap_grid %*% u_pc1)

PC1_r_mean <- setValues(vn_rast, pc1_grid_mean)
PC1_r_mean <- PC1_r_mean + vn_rast

beta0_div <- m_div$summary.fixed["intercept", "mean"]
u_div     <- m_div$summary.random$srf$mean

logdiv_grid_mean <- beta0_div + drop(Ap_grid %*% u_div)
div_grid_mean <- exp(logdiv_grid_mean)

DIV_r_mean <- setValues(vn_rast, div_grid_mean)
DIV_r_mean <- DIV_r_mean + vn_rast
plot(DIV_r_mean)

##### sample
# Draw from the latent Gaussian field (includes fixed + random effects)
# This returns samples for the whole latent vector; easiest is to draw the field and intercept.
nsamp <- 1
samp_pc1 <- inla.posterior.sample(nsamp, m_pc1)[[1]]

# Helper: find indices in the latent vector
# Fixed effects are named; random field entries are in summary.random with known names.
# Robust approach: use inla.stack and internal ordering is painful; simplest is to draw using inla.posterior.sample
# and then reconstruct via the marginal draws for fixed+random:
# practical workaround: draw the spatial field from its marginals (approx) using mean/sd
# (If you need exact joint draws, use INLA's excursion/sample utilities; this is the pragmatic version.)

u_pc1_draw <- rnorm(length(u_pc1),
                    mean = m_pc1$summary.random$srf$mean,
                    sd   = m_pc1$summary.random$srf$sd)
beta0_draw <- rnorm(1,
                    mean = m_pc1$summary.fixed["intercept", "mean"],
                    sd   = m_pc1$summary.fixed["intercept", "sd"])

pc1_grid_draw <- beta0_draw + drop(Ap_grid %*% u_pc1_draw)
PC1_r_draw <- setValues(vn_rast, pc1_grid_draw)


