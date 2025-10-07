###########################Simulation CODE###########################
run_iv_sim <- function(
    reps, ## number of replication
    n,    ## sample size
    alpha.true, beta.true, iota.true, theta.true, zeta.true,
    scenario = c("M1","M2","M3"),               
    link_delta = "tanh",
    z.type = "binary",
    quantile.z = 0.5,
    ncores = max(1, parallel::detectCores() - 1), # choose 1 for sequential
    seed = 1234
) {
  build_scenario_lib <- function(dat) {
    mm <- function(vars) as.matrix(dat[, vars, drop = FALSE])
    list(
      correct = list(
        delta   = mm(c("Intercept","X1","X2")),
        delta_z = mm(c("Intercept","X1","X2")),
        mu_d    = mm(c("Intercept","X1","X2")),
        mu_y    = mm(c("Intercept","X1","X2")),
        pi      = mm(c("Intercept","X1","X2"))
      ),
      wrong = list(
        delta   = mm(c("Intercept","X1")),  
        delta_z = mm(c("Intercept","X1")),   
        mu_d    = mm(c("Intercept","X1")), 
        mu_y    = mm(c("Intercept","X1")),  
        pi      = mm(c("Intercept","X1"))   
      )
    )
  }
  
  pick_scenario <- function(dat, scenario) {
    lib <- build_scenario_lib(dat)
    spec_map <- switch(scenario,
                       allcorrect = c(delta="correct", mu_d="correct", mu_y="correct",
                                      pi="correct", delta_z="correct"),
                       M1 = c(delta="correct", mu_d="correct", mu_y="correct",
                              pi="wrong", delta_z="wrong"),
                       M2 = c(pi="correct", mu_d="correct", delta_z="correct",
                              delta="wrong", mu_y="wrong"),
                       M3 = c(delta="correct", pi="correct",
                              mu_d="wrong", mu_y="wrong", delta_z="wrong"),
                       allwrong = c(delta="wrong", pi="wrong",
                                    mu_d="wrong", mu_y="wrong", delta_z="wrong"),
                       stop("scenario must be one of 'allcorrect', 'M1','M2','M3', 'allwrong'")
    )
    list(
      x.delta    = lib[[ spec_map[["delta"]]   ]][["delta"]],
      x.delta.z  = lib[[ spec_map[["delta_z"]] ]][["delta_z"]],
      x.mu.d     = lib[[ spec_map[["mu_d"]]    ]][["mu_d"]],
      x.mu.y     = lib[[ spec_map[["mu_y"]]    ]][["mu_y"]],
      x.pi       = lib[[ spec_map[["pi"]]      ]][["pi"]]
    )
  }
  
  # --- cluster setup (auto-cleanup) -----------------------------------------
  cl <- NULL; created <- FALSE
  if (ncores > 1) {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    created <- TRUE
    on.exit({ if (created) parallel::stopCluster(cl) }, add = TRUE)
    
    parallel::clusterEvalQ(cl, {
      source("fun_generate_data.R")
      source("fun_estimation.R")
      NULL
    })
    # NEW: export scalars used inside foreach
    parallel::clusterExport(
      cl,
      varlist = c("n","alpha.true","beta.true","iota.true","theta.true","zeta.true","link_delta"),
      envir = environment()
    )
  } else {
    doParallel::registerDoParallel(cores = 1)  # no clusterEvalQ here when cl is NULL
  }
  
  doRNG::registerDoRNG(seed)
  
  # --- one scenario ---------------------------------------------------
  run_one_scenario <- function(scn) {
    res <- foreach(i=1:reps, .combine='rbind',
                   .export=c("data_gen", "estimation", "build_scenario_lib", "pick_scenario",
                             "n", "alpha.true", "beta.true", "iota.true", "theta.true", 
                             "zeta.true", "z.type", "quantile.z"),
                   .options.RNG=seed) %dorng% {
                     
                    if(z.type == "binary"){
                      data <- data_gen_binZ(alpha.true, beta.true, iota.true, theta.true, zeta.true,
                                            n, Y.type="binary", seed=i)
                      ATE.crude <- with(data$dat, mean(Y[D==1]) - mean(Y[D==0]))
                      
                      des <- pick_scenario(data$dat, scn)
                      
                      est <- estimation_binZ(x.delta   = des$x.delta,
                                             x.delta.z = des$x.delta.z,
                                             x.mu.d    = des$x.mu.d,
                                             x.mu.y    = des$x.mu.y,
                                             x.pi      = des$x.pi,
                                             z = data$dat$Z, 
                                             d = data$dat$D, 
                                             y = data$dat$Y,
                                             link_delta = link_delta
                      )
                    } else if(z.type == "continuous"){
                      data <- data_gen(alpha.true, beta.true, iota.true, theta.true, zeta.true,
                                        n, sigma=1, Y.type="binary", seed=i)
                      
                      ATE.crude <- with(data$dat, mean(Y[D==1]) - mean(Y[D==0]))
                      
                      des <- pick_scenario(data$dat, scn)
                      
                      est <- estimation(
                        x.delta   = des$x.delta,
                        x.delta.z = des$x.delta.z,
                        x.mu.d    = des$x.mu.d,
                        x.mu.y    = des$x.mu.y,
                        x.pi      = des$x.pi,
                        z = data$dat$Z, d = data$dat$D, y = data$dat$Y,
                        link_delta = link_delta
                      )
                    } else{
                      data <- data_gen(alpha.true, beta.true, iota.true, theta.true, zeta.true,
                                       n, sigma=1, Y.type="binary", seed=i)
                      
                      dicZ <- ifelse(data$dat$Z>quantile(data$dat$Z, quantile.z), 1, 0)
                      
                      ATE.crude <- with(data$dat, mean(Y[D==1]) - mean(Y[D==0]))
                      
                      des <- pick_scenario(data$dat, scn)
                      
                      est <- estimation_binZ(
                        x.delta   = des$x.delta,
                        x.delta.z = des$x.delta.z,
                        x.mu.d    = des$x.mu.d,
                        x.mu.y    = des$x.mu.y,
                        x.pi      = des$x.pi,
                        z = dicZ, d = data$dat$D, y = data$dat$Y,
                        link_delta = link_delta
                      )
                    }
                     c(unlist(est)[c("Delta1","Delta2","Delta2.b","Delta3","Delta.mr")], ATE.crude)
                   }
    
    names <- c("Delta1","Delta2", "Delta2.b", "Delta3", "Delta.mr", "ATE.crude")
    if (ncol(res) == length(names)) colnames(res) <- names
    
    list(
      raw   = res,
      mean  = apply(res, 2, mean),
      mcse  = apply(res, 2, sd) / sqrt(reps)
    )
  }
  
  # allow vector of scenarios
  out <- setNames(lapply(scenario, run_one_scenario), scenario)
  if (length(out) == 1) out[[1]] else out
}


