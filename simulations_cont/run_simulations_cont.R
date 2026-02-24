
# Install necessary packages
if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak")
}

pak::pkg_install(c("batchtools@0.9.18", "here@1.0.2", "comets@0.2-2", 
                   "drf@1.1.0", "momentchi2@0.1.5", "ericstrobl/RCIT"))

library(batchtools)
library(here)
library(ggplot2)
library(data.table)

# Clear workspace
rm(list=ls())

# No. of iterations
k <- 2

# Sample sizes
n_vec <- c(500, 1000, 1500, 2000)
# Significance level
alpha <- 0.05

# Registry directory
reg_dir <- here("simulations_cont/registry")
# makeExperimentRegistry(file.dir = reg_dir, seed = 42)
# unlink(reg_dir, recursive = TRUE)
# Load registry
reg <- loadRegistry(reg_dir, writeable = T)
# Set cluster function (socket for Windows)
reg$cluster.functions <- makeClusterFunctionsSocket(parallel::detectCores() - 1)
reg$packages <- c(reg$packages, "RCIT")

source(here("simulations_cont/problems_cont.R"))
addProblem(name = "Null_1", fun = dgp_null1, seed = 100)
addProblem(name = "Null_2", fun = dgp_null2, seed = 200)
addProblem(name = "Alt_1",  fun = dgp_alt1,  seed = 300)
addProblem(name = "Alt_2",  fun = dgp_alt2,  seed = 400)

source(here("simulations_cont/algorithms_cont.R"))
addAlgorithm(name = "GKCM", fun = GKCM_wrapper)
addAlgorithm(name = "PCM",  fun = PCM_wrapper)
addAlgorithm(name = "GCM",  fun = GCM_wrapper)
addAlgorithm(name = "wGCM", fun = wGCM_wrapper)
addAlgorithm(name = "KCIT", fun = KCIT_wrapper)
addAlgorithm(name = "RCIT", fun = RCIT_wrapper)
addAlgorithm(name = "RCoT", fun = RCoT_wrapper)

# Define parameter grids for problems and algorithms

# Run higher sample sizes only for normal errors
prob_grid_cont <- expand.grid(n = n_vec)
prob_dsgn_cont <- list(
  Null_1  = prob_grid_cont,
  Null_2  = prob_grid_cont,
  Alt_1 = prob_grid_cont,
  Alt_2 = prob_grid_cont
)

algo_dsgn_cont <- list(
  GKCM = expand.grid(),
  PCM  = expand.grid(),
  GCM  = expand.grid(),
  wGCM = expand.grid(),
  KCIT = expand.grid(),
  RCIT = expand.grid(),
  RCoT = expand.grid()
)

# Add experiments (problem x problem parameters x 
# algorithm x hyperparameters x iterations)
addExperiments(prob_dsgn_cont, algo_dsgn_cont, repls = k)
# removeExperiments(findExperiments())

# Get an overview of jobs and test one
summarizeExperiments()
testJob(id = 216)

# Submit
submitJobs(resources = list(name = reg_dir, memory = 3000, walltime = 1800))

# Check results
getStatus()

# In case there were errors, inspect and re-run these
err_ids <- findErrors()
getErrorMessages(err_ids, missing.as.error = TRUE)
submitJobs(err_ids)          
# Rinse and repeat

# In case everything is fine save results
res <- flatten(ijoin(reduceResultsDataTable(), getJobPars()))

# Plot results 
prop_dt <- res[, .(
  r     = sum(result.1 < alpha, na.rm = TRUE),   # rejections
  m     = sum(!is.na(result.1)),                 # iterations without NA
  iters = .N                                     # num iterations (k)
), by = .(problem, n, algorithm)
][, prop_significant := r / m][]

# Wilson 95% CI for the binomial proportion
prop_dt[, c("ci_low", "ci_high") := {
  ci <- binom::binom.confint(x = r, n = m, methods = "wilson")
  list(ci$lower, ci$upper)
}, by = .(problem, n, algorithm)]

prop_dt[, problem   := factor(problem,   levels = c("Null_1", "Null_2", "Alt_1", "Alt_2"))]
prop_dt[, n         := factor(n,         levels = sort(unique(n)))]
prop_dt[, algorithm := factor(algorithm, levels = c("GCM","wGCM","PCM","KCIT","RCIT","RCoT","GKCM"))]

cols_alg <- c(
  "GCM"   = "#6F4A9F", 
  "wGCM"  = "#51127C", 
  "PCM"   = "#3A0A64",
  "KCIT"  = "#FDB8A0",
  "RCIT"  = "#FB9471", 
  "RCoT"  = "#E77D59",
  "GKCM"  = "#B63679" 
)

shapes_alg <- c(
  "GCM"   = 16,  
  "PCM"   = 16, # 18  
  "wGCM"  = 16, # 17
  "KCIT"  = 16,
  "RCIT"  = 16,
  "RCoT"  = 16,
  "GKCM"  = 16  
)

pd <- position_dodge(width = 0.5)


pbs <- unique(prop_dt$problem)
for (pb in pbs) {
  df <- prop_dt[problem == pb]
  
  p <- ggplot(
    df,
    aes(x = algorithm, y = prop_significant, group = algorithm, color = algorithm)
  ) +
    geom_hline(yintercept = 0.05, linetype = "dashed") +
    geom_errorbar(
      aes(ymin = pmax(ci_low, 0), ymax = pmin(ci_high, 1)),
      position = pd,
      linewidth = 0.6,
      width = 0.6,
      alpha = 1
    ) +
    geom_point(
      aes(shape = algorithm),
      position = pd,
      size = 2.4,
      stroke = 0.4,
      alpha = 1
    ) +
    scale_color_manual(values = cols_alg, name = "Algorithm") +
    scale_shape_manual(values = shapes_alg, name = "Algorithm") +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.02, 0.02))) +
    labs(
      x = "Algorithm",
      y = element_blank(),
      title = "Proportion significant at \u03B1 = 0.05"
    ) +
    facet_grid(
      ~ n,
      labeller = labeller(n = function(x) paste0("n = ", x))
    ) +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid.minor  = element_blank(),
      legend.position   = "none",
      panel.border      = element_rect(colour = "grey50", fill = NA, linewidth = 0.5),
      axis.text.x       = element_text(angle = 45, hjust = 1, vjust = 1)
    )
  
  
  file <- sprintf("prop_significant_%s.pdf",
                  stringr::str_replace_all(pb, "[^[:alnum:]]+", "_"))
  ggsave(filename = file, plot = p, width = 4.5, height = 2.25, device = "pdf")
}
