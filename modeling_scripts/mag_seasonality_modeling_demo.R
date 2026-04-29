# =============================================================================
# DEMO: Pseudomonas MAG Abundance Modeling in Wastewater
# Models 2, 3, and 4 from the manuscript
# NOTE: All input data here are SYNTHETIC, generated to match the structure
#       and seasonal properties of the real dataset.
# =============================================================================

library(tidyr)
library(dplyr)
library(data.table)
library(plyr)
library(ggplot2)

set.seed(42)

# =============================================================================
# SYNTHETIC DATA GENERATION
# =============================================================================
# The real dataset spans ~1 year of twice-weekly sampling.
# We simulate ~52 weeks = ~104 sampling points, aggregated to monthly x values
# (x = month index, 0.25 increments reflecting the week/4 - 0.25 transform).

n_weeks <- 52
weeks   <- 1:n_weeks

# -- Time axis (mirrors: x = ordered_week/4 - 0.25) -------------------------
x_vals <- (weeks / 4) - 0.25   # ranges ~0 to ~13 (months)

# -- Dates (twice-weekly, simplified to one date per week for demo) ----------
start_date <- as.Date("2021-01-04")
dates      <- start_date + (weeks - 1) * 7

# -- Sample IDs --------------------------------------------------------------
sample_ids <- paste0("sample_", sprintf("%03d", weeks))

# -- Influent flow (MGD): slight seasonal signal + noise --------------------
IF_vals <- 3.5 + 0.4 * cos(2 * pi / 12 * x_vals) + rnorm(n_weeks, 0, 0.15)

# -- crAssphage (min-max scaled): proxy for human fecal load -----------------
crass_vals <- 0.5 + 0.2 * cos(2 * pi / 12 * (x_vals - 2)) + rnorm(n_weeks, 0, 0.08)
crass_vals <- pmin(pmax(crass_vals, 0), 1)   # clamp to [0,1]

# -- Day of week binary (0 = Monday, 1 = Friday) ----------------------------
DW_vals <- rep(c(0L, 1L), length.out = n_weeks)

# -- Antibiotic use (SXT / folate antagonists): strong winter peak -----------
# Real SXT use peaks in winter months; we encode that here.
sxt_use_vals <- 12 + 5 * cos(2 * pi / 12 * (x_vals - 1)) + rnorm(n_weeks, 0, 1.2)

# Lag 4 weeks (k=4 in the paper) and re-align
k        <- 4
lag_idx  <- pmax(1, weeks - k)
sxt_lagged <- sxt_use_vals[lag_idx]

# Assembled antibiotic use data (merged by x)
use_data <- data.frame(x = x_vals, use = sxt_lagged)

# -- Influent flow data (merged by Date) ------------------------------------
inf_flow <- data.frame(Date = dates, IF = IF_vals)

# -- crAssphage data (merged by sample) -------------------------------------
crass_df <- data.frame(sample = sample_ids, minmax_crassphage = crass_vals)

# -- MAG abundance data: simulate 3 Pseudomonas MAGs ------------------------
# Each MAG has a distinct seasonal amplitude, phase, and response to SXT use.
mag_params <- list(
  gs43_SIM = list(amp = 0.8,  phase = 2,   B_SXT = 0.04,  B_IF = -0.05, B_Crass = 0.3,  DOW = 0.05, intercept = 2.0),
  gs77_SIM = list(amp = 0.5,  phase = 5,   B_SXT = 0.02,  B_IF =  0.03, B_Crass = 0.15, DOW = 0.02, intercept = 1.5),
  gs12_SIM = list(amp = 0.35, phase = 8.5, B_SXT = -0.01, B_IF =  0.02, B_Crass = 0.10, DOW = 0.01, intercept = 1.2)
)

omega <- 2 * pi / 12

sim_mag <- function(gnid, params) {
  p   <- params[[gnid]]
  y   <- p$amp * cos(omega * (x_vals - p$phase)) +
         p$B_SXT   * sxt_lagged     +
         p$B_IF    * IF_vals        +
         p$B_Crass * crass_vals     +
         p$DOW     * DW_vals        +
         p$intercept                +
         rnorm(n_weeks, 0, 0.12)
  y   <- pmax(y, 0)   # RPKM cannot be negative
  data.frame(
    genomeid = gnid,
    Date     = dates,
    sample   = sample_ids,
    x        = x_vals,
    y        = y,
    DW       = DW_vals
  )
}

data <- rbind(
  sim_mag("gs43_SIM", mag_params),
  sim_mag("gs77_SIM", mag_params),
  sim_mag("gs12_SIM", mag_params)
)

cat("Synthetic data generated.\n")
cat("  MAGs:    ", paste(unique(data$genomeid), collapse = ", "), "\n")
cat("  Weeks:   ", n_weeks, "\n")
cat("  Rows:    ", nrow(data), "\n\n")

# =============================================================================
# HELPER: phase conversion (keeps phases in [0, period])
# =============================================================================
convert_a_phases_func <- function(df) {
  period <- unique(df$period)[1]
  df_adj <- df
  phase_row <- which(df_adj$term == "phase")
  if (length(phase_row) > 0) {
    df_adj$estimate[phase_row] <- df_adj$estimate[phase_row] %% period
  }
  df_adj
}

# =============================================================================
# MODEL 2: Seasonal model only (no antibiotic use)
# =============================================================================
# MAG_ti ~ A * cos(omega*(ti - P)) + C
# =============================================================================

model_pseudomonas_MAGs_no_use <- function(gnid) {

  period <- 12
  omega  <- 2 * pi / period

  idf            <- data %>% subset(genomeid == gnid)
  start_intercept <- mean(idf$y)

  model <- nls(
    y ~ amplitude * cos(omega * (x - phase)) + intercept,
    start = list(amplitude = 0.1, phase = 0, intercept = start_intercept),
    data  = idf
  )

  coeffs <- summary(model)$coefficients %>%
    as_tibble(rownames = "term") %>%
    magrittr::set_colnames(c("term", "estimate", "std.error", "statistic", "p.value"))

  cints <- confint.default(model) %>%
    magrittr::set_colnames(c("ci.lower", "ci.upper")) %>%
    as_tibble(rownames = "term")

  aic <- AIC(model)

  out     <- data.frame(coeffs, cints)[, c("term","estimate","std.error","statistic","p.value","ci.lower","ci.upper")]
  out$AIC <- aic; out$period <- period
  out.adj <- convert_a_phases_func(out)
  out$p.adj <- p.adjust(out$p.value); out$raw_or_adjusted <- "raw"
  out.adj$raw_or_adjusted <- "adj"
  fo          <- rbind.fill(out, out.adj)
  fo$genomeid <- gnid; fo$period <- period

  SSE <- sum((idf$y - fitted(model))^2)
  SST <- sum((idf$y - mean(idf$y))^2)
  fo$R2 <- 1 - SSE / SST

  # store fitted values for plotting
  attr(fo, "fitted") <- data.frame(x = idf$x, y_obs = idf$y, y_fit = fitted(model), genomeid = gnid)
  return(fo)
}

mags   <- unique(data$genomeid)
wo_use <- ldply(lapply(mags, model_pseudomonas_MAGs_no_use))

cat("Model 2 (seasonal only) complete.\n")
print(wo_use %>% subset(raw_or_adjusted == "raw") %>%
        select(genomeid, term, estimate, p.value, R2) %>%
        arrange(genomeid, term))
cat("\n")

# =============================================================================
# MODEL 3: Seasonal model + antibiotic use (SXT, lagged k=4 weeks)
# =============================================================================
# MAG_ti ~ A * cos(omega*(ti - P)) + B_SXT * SXT_(ti-k) + C
# =============================================================================

model_pseudomonas_MAGs_w_use <- function(gnid) {

  period <- 12
  omega  <- 2 * pi / period

  idf             <- data %>% subset(genomeid == gnid) %>% merge(., use_data, by = "x")
  start_intercept <- mean(idf$y)

  model <- nls(
    y ~ amplitude * cos(omega * (x - phase)) + B_SXT * use + intercept,
    start = list(amplitude = 0.1, phase = 0, B_SXT = 0, intercept = start_intercept),
    data  = idf
  )

  coeffs <- summary(model)$coefficients %>%
    as_tibble(rownames = "term") %>%
    magrittr::set_colnames(c("term", "estimate", "std.error", "statistic", "p.value"))

  cints <- confint.default(model) %>%
    magrittr::set_colnames(c("ci.lower", "ci.upper")) %>%
    as_tibble(rownames = "term")

  aic <- AIC(model)

  out     <- data.frame(coeffs, cints)[, c("term","estimate","std.error","statistic","p.value","ci.lower","ci.upper")]
  out$AIC <- aic; out$period <- period
  out.adj <- convert_a_phases_func(out)
  out$p.adj <- p.adjust(out$p.value); out$raw_or_adjusted <- "raw"
  out.adj$raw_or_adjusted <- "adj"
  fo          <- rbind.fill(out, out.adj)
  fo$genomeid <- gnid; fo$period <- period

  SSE <- sum((idf$y - fitted(model))^2)
  SST <- sum((idf$y - mean(idf$y))^2)
  fo$R2 <- 1 - SSE / SST

  attr(fo, "fitted") <- data.frame(x = idf$x, y_obs = idf$y, y_fit = fitted(model), genomeid = gnid)
  return(fo)
}

w_use <- ldply(lapply(mags, model_pseudomonas_MAGs_w_use))

cat("Model 3 (seasonal + SXT use) complete.\n")
print(w_use %>% subset(raw_or_adjusted == "raw") %>%
        select(genomeid, term, estimate, p.value, R2) %>%
        arrange(genomeid, term))
cat("\n")

# =============================================================================
# MODEL 4: Full model — seasonal + SXT use + flow + crAssphage + day of week
# =============================================================================
# MAG_ti ~ A*cos(omega*(ti-P)) + B_SXT*SXT + B_IF*flow + B_Crass*crAss + DOW*wd + C
# =============================================================================

model_pseudomonas_MAGs_w_use_and_controls <- function(gnid) {

  period <- 12
  omega  <- 2 * pi / period

  idf <- data %>%
    subset(genomeid == gnid) %>%
    merge(., use_data,                        by = "x")      %>%
    merge(., inf_flow[, c("Date","IF")] %>% unique(), by = "Date") %>%
    merge(., crass_df,                        by = "sample")

  start_intercept <- mean(idf$y)

  model <- nls(
    y ~ amplitude * cos(omega * (x - phase)) +
        B_SXT   * use               +
        B_IF    * IF                +
        B_Crass * minmax_crassphage +
        DOW     * DW                +
        intercept,
    start = list(
      amplitude = 0.1, phase = 0,
      B_SXT = 0, B_Crass = 0, B_IF = 0, DOW = 0,
      intercept = start_intercept
    ),
    data = idf
  )

  coeffs <- summary(model)$coefficients %>%
    as_tibble(rownames = "term") %>%
    magrittr::set_colnames(c("term", "estimate", "std.error", "statistic", "p.value"))

  cints <- confint.default(model) %>%
    magrittr::set_colnames(c("ci.lower", "ci.upper")) %>%
    as_tibble(rownames = "term")

  aic <- AIC(model)

  out     <- data.frame(coeffs, cints)[, c("term","estimate","std.error","statistic","p.value","ci.lower","ci.upper")]
  out$AIC <- aic; out$period <- period
  out.adj <- convert_a_phases_func(out)
  out$p.adj <- p.adjust(out$p.value); out$raw_or_adjusted <- "raw"
  out.adj$raw_or_adjusted <- "adj"
  fo          <- rbind.fill(out, out.adj)
  fo$genomeid <- gnid; fo$period <- period

  SSE <- sum((idf$y - fitted(model))^2)
  SST <- sum((idf$y - mean(idf$y))^2)
  fo$R2 <- 1 - SSE / SST

  attr(fo, "fitted") <- data.frame(x = idf$x, y_obs = idf$y, y_fit = fitted(model), genomeid = gnid)
  return(fo)
}

w_use_etal <- ldply(lapply(mags, model_pseudomonas_MAGs_w_use_and_controls))

cat("Model 4 (full model) complete.\n")
print(w_use_etal %>% subset(raw_or_adjusted == "raw") %>%
        select(genomeid, term, estimate, p.value, R2) %>%
        arrange(genomeid, term))
cat("\n")

# =============================================================================
# MODEL COMPARISON: AIC across Models 2, 3, 4 per MAG
# =============================================================================

aic_compare <- rbind(
  wo_use   %>% subset(term == "intercept" & raw_or_adjusted == "raw") %>%
    select(genomeid, AIC, R2) %>% mutate(model = "Model 2: Seasonal only"),
  w_use    %>% subset(term == "intercept" & raw_or_adjusted == "raw") %>%
    select(genomeid, AIC, R2) %>% mutate(model = "Model 3: + SXT use"),
  w_use_etal %>% subset(term == "intercept" & raw_or_adjusted == "raw") %>%
    select(genomeid, AIC, R2) %>% mutate(model = "Model 4: + SXT + covariates")
)

cat("AIC comparison across models:\n")
print(aic_compare %>% arrange(genomeid, model))
cat("\n")

# =============================================================================
# OPTIONAL: Write outputs (mirrors the original write.csv calls)
# =============================================================================
# wo_use     %>% write.csv("demo-no_use-models.csv",        row.names = FALSE)
# w_use      %>% write.csv("demo-with_use-only.csv",        row.names = FALSE)
# w_use_etal %>% write.csv("demo-with_use_models_etal.csv", row.names = FALSE)

# =============================================================================
# ADDITIONAL LIBRARIES for plotting
# =============================================================================
library(ggpubr)
library(ggthemes)

# =============================================================================
# PLOTTING SETUP
# =============================================================================

# -- cleanup: minimal ggplot theme matching original style -------------------
cleanup <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(linewidth = 0.3, color = "grey90"),
    strip.background  = element_rect(fill = "#E6F1FB"),
    strip.text        = element_text(face = "bold"),
    axis.text.x       = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.y      = element_blank()
  )

# -- model_data: full merged dataset needed by getAdjustedXcorrs ------------
# (mirrors the merge pipeline in the original plotting code)
model_data <- data %>%
  merge(., use_data,                             by = "x")    %>%
  merge(., inf_flow[, c("Date", "IF")] %>% unique(), by = "Date") %>%
  merge(., crass_df,                             by = "sample")

# -- model3: alias for w_use (Model 3 results) --------------------------------
# The plotting function references 'model3' by name.
model3 <- w_use

# -- m3sig: all MAGs from Model 3, regardless of significance ----------------
# In the real analysis this would filter to MAGs where B_SXT p < 0.05.
# For the demo we include all simulated MAGs so all panels are visible.
# To apply the significance filter on real data, uncomment the line below:
# m3sig <- w_use %>% subset(term == "B_SXT" & raw_or_adjusted == "raw" & p.value < 0.05)
m3sig <- w_use %>% subset(term == "B_SXT" & raw_or_adjusted == "raw")

# =============================================================================
# PLOTTING: getAdjustedXcorrs
# Seasonal deviate plot for each significant MAG from Model 3.
# Shows: observed deviates (open points), model-predicted deviates (filled),
#        cosine fit (line), and amplitude confidence ribbon.
# =============================================================================

cos_func <- function(month, amplitude, phase, omega, intercept) {
  amplitude * cos(omega * (month - phase)) + intercept
}

getAdjustedXcorrs <- function(gnid) {

  period <- 12
  omega  <- 2 * pi / period

  working_model <- model3 %>% subset(genomeid == gnid)
  working_data  <- model_data %>% subset(genomeid == gnid)

  # Extract model parameters
  phase     <- (working_model %>% subset(raw_or_adjusted == "adj" & term == "phase"))$estimate
  amplitude <- (working_model %>% subset(raw_or_adjusted == "adj" & term == "amplitude"))$estimate
  alow      <- (working_model %>% subset(raw_or_adjusted == "adj" & term == "amplitude"))$ci.lower
  ahigh     <- (working_model %>% subset(raw_or_adjusted == "adj" & term == "amplitude"))$ci.upper
  intercept <- (working_model %>% subset(raw_or_adjusted == "raw" & term == "intercept"))$estimate
  BSXT      <- (working_model %>% subset(raw_or_adjusted == "raw" & term == "B_SXT"))$estimate

  # Compute model predictions and component contributions
  working_data$model_predicted_data        <- amplitude * cos(omega * (working_data$x - phase)) +
                                              (BSXT * working_data$use) + intercept
  working_data$model_predicted_data_no_use <- amplitude * cos(omega * (working_data$x - phase)) + intercept
  working_data$seasonality_val             <- amplitude * cos(omega * (working_data$x - phase)) + intercept
  working_data$bsxt_contribution           <- (BSXT * working_data$use) + intercept

  # Confidence ribbon from amplitude CI (fine-grained month grid)
  ci <- data.frame(month = seq(0, 12, 0.01)) %>%
    mutate(
      lower_ci = purrr::map_dbl(month, ~ cos_func(., alow,  phase, omega, 0)),
      upper_ci = purrr::map_dbl(month, ~ cos_func(., ahigh, phase, omega, 0))
    )

  res <- working_data %>%
    ggplot(aes(x = x)) +
    # Observed deviates (mean-centred)
    geom_point(aes(x = x, y = y - mean(y)),
               shape = 21, size = 2) +
    # Model-predicted deviates (filled, mean-centred)
    geom_point(aes(x = x, y = model_predicted_data - mean(model_predicted_data)),
               shape = 21, fill = "black", size = 2) +
    # Cosine fit line
    stat_function(fun  = cos_func,
                  args = list(amplitude = amplitude, phase = phase,
                              omega = omega, intercept = 0),
                  linewidth = 0.7) +
    # Amplitude CI ribbon
    geom_ribbon(data = ci,
                aes(x = month, ymin = lower_ci, ymax = upper_ci),
                alpha = 0.3) +
    cleanup +
    labs(
      y     = "Seasonal deviate",
      x     = "Month",
      title = gnid
    ) +
    scale_x_continuous(
      breaks = c(1, 3, 5, 7, 9, 11),
      labels = c("09-20", "11-20", "01-21", "03-21", "05-21", "07-21")
    )

  return(res)
}

# Generate one panel per MAG and arrange into a multi-panel figure.
# ncol=2 fits panels side by side; adjust ncol=1 for a single-column layout.
n_plots        <- length(unique(m3sig$genomeid))
list_of_plots  <- lapply(unique(m3sig$genomeid), getAdjustedXcorrs)
gp <- ggarrange(
  plotlist = list_of_plots,
  ncol     = min(2, n_plots),
  nrow     = ceiling(n_plots / min(2, n_plots))
)

print(gp)
# ggsave("demo-seasonal-deviates.pdf", gp, width = 10, height = 4 * ceiling(n_plots / 2))
