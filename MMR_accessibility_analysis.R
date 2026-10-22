############################################################
# Spatial accessibility to critical maternal care and maternal mortality
# Models: negative binomial regression, mixed-effects NB, scenarios
# Input: integrated county-level Excel dataset
############################################################

# 1. Packages ---------------------------------------------------------------
required_packages <- c("readxl", "dplyr", "MASS", "glmmTMB", "ggplot2", "writexl")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Please install the following packages first: ", paste(missing_packages, collapse = ", "),
       "\nRun: install.packages(c(", paste(sprintf('"%s"', missing_packages), collapse = ", "), "))")
}

library(dplyr)
library(ggplot2)

# 2. Input file and output folder ------------------------------------------
# If the file is not in the current working directory, a file-selection window opens.
data_file <- "data.xlsx"  # Or replace with the full path to your Excel file.
if (!file.exists(data_file)) {
  if (interactive()) {
    message("Please select the integrated county-level Excel dataset.")
    data_file <- file.choose()
  } else {
    stop("Input file not found. Set 'data_file' to the Excel file path before running.")
  }
}
output_dir <- dirname(normalizePath(data_file, winslash = "/", mustWork = TRUE))
data <- readxl::read_excel(data_file)

# Map original Excel headers to English internal variable names.
# Unicode escapes preserve compatibility with the original Chinese column names.

for (source_name in names(source_column_map)) {
  if (source_name %in% names(data)) {
    names(data)[names(data) == source_name] <- unname(source_column_map[[source_name]])
  }
}


# 3. Validate and prepare variables ----------------------------------------
coverage_var <- "population_coverage_120min_pct"
required_columns <- c(
  "maternal_deaths_2022_2024", "live_births_2022_2024", "maternal_deaths_2024", "live_births_2024",
  coverage_var, "GDP_10k_RMB", "doctor_per_10000_women",
  "bed_per_10000_women", "province_n.x"
)
missing_columns <- setdiff(required_columns, names(data))
if (length(missing_columns) > 0) {
  stop("The Excel file is missing required columns: ", paste(missing_columns, collapse = ", "))
}

to_numeric <- function(x) {
  suppressWarnings(as.numeric(gsub("[,%[:space:]]", "", trimws(as.character(x)))))
}
safe_log <- function(x) {
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x) & x > 0
  out[ok] <- log(x[ok])
  out
}
safe_log1p <- function(x) {
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x) & x >= 0
  out[ok] <- log1p(x[ok])
  out
}

numeric_columns <- c(
  "maternal_deaths_2022_2024", "live_births_2022_2024", "maternal_deaths_2024", "live_births_2024",
  coverage_var, "GDP_10k_RMB", "doctor_per_10000_women", "bed_per_10000_women"
)
data <- data %>%
  mutate(across(all_of(numeric_columns), to_numeric)) %>%
  mutate(
    log_GDP = safe_log(`GDP_10k_RMB`),
    log_doctor = safe_log1p(doctor_per_10000_women),
    log_bed = safe_log1p(bed_per_10000_women),
    province_n.x = factor(province_n.x),
    access_group5 = case_when(
      .data[[coverage_var]] < 85 ~ "<85%",
      .data[[coverage_var]] < 90 ~ "85-<90%",
      .data[[coverage_var]] < 95 ~ "90-<95%",
      .data[[coverage_var]] < 100 ~ "95-<100%",
      .data[[coverage_var]] == 100 ~ "100%",
      TRUE ~ NA_character_
    )
  )

access_levels <- c("100%", "95-<100%", "90-<95%", "85-<90%", "<85%")
data$access_group5 <- factor(data$access_group5, levels = access_levels)

if (any(data[[coverage_var]] < 0 | data[[coverage_var]] > 100, na.rm = TRUE)) {
  stop("Accessibility coverage values must be on a 0–100 scale. Please check the coverage column.")
}

check_access_groups <- function(dat, label, require_all = TRUE) {
  observed <- unique(as.character(stats::na.omit(dat$access_group5)))
  missing <- setdiff(access_levels, observed)
  if (!"100%" %in% observed) stop(label, ": the reference group (100%) is absent.")
  if (length(observed) < 2) stop(label, ": fewer than two accessibility groups are present.")
  if (require_all && length(missing) > 0) {
    stop(label, ": these accessibility groups are absent: ", paste(missing, collapse = ", "))
  }
}

# Common complete-case sample for the three-year main analysis.
data_main <- data %>%
  filter(
    is.finite(`maternal_deaths_2022_2024`), `maternal_deaths_2022_2024` >= 0,
    is.finite(`live_births_2022_2024`), `live_births_2022_2024` > 0,
    !is.na(access_group5), is.finite(log_GDP), is.finite(log_doctor),
    !is.na(province_n.x)
  )
check_access_groups(data_main, "Main analysis")
message("Main-analysis sample size: ", nrow(data_main))
message("Counties by accessibility group:")
print(table(data_main$access_group5))

# 4. Model functions --------------------------------------------------------
make_formula <- function(outcome, births, terms, random_intercept = FALSE) {
  rhs <- c(terms, sprintf("offset(log(`%s`))", births))
  if (random_intercept) rhs <- c(rhs, "(1 | province_n.x)")
  stats::as.formula(paste0("`", outcome, "` ~ ", paste(rhs, collapse = " + ")))
}

fit_nb_models <- function(dat, outcome, births, access_type = c("categorical", "continuous"), density_var = "log_doctor") {
  access_type <- match.arg(access_type)
  access_term <- if (access_type == "continuous") paste0("`", coverage_var, "`") else "access_group5"
  f1 <- make_formula(outcome, births, access_term)
  f2 <- make_formula(outcome, births, c(access_term, "log_GDP"))
  f3 <- make_formula(outcome, births, c(access_term, "log_GDP", density_var))
  f4 <- make_formula(outcome, births, c(access_term, "log_GDP", density_var), random_intercept = TRUE)
  list(
    "Model 1" = MASS::glm.nb(f1, data = dat),
    "Model 2" = MASS::glm.nb(f2, data = dat),
    "Model 3" = MASS::glm.nb(f3, data = dat),
    "Model 4" = glmmTMB::glmmTMB(f4, family = glmmTMB::nbinom2(), data = dat)
  )
}

extract_coefficients <- function(model, model_name) {
  x <- if (inherits(model, "glmmTMB")) summary(model)$coefficients$cond else coef(summary(model))
  data.frame(
    term = rownames(x), estimate = x[, "Estimate"], se = x[, "Std. Error"],
    p = x[, "Pr(>|z|)"], Model = model_name, stringsAsFactors = FALSE
  )
}

label_term <- function(term) {
  if (term == coverage_var) return("120-min population coverage (%)")
  if (term == "log_GDP") return("GDP (log-transformed)")
  if (term == "log_doctor") return("Physician density (log-transformed)")
  if (term == "log_bed") return("Bed density (log-transformed)")
  if (startsWith(term, "access_group5")) {
    level <- substring(term, nchar("access_group5") + 1L)
    labels <- c("95-<100%" = "95–<100%", "90-<95%" = "90–<95%",
                "85-<90%" = "85–<90%", "<85%" = "<85%")
    if (level %in% names(labels)) return(unname(labels[[level]]))
  }
  term
}

format_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", paste0("=", sprintf("%.3f", p))))
}

make_results_table <- function(models, variable_order, add_reference = FALSE) {
  pieces <- lapply(names(models), function(model_name) {
    x <- extract_coefficients(models[[model_name]], model_name)
    x <- x[x$term != "(Intercept)", , drop = FALSE]
    x$Variable <- vapply(x$term, label_term, character(1))
    x$IRR <- exp(x$estimate)
    x$Lower95 <- exp(x$estimate - 1.96 * x$se)
    x$Upper95 <- exp(x$estimate + 1.96 * x$se)
    x$result <- paste0(sprintf("%.3f", x$IRR), " (", sprintf("%.3f", x$Lower95),
                        "–", sprintf("%.3f", x$Upper95), "); p", format_p(x$p))
    out <- x[, c("Variable", "result"), drop = FALSE]
    names(out)[2] <- model_name
    out
  })
  out <- Reduce(function(x, y) dplyr::full_join(x, y, by = "Variable"), pieces)

  if (add_reference) {
    ref <- data.frame(Variable = "100%", check.names = FALSE)
    for (nm in names(models)) ref[[nm]] <- "1.000 (Reference)"
    out <- dplyr::bind_rows(ref, out)
  }
  out$sort_order <- match(out$Variable, variable_order)
  out <- out[order(out$sort_order, na.last = TRUE), , drop = FALSE]
  out$sort_order <- NULL
  out <- dplyr::select(out, Variable, dplyr::all_of(names(models)))

  info <- data.frame(Variable = c("N", "Province-level random intercept"), check.names = FALSE)
  for (nm in names(models)) {
    info[[nm]] <- c(as.character(stats::nobs(models[[nm]])), ifelse(nm == "Model 4", "Yes", "No"))
  }
  dplyr::bind_rows(out, info)
}

save_xlsx <- function(x, filename) {
  writexl::write_xlsx(x, path = file.path(output_dir, filename))
}

# 5. Main analysis: continuous accessibility -------------------------------
models_continuous <- fit_nb_models(data_main, "maternal_deaths_2022_2024", "live_births_2022_2024", "continuous", "log_doctor")
table_continuous <- make_results_table(
  models_continuous,
  c("120-min population coverage (%)", "GDP (log-transformed)", "Physician density (log-transformed)")
)
save_xlsx(table_continuous, "table_continuous_accessibility_models.xlsx")

# 6. Main analysis: categorical accessibility ------------------------------
models_group <- fit_nb_models(data_main, "maternal_deaths_2022_2024", "live_births_2022_2024", "categorical", "log_doctor")
table_group <- make_results_table(
  models_group,
  c("100%", "95–<100%", "90–<95%", "85–<90%", "<85%",
    "GDP (log-transformed)", "Physician density (log-transformed)"),
  add_reference = TRUE
)
save_xlsx(table_group, "table_categorical_accessibility_models.xlsx")

# Forest plot for categorical models.
forest_df <- dplyr::bind_rows(lapply(names(models_group), function(nm) extract_coefficients(models_group[[nm]], nm))) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    Variable = vapply(term, label_term, character(1)),
    Variable = case_when(
      Variable %in% c("95–<100%", "90–<95%", "85–<90%", "<85%") ~ paste("Population coverage", Variable),
      Variable == "GDP (log-transformed)" ~ "GDP",
      Variable == "Physician density (log-transformed)" ~ "Physician density",
      TRUE ~ Variable
    ),
    IRR = exp(estimate),
    Lower95 = exp(estimate - 1.96 * se),
    Upper95 = exp(estimate + 1.96 * se),
    significance = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ ""),
    Model = factor(Model, levels = names(models_group))
  ) %>%
  filter(is.finite(IRR), is.finite(Lower95), is.finite(Upper95))

forest_order <- c("Physician density", "GDP", "Population coverage <85%", "Population coverage 85–<90%",
                  "Population coverage 90–<95%", "Population coverage 95–<100%")
forest_df$Variable <- factor(forest_df$Variable, levels = forest_order)
forest_df <- forest_df %>% filter(!is.na(Variable))

if (nrow(forest_df) > 0) {
  x_min <- min(0.5, 2^floor(log(min(forest_df$Lower95), base = 2)))
  x_max <- max(5.5, 2^ceiling(log(max(forest_df$Upper95) * 1.15, base = 2)))
  x_breaks <- 2^(seq(floor(log(x_min, base = 2)), ceiling(log(x_max, base = 2)), by = 1))
  dodge <- position_dodge(width = 0.72)
  p_forest <- ggplot(forest_df, aes(x = IRR, y = Variable, colour = Model)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45", linewidth = 0.65) +
    geom_segment(aes(x = Lower95, xend = Upper95, y = Variable, yend = Variable),
                 linewidth = 0.8, position = dodge) +
    geom_point(size = 2.7, position = dodge) +
    geom_text(aes(x = Upper95 * 1.04, label = significance), position = dodge,
              hjust = 0, size = 3.5, fontface = "bold", show.legend = FALSE) +
    scale_colour_manual(values = c("Model 1" = "#00468B", "Model 2" = "#925E9F",
                                   "Model 3" = "#42B540", "Model 4" = "#ED0000"), name = NULL) +
    scale_x_log10(breaks = x_breaks, labels = function(x) format(x, trim = TRUE, scientific = FALSE)) +
    coord_cartesian(xlim = c(x_min, x_max), clip = "off") +
    labs(x = "Incidence rate ratio (IRR)", y = NULL,
         caption = "* p<0.05; ** p<0.01; *** p<0.001") +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom", legend.direction = "horizontal",
          axis.text = element_text(colour = "black"),
          axis.title.x = element_text(face = "bold", margin = margin(t = 8)),
          plot.caption = element_text(hjust = 0.5, margin = margin(t = 6)))
  ggsave(file.path(output_dir, "forest_plot_categorical_models.png"), p_forest,
         width = 20, height = 15, units = "cm", dpi = 600, bg = "white")
}

# 7. Counterfactual scenario function --------------------------------------
run_scenario <- function(dat, model, new_group, scenario, policy, births_var) {
  original_group <- as.character(dat$access_group5)
  if (length(new_group) != nrow(dat) || anyNA(new_group) || any(!new_group %in% access_levels)) {
    stop("Scenario contains missing, invalid, or incorrectly sized accessibility groups.")
  }
  scenario_data <- dat
  scenario_data$access_group5 <- factor(new_group, levels = access_levels)
  changed <- original_group != new_group
  baseline <- predict(model, newdata = dat, type = "response", re.form = NULL, allow.new.levels = FALSE)
  counterfactual <- predict(model, newdata = scenario_data, type = "response", re.form = NULL, allow.new.levels = FALSE)
  baseline_all <- sum(baseline)
  counterfactual_all <- sum(counterfactual)
  baseline_affected <- sum(baseline[changed])
  counterfactual_affected <- sum(counterfactual[changed])
  avoided_all <- baseline_all - counterfactual_all
  avoided_affected <- baseline_affected - counterfactual_affected

  summary <- data.frame(
    Scenario = scenario, Policy = policy,
    Counties_affected = sum(changed),
    Live_births_affected = sum(dat[[births_var]][changed], na.rm = TRUE),
    Current_predicted_deaths_all = baseline_all,
    Policy_predicted_deaths_all = counterfactual_all,
    Deaths_avoided = avoided_all,
    National_reduction_percent = ifelse(baseline_all > 0, avoided_all / baseline_all * 100, NA_real_),
    Current_predicted_deaths_affected = baseline_affected,
    Policy_predicted_deaths_affected = counterfactual_affected,
    Deaths_avoided_affected = avoided_affected,
    Reduction_percent_affected = ifelse(baseline_affected > 0, avoided_affected / baseline_affected * 100, NA_real_),
    check.names = FALSE
  )
  county_data <- data.frame(dat, baseline_predicted_deaths = baseline,
                            counterfactual_predicted_deaths = counterfactual,
                            deaths_avoided = baseline - counterfactual,
                            changed = changed)
  list(summary = summary, county_data = county_data)
}

format_scenario_results <- function(results, full_benefit) {
  out <- dplyr::bind_rows(lapply(results, `[[`, "summary")) %>%
    mutate(Proportion_of_full_coverage_benefit = ifelse(full_benefit != 0, Deaths_avoided / full_benefit * 100, NA_real_))
  numeric_cols <- names(out)[vapply(out, is.numeric, logical(1))]
  out[numeric_cols] <- lapply(out[numeric_cols], round, digits = 2)
  if ("Proportion_of_full_coverage_benefit" %in% names(out)) {
    out$Proportion_of_full_coverage_benefit <- round(out$Proportion_of_full_coverage_benefit, 1)
  }
  out
}

# 8. Main counterfactual scenarios (three-year outcome) --------------------
check_access_groups(data_main, "Three-year counterfactual scenarios", require_all = TRUE)
original_group <- as.character(data_main$access_group5)
group_s1 <- ifelse(original_group == "<85%", "85-<90%", original_group)
group_s2 <- c("<85%" = "85-<90%", "85-<90%" = "90-<95%", "90-<95%" = "95-<100%",
              "95-<100%" = "100%", "100%" = "100%")[original_group]
group_s3 <- rep("100%", length(original_group))

S1 <- run_scenario(data_main, models_group[["Model 4"]], group_s1, "Scenario 1", "<85% to 85-<90%", "live_births_2022_2024")
S2 <- run_scenario(data_main, models_group[["Model 4"]], group_s2, "Scenario 2", "Stepwise improvement by one category", "live_births_2022_2024")
S3 <- run_scenario(data_main, models_group[["Model 4"]], group_s3, "Scenario 3", "All counties below 100% to 100%", "live_births_2022_2024")
scenario_results <- format_scenario_results(list(S1, S2, S3), S3$summary$Deaths_avoided)
save_xlsx(scenario_results, "counterfactual_scenarios_three_year.xlsx")

# 9. Sensitivity analysis 1: 2024 deaths and live births --------------------
data_2024 <- data %>%
  filter(
    is.finite(`maternal_deaths_2024`), `maternal_deaths_2024` >= 0,
    is.finite(`live_births_2024`), `live_births_2024` > 0,
    !is.na(access_group5), is.finite(log_GDP), is.finite(log_doctor),
    !is.na(province_n.x)
  )
check_access_groups(data_2024, "2024 sensitivity analysis", require_all = TRUE)
message("2024 sensitivity-analysis sample size: ", nrow(data_2024))
models_2024 <- fit_nb_models(data_2024, "maternal_deaths_2024", "live_births_2024", "categorical", "log_doctor")
table_2024 <- make_results_table(
  models_2024,
  c("100%", "95–<100%", "90–<95%", "85–<90%", "<85%",
    "GDP (log-transformed)", "Physician density (log-transformed)"),
  add_reference = TRUE
)
save_xlsx(table_2024, "sensitivity_2024_models.xlsx")

# 2024 counterfactual scenarios use exactly the same sample as models_2024.
original_group_2024 <- as.character(data_2024$access_group5)
group_2024_s1 <- ifelse(original_group_2024 == "<85%", "85-<90%", original_group_2024)
group_2024_s2 <- c("<85%" = "85-<90%", "85-<90%" = "90-<95%", "90-<95%" = "95-<100%",
                   "95-<100%" = "100%", "100%" = "100%")[original_group_2024]
group_2024_s3 <- rep("100%", nrow(data_2024))
S1_2024 <- run_scenario(data_2024, models_2024[["Model 4"]], group_2024_s1,
                        "Scenario 1", "<85% to 85-<90%", "live_births_2024")
S2_2024 <- run_scenario(data_2024, models_2024[["Model 4"]], group_2024_s2,
                        "Scenario 2", "Stepwise improvement by one category", "live_births_2024")
S3_2024 <- run_scenario(data_2024, models_2024[["Model 4"]], group_2024_s3,
                        "Scenario 3", "All counties below 100% to 100%", "live_births_2024")
scenario_results_2024 <- format_scenario_results(list(S1_2024, S2_2024, S3_2024), S3_2024$summary$Deaths_avoided)
save_xlsx(scenario_results_2024, "counterfactual_scenarios_2024.xlsx")

# 10. Sensitivity analysis 2: replace physician density with bed density ---
data_bed <- data_main %>% filter(is.finite(log_bed))
check_access_groups(data_bed, "Bed-density sensitivity analysis", require_all = FALSE)
message("Bed-density sensitivity-analysis sample size: ", nrow(data_bed))
models_bed <- fit_nb_models(data_bed, "maternal_deaths_2022_2024", "live_births_2022_2024", "categorical", "log_bed")
table_bed <- make_results_table(
  models_bed,
  c("100%", "95–<100%", "90–<95%", "85–<90%", "<85%",
    "GDP (log-transformed)", "Bed density (log-transformed)"),
  add_reference = TRUE
)
save_xlsx(table_bed, "sensitivity_bed_density_models.xlsx")

message("Analysis completed. Output files saved to: ", output_dir)
