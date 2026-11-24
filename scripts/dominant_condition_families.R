# Condition vocabulary for the step scripts. The definitions live in
# R/contrasts.R (dcs_* functions, tested); this file keeps the dominant_*
# names the scripts use.
suppressPackageStartupMessages({
  library(data.table)
})

source(file.path(Sys.getenv("DOMINANT_CELL_STATES_DIR", "."), "R", "contrasts.R"),
       local = TRUE)

dominant_clean_level <- dcs_clean_level
dominant_normalise_label <- dcs_normalise_label
dominant_is_control_condition <- dcs_is_control_condition
dominant_classify_design_family <- dcs_classify_design_family
dominant_classify_perturbation_class <- dcs_classify_perturbation_class
