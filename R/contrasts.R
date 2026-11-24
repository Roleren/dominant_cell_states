# Contrast units: matched case/control sample sets.
#
# A contrast unit is one case condition inside one stratum (a set of metadata
# fields such as study, author, cell line, tissue and gene that must match),
# compared with the control samples of the same stratum. The pipeline scripts
# (next model, grouped branch usage, gene perturbation screen, control-aware
# enrichment) each rebuilt these pieces; this file holds the shared,
# tested definitions. Regulation v2's group registry and placebo contrasts
# build on the same units.
#
# Condition vocabulary (moved unchanged from scripts/dominant_condition_families.R,
# which now sources this file and keeps the dominant_* names).

# Base R >= 4.4 provides %||%; define it for older R.
if (!exists("%||%", envir = baseenv(), inherits = FALSE)) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}

#' Clean a metadata level
#'
#' Trims whitespace and maps empty or "NA"/"NULL"/"None" spellings to `NA`.
#' @param x Character vector.
#' @export
dcs_clean_level <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "na", "n/a", "NULL", "null",
             "None", "none")] <- NA_character_
  x
}

#' Normalise a label to lower-case snake case
#' @param x Character vector.
#' @export
dcs_normalise_label <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

#' Is a condition label a control?
#'
#' Wild type, untreated, vehicle, mock, scrambled, non-targeting and empty
#' vector spellings.
#' @param x Condition labels.
#' @export
dcs_is_control_condition <- function(x) {
  y <- dcs_normalise_label(x)
  y %in% c(
    "wt", "wildtype", "wild_type", "control", "ctrl", "mock", "vehicle",
    "dmso", "untreated", "unperturbed", "baseline", "none", "empty_vector",
    "ev", "scramble", "scrambled", "scr", "shctrl", "sh_control",
    "non_targeting", "nontargeting", "nt", "ntc", "vector"
  ) |
    grepl("^(wt|ctrl|control|mock|vehicle|dmso)$", y) |
    grepl("(^|_)(wild_type|wildtype|untreated|non_targeting|nontargeting)($|_)", y)
}

.dcs_design_family_rules <- list(
  viral_infection = function(text) {
    grepl("infect|virus|viral|influenza|adenovirus|reovirus", text) |
      grepl("(^|_)(sars|covid|cov_2|sars_cov|ebv|hsv|hiv|hcv|iav|rsv|zika|denv)($|_)", text)
  },
  acute_interferon_cytokine = function(text) {
    grepl("ifn|interferon|ifng|ifna|ifnb|cytokine|tnf|il1|il6|poly_i_c|poly_ic|lps|tlr", text)
  },
  isr_er_translation_stress = function(text) {
    grepl("arsen|tunicamycin|thapsigargin|dtt|er_stress|upr|azc|puro|harringtonine|lactimidomycin|cycloheximide|translation_stress|proteasome|bortez", text)
  },
  nutrient_starvation = function(text) {
    grepl("starv|depriv|no_b?caa|nobcaa|no_cys|nocys|no_glu|noglu|no_gln|nogln|noq|nog|glucose|amino_acid|nutrient|serum_free", text)
  },
  hypoxia_oxidative_stress = function(text) {
    grepl("hypox|hif|oxygen|anoxia|ros|oxid|h2o2|menadione", text)
  },
  heat_osmotic_recovery = function(text) {
    grepl("heat|\\bhs\\b|hs_recovery|recovery|wash|nacl|osmotic|salt", text)
  },
  dna_damage_senescence = function(text) {
    grepl("senesc|irradiat|radiation|doxo|etoposide|dna_damage|p53|uv", text)
  }
)

#' Classify a sample's design family
#'
#' First matching rule wins, in this order: control, viral infection,
#' interferon/cytokine, ISR/ER/translation stress, nutrient starvation,
#' hypoxia/oxidative stress, heat/osmotic/recovery, DNA damage/senescence,
#' tumour/cancer context, loss of function, gain of function, genetic variant,
#' drug or stimulus, immune/infection model context, fraction or protocol,
#' other. Text is the normalised concatenation of all arguments.
#' @param condition,fraction,inhibitor,cell_line,tissue,gene,cancer_type
#'   Metadata columns (vectors of equal length or length one).
#' @export
dcs_classify_design_family <- function(condition = NA_character_,
                                       fraction = NA_character_,
                                       inhibitor = NA_character_,
                                       cell_line = NA_character_,
                                       tissue = NA_character_,
                                       gene = NA_character_,
                                       cancer_type = NA_character_) {
  text <- paste(condition, fraction, inhibitor, cell_line, tissue, gene,
                cancer_type, sep = " ")
  text <- dcs_normalise_label(text)
  cancer_clean <- dcs_clean_level(cancer_type)
  cancer_present <- !is.na(cancer_clean) & nzchar(cancer_clean)
  n <- length(text)
  out <- rep(NA_character_, n)
  assign_rule <- function(hit, label) {
    hit <- rep_len(hit, n)
    take <- is.na(out) & hit
    out[take] <<- label
  }
  assign_rule(dcs_is_control_condition(condition), "control")
  for (label in names(.dcs_design_family_rules)) {
    assign_rule(.dcs_design_family_rules[[label]](text), label)
  }
  assign_rule(grepl("tumou?r|cancer|malignan|metasta|carcinoma|adenoma|normal_adjacent|primary_tumou?r", text) |
                cancer_present, "tumor_cancer_context")
  assign_rule(grepl("(^|_)(ko|knockout|knock_out|null|delet|deplet|kd|knockdown|knock_down|sirna|shrna|crispr)($|_)", text),
              "loss_of_function")
  assign_rule(grepl("(^|_)(oe|overexpression|over_expression|rescue|reconstitution|addback|complement)($|_)", text),
              "gain_of_function")
  assign_rule(grepl("mutant|mutation|variant|mut\\b|snp|allele", text),
              "genetic_variant_or_mutant")
  assign_rule(grepl("inhib|torin|rapa|rapamycin|drug|treated|treat|stim|compound|ver|pes|ga|mg132|chloroquine|actd", text),
              "drug_or_stimulus")
  assign_rule(grepl("lcl|lymphoblast|macrophage|monocyte|cd14|pbmc|blood|b_cell|t_cell|immune|hbec|a549|huh7", text),
              "immune_or_infection_model_context")
  assign_rule(grepl("fraction|polysome|monosome|ribosome|input|lysate|cyto|nuc|membrane", text),
              "fraction_or_protocol")
  out[is.na(out)] <- "other_non_control"
  out
}

#' Perturbation class of a design family
#' @param condition Condition labels (used when `design_family` is NULL).
#' @param design_family Output of [dcs_classify_design_family()].
#' @export
dcs_classify_perturbation_class <- function(condition = NA_character_,
                                            design_family = NULL) {
  if (is.null(design_family)) {
    design_family <- dcs_classify_design_family(condition)
  }
  out <- rep("treatment_or_stimulus", length(design_family))
  out[design_family == "other_non_control"] <- "other_non_control"
  out[design_family %in% c("tumor_cancer_context",
                           "immune_or_infection_model_context",
                           "fraction_or_protocol")] <- "context_or_protocol"
  out[design_family == "genetic_variant_or_mutant"] <- "genetic_perturbation"
  out[design_family == "gain_of_function"] <- "gain_of_function"
  out[design_family == "loss_of_function"] <- "loss_of_function"
  out[design_family == "control"] <- "control"
  out[is.na(design_family)] <- NA_character_
  out
}

#' Stratum identifier
#'
#' Pastes the stratum fields with " | ". Fields absent from `dt` are skipped.
#' @param dt data.frame / data.table.
#' @param fields Metadata columns that define a stratum.
#' @param missing Optional replacement for `NA` values (e.g. "MISSING").
#' @export
dcs_stratum_id <- function(dt, fields, missing = NULL) {
  fields <- intersect(fields, names(dt))
  cols <- lapply(fields, function(f) {
    v <- as.character(dt[[f]])
    if (!is.null(missing)) v[is.na(v)] <- missing
    v
  })
  do.call(paste, c(cols, sep = " | "))
}

#' Two-tier matched controls
#'
#' Within each group (e.g. study and cell line) a procedure-matched "tight"
#' control (scrambled, non-targeting, mock, empty vector) is preferred; a
#' generic "loose" baseline (WT, control, parental) counts only when the
#' group has no tight control.
#' @param condition Condition labels.
#' @param group List of grouping vectors (or one vector).
#' @param tight_pattern,loose_pattern Regular expressions on the lower-cased,
#'   trimmed condition.
#' @return Data frame with logical columns `is_tight_control`,
#'   `is_loose_control`, `group_has_tight` and `is_control`.
#' @export
dcs_control_tiers <- function(condition, group,
                              tight_pattern = "^(scram|scr|scrambled|ntc|non.?targeting|mock|ev|empty.?vector)$",
                              loose_pattern = "^(wt|control|ctrl|parental)$") {
  cond <- tolower(trimws(condition))
  tight <- grepl(tight_pattern, cond)
  loose <- grepl(loose_pattern, cond)
  if (!is.list(group)) group <- list(group)
  key <- do.call(paste, c(lapply(group, as.character), sep = "\r"))
  has_tight <- ave(tight, key, FUN = any)
  data.frame(is_tight_control = tight, is_loose_control = loose,
             group_has_tight = has_tight,
             is_control = tight | (loose & !has_tight))
}

#' Contrast units of a sample table
#'
#' One row per (stratum, case condition) that has at least one control sample
#' in the same stratum.
#' @param dt Sample-level table.
#' @param stratum_col Column holding the stratum id.
#' @param condition_col Condition column.
#' @param control_col Logical control column (`TRUE` = control; `NA` is
#'   neither case nor control).
#' @param sample_col Sample (run) id column.
#' @param carry Columns copied from the first case row (e.g. design family).
#' @return data.frame with stratum_id, case_condition, the carried columns,
#'   n_case_samples, n_control_samples, control_conditions (sorted, ";"),
#'   case_runs and control_runs (sorted, ";") and sample_set_signature.
#' @export
dcs_contrast_units <- function(dt, stratum_col = "stratum_id",
                               condition_col = "CONDITION",
                               control_col = "is_control_condition",
                               sample_col = "Run", carry = character()) {
  dt <- as.data.frame(dt, stringsAsFactors = FALSE)
  stratum <- as.character(dt[[stratum_col]])
  cond <- as.character(dt[[condition_col]])
  ctrl <- dt[[control_col]]
  run <- as.character(dt[[sample_col]])
  is_ctrl <- !is.na(ctrl) & ctrl
  is_case <- !is.na(ctrl) & !ctrl
  keep <- is_case & !is.na(stratum) & !is.na(cond)
  case_keys <- unique(data.frame(stratum = stratum[keep],
                                 condition = cond[keep],
                                 stringsAsFactors = FALSE))
  ctrl_strata <- unique(stratum[is_ctrl])
  case_keys <- case_keys[case_keys$stratum %in% ctrl_strata, , drop = FALSE]
  if (!nrow(case_keys)) {
    return(data.frame(stratum_id = character(), case_condition = character(),
                      n_case_samples = integer(), n_control_samples = integer(),
                      control_conditions = character(), case_runs = character(),
                      control_runs = character(),
                      sample_set_signature = character(),
                      stringsAsFactors = FALSE))
  }
  collapse <- function(x) paste(sort(unique(x[!is.na(x)])), collapse = ";")
  rows <- lapply(seq_len(nrow(case_keys)), function(i) {
    s <- case_keys$stratum[[i]]
    k <- case_keys$condition[[i]]
    case_i <- which(is_case & !is.na(stratum) & stratum == s &
                      !is.na(cond) & cond == k)
    ctrl_i <- which(is_ctrl & !is.na(stratum) & stratum == s)
    out <- data.frame(stratum_id = s, case_condition = k,
                      stringsAsFactors = FALSE)
    for (col in carry) out[[col]] <- dt[[col]][case_i[1]]
    out$n_case_samples <- length(unique(run[case_i]))
    out$n_control_samples <- length(unique(run[ctrl_i]))
    out$control_conditions <- collapse(cond[ctrl_i])
    out$case_runs <- collapse(run[case_i])
    out$control_runs <- collapse(run[ctrl_i])
    out$sample_set_signature <- paste(s, k, out$case_runs, out$control_runs,
                                      sep = " || ")
    out
  })
  do.call(rbind, rows)
}

#' Row indices of a contrast unit
#'
#' @param stratum,condition,is_control Sample-level vectors.
#' @param unit_stratum,unit_condition The unit.
#' @return Integer row indices of the unit's case and control samples, in
#'   table order.
#' @export
dcs_contrast_rows <- function(stratum, condition, is_control,
                              unit_stratum, unit_condition) {
  # Same rule as data.table's i = stratum == s & (is_control == TRUE |
  # CONDITION == k): comparisons with NA select nothing.
  in_stratum <- !is.na(stratum) & stratum == unit_stratum
  is_ctrl <- !is.na(is_control) & is_control
  is_case <- !is.na(condition) & condition == unit_condition
  which(in_stratum & (is_ctrl | is_case))
}
