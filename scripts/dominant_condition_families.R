suppressPackageStartupMessages({
  library(data.table)
})

dominant_clean_level <- function(x) {
  x <- trimws(as.character(x))
  x[x %chin% c("", "NA", "N/A", "na", "n/a", "NULL", "null",
               "None", "none")] <- NA_character_
  x
}

dominant_normalise_label <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

dominant_is_control_condition <- function(x) {
  y <- dominant_normalise_label(x)
  y %chin% c(
    "wt", "wildtype", "wild_type", "control", "ctrl", "mock", "vehicle",
    "dmso", "untreated", "unperturbed", "baseline", "none", "empty_vector",
    "ev", "scramble", "scrambled", "scr", "shctrl", "sh_control",
    "non_targeting", "nontargeting", "nt", "ntc", "vector"
  ) |
    grepl("^(wt|ctrl|control|mock|vehicle|dmso)$", y) |
    grepl("(^|_)(wild_type|wildtype|untreated|non_targeting|nontargeting)($|_)", y)
}

dominant_classify_design_family <- function(condition = NA_character_,
                                            fraction = NA_character_,
                                            inhibitor = NA_character_,
                                            cell_line = NA_character_,
                                            tissue = NA_character_,
                                            gene = NA_character_,
                                            cancer_type = NA_character_) {
  text <- paste(condition, fraction, inhibitor, cell_line, tissue, gene,
                cancer_type, sep = " ")
  text <- dominant_normalise_label(text)
  cancer_clean <- dominant_clean_level(cancer_type)
  cancer_present <- !is.na(cancer_clean) & nzchar(cancer_clean)

  fifelse(
    dominant_is_control_condition(condition),
    "control",
    fifelse(
      grepl("infect|virus|viral|influenza|adenovirus|reovirus", text) |
        grepl("(^|_)(sars|covid|cov_2|sars_cov|ebv|hsv|hiv|hcv|iav|rsv|zika|denv)($|_)", text),
      "viral_infection",
      fifelse(
        grepl("ifn|interferon|ifng|ifna|ifnb|cytokine|tnf|il1|il6|poly_i_c|poly_ic|lps|tlr", text),
        "acute_interferon_cytokine",
        fifelse(
          grepl("arsen|tunicamycin|thapsigargin|dtt|er_stress|upr|azc|puro|harringtonine|lactimidomycin|cycloheximide|translation_stress|proteasome|bortez", text),
          "isr_er_translation_stress",
          fifelse(
            grepl("starv|depriv|no_b?caa|nobcaa|no_cys|nocys|no_glu|noglu|no_gln|nogln|noq|nog|glucose|amino_acid|nutrient|serum_free", text),
            "nutrient_starvation",
            fifelse(
              grepl("hypox|hif|oxygen|anoxia|ros|oxid|h2o2|menadione", text),
              "hypoxia_oxidative_stress",
              fifelse(
                grepl("heat|\\bhs\\b|hs_recovery|recovery|wash|nacl|osmotic|salt", text),
                "heat_osmotic_recovery",
                fifelse(
                  grepl("senesc|irradiat|radiation|doxo|etoposide|dna_damage|p53|uv", text),
                  "dna_damage_senescence",
                  fifelse(
                    grepl("tumou?r|cancer|malignan|metasta|carcinoma|adenoma|normal_adjacent|primary_tumou?r", text) |
                      cancer_present,
                    "tumor_cancer_context",
                    fifelse(
                      grepl("(^|_)(ko|knockout|knock_out|null|delet|deplet|kd|knockdown|knock_down|sirna|shrna|crispr)($|_)", text),
                      "loss_of_function",
                      fifelse(
                        grepl("(^|_)(oe|overexpression|over_expression|rescue|reconstitution|addback|complement)($|_)", text),
                        "gain_of_function",
                        fifelse(
                          grepl("mutant|mutation|variant|mut\\b|snp|allele", text),
                          "genetic_variant_or_mutant",
                          fifelse(
                            grepl("inhib|torin|rapa|rapamycin|drug|treated|treat|stim|compound|ver|pes|ga|mg132|chloroquine|actd", text),
                            "drug_or_stimulus",
                            fifelse(
                              grepl("lcl|lymphoblast|macrophage|monocyte|cd14|pbmc|blood|b_cell|t_cell|immune|hbec|a549|huh7", text),
                              "immune_or_infection_model_context",
                              fifelse(
                                grepl("fraction|polysome|monosome|ribosome|input|lysate|cyto|nuc|membrane", text),
                                "fraction_or_protocol",
                                "other_non_control"
                              )
                            )
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  )
}

dominant_classify_perturbation_class <- function(condition = NA_character_,
                                                 design_family = NULL) {
  if (is.null(design_family)) {
    design_family <- dominant_classify_design_family(condition)
  }
  fifelse(
    design_family == "control",
    "control",
    fifelse(
      design_family == "loss_of_function",
      "loss_of_function",
      fifelse(
        design_family == "gain_of_function",
        "gain_of_function",
        fifelse(
          design_family %chin% c("genetic_variant_or_mutant"),
          "genetic_perturbation",
          fifelse(
            design_family %chin% c("tumor_cancer_context",
                                   "immune_or_infection_model_context",
                                   "fraction_or_protocol"),
            "context_or_protocol",
            fifelse(
              design_family == "other_non_control",
              "other_non_control",
              "treatment_or_stimulus"
            )
          )
        )
      )
    )
  )
}
