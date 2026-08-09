#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

stopf <- function(...) {
  stop(sprintf(...), call. = FALSE)
}

messagef <- function(...) {
  cat(sprintf(...), "\n", sep = "")
}

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) == 0) {
    return(list())
  }

  parsed <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stopf("Unexpected argument: %s", key)
    }

    key <- sub("^--", "", key)
    if (i == length(args) || startsWith(args[[i + 1]], "--")) {
      parsed[[key]] <- TRUE
      i <- i + 1
      next
    }

    parsed[[key]] <- args[[i + 1]]
    i <- i + 2
  }

  parsed
}

require_args <- function(opts, required) {
  missing <- required[!required %in% names(opts)]
  if (length(missing) > 0) {
    stopf("Missing required arguments: %s", paste(missing, collapse = ", "))
  }
}

trim_names <- function(x) {
  gsub("^\\s+|\\s+$", "", x)
}

detect_scale_column <- function(ref) {
  if ("Scale" %in% names(ref)) {
    return("Scale")
  }
  if ("SD" %in% names(ref)) {
    return("SD")
  }
  if ("SE" %in% names(ref)) {
    messagef("Warning: profile reference uses column 'SE'. The software will treat it as a scale term.")
    return("SE")
  }
  stopf("Profile reference must contain one of: Scale, SD, or SE.")
}

load_profile_reference <- function(path) {
  ref <- fread(path, data.table = FALSE)
  names(ref) <- trim_names(names(ref))

  if (!"File" %in% names(ref)) {
    stopf("Profile reference is missing the 'File' column: %s", path)
  }
  if (!"Mean" %in% names(ref)) {
    stopf("Profile reference is missing the 'Mean' column: %s", path)
  }

  scale_col <- detect_scale_column(ref)
  keep_cols <- unique(c("File", "Mean", scale_col, "Max", "Min"))
  ref <- ref[, keep_cols[keep_cols %in% names(ref)], drop = FALSE]

  names(ref)[names(ref) == "File"] <- "feature_id"
  names(ref)[names(ref) == "Mean"] <- "ref_mean"
  names(ref)[names(ref) == scale_col] <- "ref_scale"
  names(ref)[names(ref) == "Max"] <- "ref_max"
  names(ref)[names(ref) == "Min"] <- "ref_min"

  ref$feature_id <- as.character(ref$feature_id)
  ref$ref_mean <- as.numeric(ref$ref_mean)
  ref$ref_scale <- as.numeric(ref$ref_scale)
  if ("ref_max" %in% names(ref)) {
    ref$ref_max <- as.numeric(ref$ref_max)
  }
  if ("ref_min" %in% names(ref)) {
    ref$ref_min <- as.numeric(ref$ref_min)
  }

  ref
}

standardize_feature <- function(values, mean_value, scale_value) {
  if (is.na(scale_value) || scale_value <= 0) {
    return(rep(0, length(values)))
  }
  (values - mean_value) / scale_value
}

read_score_file <- function(path) {
  dat <- fread(path, data.table = FALSE)
  name_map <- setNames(names(dat), names(dat))

  if ("#IID" %in% names(dat) && !"IID" %in% names(dat)) {
    name_map["#IID"] <- "IID"
  }
  names(dat) <- unname(name_map)

  if (!"IID" %in% names(dat)) {
    stopf("Score file is missing IID column: %s", path)
  }

  score_candidates <- c("SCORESUM", "SCORE1_SUM", "SCORE")
  score_col <- score_candidates[score_candidates %in% names(dat)]
  if (length(score_col) == 0) {
    stopf("Could not find score column in %s. Expected one of: %s",
          path, paste(score_candidates, collapse = ", "))
  }

  data.frame(
    IID = as.character(dat$IID),
    raw_score = as.numeric(dat[[score_col[[1]]]]),
    stringsAsFactors = FALSE
  )
}

merge_score_tables <- function(score_tables) {
  merged <- score_tables[[1]]
  if (length(score_tables) == 1) {
    return(merged)
  }

  for (i in 2:length(score_tables)) {
    merged <- merge(merged, score_tables[[i]], by = "IID", all = TRUE, sort = FALSE)
  }
  merged
}

coerce_binary_outcome <- function(x) {
  if (is.numeric(x)) {
    uniq <- sort(unique(stats::na.omit(x)))
    if (all(uniq %in% c(0, 1))) {
      return(as.numeric(x))
    }
  }

  x_chr <- tolower(trimws(as.character(x)))
  case_values <- c("1", "case", "yes", "true", "disease", "affected")
  ctrl_values <- c("0", "control", "no", "false", "healthy", "unaffected")

  out <- rep(NA_real_, length(x_chr))
  out[x_chr %in% ctrl_values] <- 0
  out[x_chr %in% case_values] <- 1

  if (any(is.na(out))) {
    stopf("Outcome column must be binary and encoded as 0/1 or control/case style labels.")
  }

  out
}

safe_numeric_frame <- function(dat, cols) {
  for (col in cols) {
    dat[[col]] <- as.numeric(dat[[col]])
    if (all(is.na(dat[[col]]))) {
      stopf("Predictor column became all NA after numeric coercion: %s", col)
    }
    missing_idx <- is.na(dat[[col]])
    if (any(missing_idx)) {
      dat[[col]][missing_idx] <- stats::median(dat[[col]], na.rm = TRUE)
    }
  }
  dat
}

sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

write_tsv <- function(dat, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(dat, file = path, sep = "\t", quote = FALSE, na = "NA")
}

load_model_coefficients <- function(path) {
  coefs <- fread(path, data.table = FALSE)
  names(coefs) <- trim_names(names(coefs))
  required <- c("feature", "coefficient")
  if (!all(required %in% names(coefs))) {
    stopf("Model coefficient file must contain columns: %s", paste(required, collapse = ", "))
  }
  coefs$feature <- as.character(coefs$feature)
  coefs$coefficient <- as.numeric(coefs$coefficient)
  coefs
}

usage <- function() {
  cat(
    paste(
      "Usage:",
      "Rscript AD-GMRS.script/AD-GMRS.R predict --score-dir <dir> --profile-reference <AD-GMRS.profile.txt> --model <coefficients.tsv> --out-prefix <prefix> --risk-output <risk_scores.tsv>",
      "Rscript AD-GMRS.script/AD-GMRS.R train-model --input <training_matrix.tsv> --outcome-col <label_column> --model-out <coefficients.tsv>",
      sep = "\n"
    ),
    "\n"
  )
}

command_predict <- function(opts) {
  require_args(opts, c("score-dir", "profile-reference", "model", "out-prefix", "risk-output"))

  score_dir <- opts[["score-dir"]]
  profile_reference <- opts[["profile-reference"]]
  model_path <- opts[["model"]]
  out_prefix <- opts[["out-prefix"]]
  risk_output <- opts[["risk-output"]]
  id_col <- if ("id-col" %in% names(opts)) opts[["id-col"]] else "IID"

  score_files <- list.files(score_dir, pattern = "\\.(profile|sscore)$", full.names = TRUE)
  if (length(score_files) == 0) {
    stopf("No .profile or .sscore files found in %s", score_dir)
  }

  ref <- load_profile_reference(profile_reference)
  score_tables <- vector("list", length(score_files))

  for (i in seq_along(score_files)) {
    feature_id <- sub("\\.(profile|sscore)$", "", basename(score_files[[i]]))
    score_dat <- read_score_file(score_files[[i]])
    names(score_dat)[names(score_dat) == "raw_score"] <- feature_id
    score_tables[[i]] <- score_dat
  }

  raw_matrix <- merge_score_tables(score_tables)
  feature_cols <- setdiff(names(raw_matrix), "IID")
  ordered_feature_cols <- ref$feature_id[ref$feature_id %in% feature_cols]
  missing_in_reference <- setdiff(feature_cols, ref$feature_id)
  if (length(missing_in_reference) > 0) {
    messagef(
      "Warning: %d score files are missing from the profile reference and will be dropped.",
      length(missing_in_reference)
    )
  }
  if (length(ordered_feature_cols) == 0) {
    stopf("None of the score files matched the profile reference.")
  }

  raw_matrix <- raw_matrix[, c("IID", ordered_feature_cols), drop = FALSE]
  std_matrix <- raw_matrix

  for (feature_id in ordered_feature_cols) {
    ref_row <- ref[ref$feature_id == feature_id, , drop = FALSE]
    std_matrix[[feature_id]] <- standardize_feature(
      values = as.numeric(raw_matrix[[feature_id]]),
      mean_value = ref_row$ref_mean[[1]],
      scale_value = ref_row$ref_scale[[1]]
    )
  }

  feature_summary <- merge(
    data.frame(feature_id = ordered_feature_cols, stringsAsFactors = FALSE),
    ref,
    by = "feature_id",
    all.x = TRUE,
    sort = FALSE
  )

  write_tsv(raw_matrix, paste0(out_prefix, ".raw.tsv"))
  write_tsv(std_matrix, paste0(out_prefix, ".standardized.tsv"))
  write_tsv(feature_summary, paste0(out_prefix, ".feature_summary.tsv"))

  dat <- std_matrix
  if (!id_col %in% names(dat)) {
    stopf("ID column not found in standardized matrix: %s", id_col)
  }

  coefs <- load_model_coefficients(model_path)
  intercept_idx <- coefs$feature %in% c("(Intercept)", "Intercept")
  intercept <- if (any(intercept_idx)) sum(coefs$coefficient[intercept_idx]) else 0
  feature_coefs <- coefs[!intercept_idx & coefs$coefficient != 0, , drop = FALSE]
  if (nrow(feature_coefs) == 0) {
    stopf("Model coefficient file contains no non-zero feature weights.")
  }

  missing_features <- setdiff(feature_coefs$feature, names(dat))
  if (length(missing_features) > 0) {
    messagef(
      "Warning: %d model features are missing from the input matrix. They will be filled with 0.",
      length(missing_features)
    )
    for (feature_id in missing_features) {
      dat[[feature_id]] <- 0
    }
  }

  model_feature_order <- feature_coefs$feature
  dat <- safe_numeric_frame(dat, model_feature_order)
  x <- as.matrix(dat[, model_feature_order, drop = FALSE])
  beta <- feature_coefs$coefficient
  linear_score <- intercept + as.numeric(x %*% beta)
  risk_probability <- sigmoid(linear_score)

  risk_dat <- data.frame(
    IID = dat[[id_col]],
    AD_GMRS_linear_score = linear_score,
    AD_GMRS_risk_probability = risk_probability,
    stringsAsFactors = FALSE
  )

  write_tsv(risk_dat, risk_output)
  messagef("Wrote raw matrix: %s.raw.tsv", out_prefix)
  messagef("Wrote standardized matrix: %s.standardized.tsv", out_prefix)
  messagef("Risk scores written to %s", risk_output)
}

command_train_model <- function(opts) {
  require_args(opts, c("input", "outcome-col", "model-out"))

  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stopf("Package 'glmnet' is required only when rebuilding the training model.")
  }

  input_path <- opts[["input"]]
  outcome_col <- opts[["outcome-col"]]
  model_out <- opts[["model-out"]]
  id_col <- if ("id-col" %in% names(opts)) opts[["id-col"]] else "IID"
  exclude_cols <- if ("exclude-cols" %in% names(opts)) {
    trimws(strsplit(opts[["exclude-cols"]], ",", fixed = TRUE)[[1]])
  } else {
    character()
  }
  alpha <- if ("alpha" %in% names(opts)) as.numeric(opts[["alpha"]]) else 1
  nfolds <- if ("nfolds" %in% names(opts)) as.integer(opts[["nfolds"]]) else 10L
  seed <- if ("seed" %in% names(opts)) as.integer(opts[["seed"]]) else 20260803L
  lambda_rule <- if ("lambda-rule" %in% names(opts)) opts[["lambda-rule"]] else "1se"

  if (!lambda_rule %in% c("1se", "min")) {
    stopf("--lambda-rule must be either '1se' or 'min'.")
  }

  train_dat <- fread(input_path, data.table = FALSE)
  if (!outcome_col %in% names(train_dat)) {
    stopf("Outcome column not found: %s", outcome_col)
  }

  y <- coerce_binary_outcome(train_dat[[outcome_col]])
  drop_cols <- unique(c(outcome_col, exclude_cols, id_col))
  predictor_cols <- setdiff(names(train_dat), drop_cols)
  if (length(predictor_cols) == 0) {
    stopf("No predictor columns remain after excluding id/outcome columns.")
  }

  train_dat <- safe_numeric_frame(train_dat, predictor_cols)
  sd_values <- vapply(train_dat[predictor_cols], stats::sd, numeric(1), na.rm = TRUE)
  zero_var_cols <- names(sd_values)[is.na(sd_values) | sd_values == 0]
  predictor_cols <- setdiff(predictor_cols, zero_var_cols)
  if (length(predictor_cols) == 0) {
    stopf("All predictors were removed because of zero variance.")
  }

  x <- as.matrix(train_dat[, predictor_cols, drop = FALSE])
  set.seed(seed)
  fit <- glmnet::cv.glmnet(
    x = x,
    y = y,
    family = "binomial",
    alpha = alpha,
    nfolds = nfolds,
    standardize = FALSE
  )

  chosen_lambda <- if (lambda_rule == "min") fit$lambda.min else fit$lambda.1se
  coef_matrix <- as.matrix(stats::coef(fit, s = chosen_lambda))
  coef_dt <- data.frame(
    feature = rownames(coef_matrix),
    coefficient = as.numeric(coef_matrix[, 1]),
    stringsAsFactors = FALSE
  )

  selected_dt <- coef_dt[coef_dt$feature == "(Intercept)" | coef_dt$coefficient != 0, , drop = FALSE]
  summary_dt <- data.frame(
    item = c(
      "n_samples",
      "n_predictors_input",
      "n_predictors_zero_variance_removed",
      "n_predictors_selected",
      "alpha",
      "nfolds",
      "lambda_rule",
      "lambda_value"
    ),
    value = c(
      nrow(train_dat),
      length(setdiff(names(train_dat), drop_cols)),
      length(zero_var_cols),
      sum(selected_dt$feature != "(Intercept)"),
      alpha,
      nfolds,
      lambda_rule,
      chosen_lambda
    ),
    stringsAsFactors = FALSE
  )

  write_tsv(selected_dt, model_out)
  write_tsv(summary_dt, sub("\\.tsv$", ".summary.tsv", model_out))
  messagef("Selected %d predictors.", sum(selected_dt$feature != "(Intercept)"))
  messagef("Model weights written to %s", model_out)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || args[[1]] %in% c("help", "-h", "--help")) {
  usage()
  quit(save = "no", status = 0)
}

subcommand <- args[[1]]
opts <- parse_args(args[-1])

if (subcommand == "predict") {
  command_predict(opts)
} else if (subcommand == "train-model") {
  command_train_model(opts)
} else {
  stopf("Unknown subcommand: %s", subcommand)
}
