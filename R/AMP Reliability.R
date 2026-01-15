#' Calculate AMP Reliability
#'
#' Calculate split-half reliability for Affect Misattribution Procedure (AMP) data
#' using odd-even, first-second half, or random splits with Spearman-Brown correction.
#'
#' @param data A data frame containing AMP data
#' @param columns Column indices or names containing AMP CSV data
#' @param id_col Name of ID column (default: "ID")
#' @param target1 First target category (e.g., "Physical Activity")
#' @param target2 Second target category (e.g., "Sedentary Behaviour")
#' @param positive Positive response labels (e.g., c("Pleasant", "Miellyttävä"))
#' @param negative Negative response labels (e.g., c("Unpleasant", "Epämiellyttävä"))
#' @param method Split method: "odd_even", "first_second", or "random" (default: "odd_even")
#' @param n_splits Number of random splits when method = "random" (default: 100)
#' @param seed Random seed for reproducibility (default: NULL)
#'
#' @return Data frame with split-half correlations, Spearman-Brown corrected reliability,
#'   sample size, and method. If multiple columns are provided, includes per-column reliability.
#' @export
#'
#' @examples
#' \dontrun{
#' # Odd-even split
#' rel1 <- AMPReliability(data, c(2:5),
#'                        target1 = "Physical Activity",
#'                        target2 = "Sedentary Behaviour",
#'                        positive = c("Pleasant", "Miellyttävä"),
#'                        negative = c("Unpleasant", "Epämiellyttävä"),
#'                        method = "odd_even")
#'
#' # Random splits (most robust)
#' rel2 <- AMPReliability(data, c(2:5),
#'                        target1 = "Physical Activity",
#'                        target2 = "Sedentary Behaviour",
#'                        positive = c("Pleasant", "Miellyttävä"),
#'                        negative = c("Unpleasant", "Epämiellyttävä"),
#'                        method = "random",
#'                        n_splits = 100,
#'                        seed = 123)
#' }
AMPReliability <- function(data, columns, id_col = "ID",
                           target1 = NULL, target2 = NULL,
                           positive = NULL, negative = NULL,
                           method = "odd_even", n_splits = 100, seed = NULL) {

  # Check that all required parameters are provided
  if (is.null(target1) || is.null(target2)) {
    stop("target1 and target2 must be specified")
  }
  if (is.null(positive) || is.null(negative)) {
    stop("positive and negative response labels must be specified")
  }

  # Check method
  valid_methods <- c("odd_even", "first_second", "random")
  if (!method %in% valid_methods) {
    stop(paste("method must be one of:", paste(valid_methods, collapse = ", ")))
  }

  if (!is.null(seed)) set.seed(seed)

  # Collect all trial data across participants
  all_trials <- list()

  for (iRow in 1:nrow(data)) {
    # Find the non-NA column for this row
    row_data <- data[iRow, columns, drop = TRUE]
    non_na_cols <- which(!is.na(row_data) & row_data != '')

    if (length(non_na_cols) == 0) next

    # Get the column to process
    col_idx <- columns[non_na_cols[1]]
    col_name <- names(data)[col_idx]
    csv_string <- as.character(data[[col_idx]][iRow])

    if (!grepl('block', csv_string)) next

    # Parse CSV
    csv_string <- gsub('""', '"', csv_string)
    df2 <- tryCatch({
      read.csv(text = csv_string, stringsAsFactors = FALSE)
    }, error = function(err) NULL)

    if (is.null(df2) || nrow(df2) == 0) next
    if (!all(c("cond", "resp", "trial") %in% names(df2))) next

    # Process AMP data
    df2$prime <- ifelse(grepl(target1, df2$cond, fixed = TRUE), target1,
                        ifelse(grepl(target2, df2$cond, fixed = TRUE), target2, NA))

    # Create regex patterns for matching
    positive_pattern <- paste(positive, collapse = "|")
    negative_pattern <- paste(negative, collapse = "|")

    df2$rating <- ifelse(grepl(negative_pattern, df2$resp), 0,
                         ifelse(grepl(positive_pattern, df2$resp), 1, NA))

    # Filter to correct trials
    df2.correct <- df2[which(df2$prime %in% c(target1, target2)), ]

    if (nrow(df2.correct) == 0) next

    df2.correct$ID <- data[[id_col]][iRow]
    df2.correct$source_col <- col_name
    all_trials[[iRow]] <- df2.correct
  }

  # Combine all data
  amp_data <- do.call(rbind, all_trials)

  # Spearman-Brown correction function
  sb_reliability <- function(r) {
    (2 * r) / (1 + r)
  }

  # Function to calculate reliability based on method
  calc_reliability <- function(data_subset) {

    if (method == "odd_even") {
      # Split trials into odd and even
      data_subset$split <- ifelse(data_subset$trial %% 2 == 1, "odd", "even")
      split_names <- c("odd", "even")

    } else if (method == "first_second") {
      # For each person, split their trials into first half and second half
      data_subset <- data_subset[order(data_subset$ID, data_subset$trial), ]
      data_subset <- data_subset %>%
        group_by(ID) %>%
        mutate(
          trial_rank = rank(trial),
          n_trials = n(),
          split = ifelse(trial_rank <= n_trials/2, "first", "second")
        ) %>%
        ungroup()
      split_names <- c("first", "second")

    } else if (method == "random") {
      # Run multiple random splits and average
      cors <- numeric(n_splits)

      for (i in 1:n_splits) {
        # For each person, randomly assign their trials to half A or half B
        data_subset_split <- data_subset %>%
          group_by(ID) %>%
          mutate(
            split = sample(rep(c("A", "B"), length.out = n()))
          ) %>%
          ungroup()

        # Calculate means for each split
        split.means <- summaryBy(rating ~ ID + prime + split, data = as.data.frame(data_subset_split))

        # Reshape to wide format
        wide <- split.means %>%
          tidyr::pivot_wider(names_from = split,
                      values_from = rating.mean,
                      id_cols = c(ID, prime))

        # Calculate D-score reliability
        wide_diff <- wide %>%
          group_by(ID) %>%
          summarize(A_diff = A[prime == target1] - A[prime == target2],
                    B_diff = B[prime == target1] - B[prime == target2])

        cors[i] <- cor(wide_diff$A_diff, wide_diff$B_diff, use = "complete.obs")
      }

      # Average correlation across all random splits
      mean_cor <- mean(cors, na.rm = TRUE)
      rel_diff <- sb_reliability(mean_cor)

      return(data.frame(
        Split_Half_r = round(mean_cor, 3),
        Spearman_Brown = round(rel_diff, 3),
        N = nrow(wide_diff),
        Method = method,
        N_splits = n_splits
      ))
    }

    # For odd_even and first_second methods
    if (method != "random") {
      # Calculate means for each split
      split.means <- summaryBy(rating ~ ID + prime + split, data = as.data.frame(data_subset))

      # Reshape to wide format
      wide <- split.means %>%
        tidyr::pivot_wider(names_from = split,
                    values_from = rating.mean,
                    id_cols = c(ID, prime))

      # Calculate D-score reliability
      wide_diff <- wide %>%
        group_by(ID) %>%
        summarize(split1_diff = .data[[split_names[1]]][prime == target1] - .data[[split_names[1]]][prime == target2],
                  split2_diff = .data[[split_names[2]]][prime == target1] - .data[[split_names[2]]][prime == target2])

      cor_diff <- cor(wide_diff$split1_diff, wide_diff$split2_diff, use = "complete.obs")
      rel_diff <- sb_reliability(cor_diff)

      return(data.frame(
        Split_Half_r = round(cor_diff, 3),
        Spearman_Brown = round(rel_diff, 3),
        N = nrow(wide_diff),
        Method = method
      ))
    }
  }

  # Calculate overall reliability
  overall <- calc_reliability(amp_data)
  overall$Source <- "Overall"

  results <- overall

  # If multiple columns, calculate per-column reliability
  if (length(columns) > 1) {
    source_cols <- unique(amp_data$source_col)

    for (col in source_cols) {
      col_data <- amp_data[amp_data$source_col == col, ]
      if (nrow(col_data) > 0) {
        col_rel <- calc_reliability(col_data)
        col_rel$Source <- col
        results <- rbind(results, col_rel)
      }
    }
  }

  # Reorder columns
  if (method == "random") {
    results <- results[, c("Source", "Split_Half_r", "Spearman_Brown", "N", "Method", "N_splits")]
  } else {
    results <- results[, c("Source", "Split_Half_r", "Spearman_Brown", "N", "Method")]
  }
  rownames(results) <- NULL

  return(results)
}
