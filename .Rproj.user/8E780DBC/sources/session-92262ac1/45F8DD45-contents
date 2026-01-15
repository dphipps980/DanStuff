#' Score AMP Data
#'
#' @param data A data frame containing AMP data
#' @param columns Column indices or names containing AMP CSV data
#' @param id_col Name of ID column (default: "ID")
#' @param prefix Prefix for output column names (default: "")
#' @param quality Include quality checks (default: TRUE)
#' @param debug Print debug messages (default: FALSE)
#' @param target1 First target category (e.g., "Physical Activity")
#' @param target2 Second target category (e.g., "Sedentary Behaviour")
#' @param positive Positive response labels (e.g., c("Pleasant", "Miellyttävä"))
#' @param negative Negative response labels (e.g., c("Unpleasant", "Epämiellyttävä"))
#'
#' @return Data frame with AMP scores
#' @export
#'
#' @examples
#' \dontrun{
#' scores <- AMPScore(data, c(2:5),
#'                    target1 = "Physical Activity",
#'                    target2 = "Sedentary Behaviour",
#'                    positive = c("Pleasant", "Miellyttävä"),
#'                    negative = c("Unpleasant", "Epämiellyttävä"))
#' }

AMPScore <- function(data, columns, id_col = "ID", prefix = "", quality = TRUE, debug = FALSE,
                     target1 = NULL, target2 = NULL,
                     positive = NULL, negative = NULL) {

  # Check that all required parameters are provided
  if (is.null(target1) || is.null(target2)) {
    stop("target1 and target2 must be specified")
  }
  if (is.null(positive) || is.null(negative)) {
    stop("positive and negative response labels must be specified")
  }

  # Create short labels from target names (first 3 characters, uppercase)
  t1_label <- toupper(substr(gsub(" ", "", target1), 1, 3))
  t2_label <- toupper(substr(gsub(" ", "", target2), 1, 3))

  # Create column names with optional prefix
  col_t1 <- paste0(prefix, "AMP.", t1_label)
  col_t2 <- paste0(prefix, "AMP.", t2_label)
  col_diff <- paste0(prefix, "AMP.diff")
  col_ntrial <- paste0(prefix, "AMP.NTrial")
  col_source <- paste0(prefix, "AMP.source")
  col_rt <- paste0(prefix, "AMP.RT")
  col_samekey <- paste0(prefix, "AMP.SameKey")

  # Initialize output dataframe with just ID and scores
  output <- data.frame(id = data[[id_col]])
  names(output)[1] <- id_col
  output[[col_t1]] <- NA
  output[[col_t2]] <- NA
  output[[col_diff]] <- NA
  output[[col_ntrial]] <- NA

  # Add quality columns if requested
  if (quality) {
    output[[col_rt]] <- NA
    output[[col_samekey]] <- NA
  }

  # Only add source column if multiple columns are provided
  include_source <- length(columns) > 1
  if (include_source) {
    output[[col_source]] <- NA
  }

  # Create regex patterns for matching
  positive_pattern <- paste(positive, collapse = "|")
  negative_pattern <- paste(negative, collapse = "|")

  # Process each row
  for (iRow in 1:nrow(data)) {
    # Find the non-NA column for this row
    row_data <- data[iRow, columns, drop = TRUE]
    non_na_cols <- which(!is.na(row_data) & row_data != '')

    # Skip if no valid data
    if (length(non_na_cols) == 0) {
      next
    }
    if (length(non_na_cols) > 1) {
      warning(paste("Row", iRow, "has multiple non-NA columns. Using first one."))
    }

    # Get the column to process
    col_idx <- columns[non_na_cols[1]]
    col_name <- names(data)[col_idx]

    # Record source column if multiple columns provided
    if (include_source) {
      output[[col_source]][iRow] <- col_name
    }

    # Extract data - properly handle tibble
    csv_string <- as.character(data[[col_idx]][iRow])

    # Skip if doesn't contain 'block'
    if (!grepl('block', csv_string)) {
      if (debug) message(paste('Row', iRow, '- no "block" found'))
      next
    }

    # Clean up quotes
    csv_string <- gsub('""', '"', csv_string)

    # Try to parse CSV
    df2 <- tryCatch({
      read.csv(text = csv_string, stringsAsFactors = FALSE)
    },
    error = function(err) {
      message(paste('Row', iRow, '- malformed CSV in column', col_name, ':', err$message))
      return(NULL)
    })

    # Skip if parsing failed or empty
    if (is.null(df2) || nrow(df2) == 0) {
      if (debug) message(paste('Row', iRow, '- empty dataframe'))
      next
    }

    if (debug) message(paste('Row', iRow, '- parsed', nrow(df2), 'rows'))

    # Check if required columns exist
    required_cols <- c("cond", "resp")
    if (quality) required_cols <- c(required_cols, "rt")

    if (!all(required_cols %in% names(df2))) {
      message(paste('Row', iRow, '- missing required columns. Has:', paste(names(df2), collapse=", ")))
      next
    }

    # Process AMP data - identify targets
    df2$prime <- ifelse(grepl(target1, df2$cond, fixed = TRUE), target1,
                        ifelse(grepl(target2, df2$cond, fixed = TRUE), target2, NA))

    if (debug) message(paste('Row', iRow, '- prime values:', paste(table(df2$prime, exclude=NULL), collapse=", ")))

    # Rating: 1 = positive, 0 = negative
    df2$rating <- ifelse(grepl(negative_pattern, df2$resp), 0,
                         ifelse(grepl(positive_pattern, df2$resp), 1, NA))

    if (debug) message(paste('Row', iRow, '- ratings:', paste(table(df2$rating, exclude=NULL), collapse=", ")))

    # Filter to correct trials
    df2.correct <- df2[which(df2$prime %in% c(target1, target2)), ]

    # Skip if no valid trials
    if (nrow(df2.correct) == 0) {
      message(paste('Row', iRow, '- no valid trials found after filtering'))
      next
    }

    # Store number of trials
    output[[col_ntrial]][iRow] <- nrow(df2.correct)

    if (debug) message(paste('Row', iRow, '- valid trials:', nrow(df2.correct)))

    # Quality checks
    if (quality) {
      # Average RT for valid trials
      output[[col_rt]][iRow] <- mean(df2.correct$rt, na.rm = TRUE)

      # Check if always same key (all 0s or all 1s)
      unique_ratings <- unique(df2.correct$rating[!is.na(df2.correct$rating)])
      output[[col_samekey]][iRow] <- length(unique_ratings) == 1
    }

    # Calculate means for this person
    df2.correct$id <- data[[id_col]][iRow]
    amp.means <- summaryBy(formula = rating ~ prime, data = df2.correct)

    if (debug) {
      message(paste('Row', iRow, '- amp.means:'))
      print(amp.means)
    }

    # Store results directly from amp.means
    t1_idx <- which(amp.means$prime == target1)
    t2_idx <- which(amp.means$prime == target2)

    if (length(t1_idx) > 0) {
      output[[col_t1]][iRow] <- amp.means$rating.mean[t1_idx]
    }
    if (length(t2_idx) > 0) {
      output[[col_t2]][iRow] <- amp.means$rating.mean[t2_idx]
    }

    # Calculate difference
    if (!is.na(output[[col_t1]][iRow]) && !is.na(output[[col_t2]][iRow])) {
      output[[col_diff]][iRow] <- output[[col_t1]][iRow] - output[[col_t2]][iRow]
    }
  }

  return(output)
}

