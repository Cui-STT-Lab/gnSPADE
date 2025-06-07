weightedLDA <- function(docs, Weight_Mat, neiPath, model,
                        number_of_topics, model_settings = list(),
                        priors = list(), options = list(), keep = c(), Z)
{
  stData <- keyATM_read(texts = quanteda::as.dfm(docs), keep_docnames = TRUE)

  # Check type
  if (length(keep) != 0)
    check_arg_type(keep, "character")

  model_name <- full_model_name(model, type = "lda")
  if (is.null(options$seed))
    options$seed <- floor(stats::runif(1) * 1e5)
  set.seed(options$seed)

 # Check if there is a resume object
  if ("resume" %in% names(options) && fs::file_exists(options$resume)) {
    resume <- TRUE
    fitted <- fitted_load(options$resume)
    exists_iter <- fitted$used_iter
    fitted <- fitted_update_iterations(fitted, options)
    fitted <- keyATM_fit(fitted, resume = TRUE)
    used_iter <- get_used_iter(fitted, resume, exists = exists_iter)
  } else {
    resume <- FALSE
    initialized <- keyATM_initialize(
      stData, model_name, number_of_topics,
      keywords = list(), model_settings = model_settings,
      priors = priors, options = options, Z = Z
    )

    initialized$model$Weight_Mat = Weight_Mat[,initialized$model$vocab]
    word_neighbors = load_wordid_neighcnt(neiPath, initialized$model$vocab)
    initialized$model$word_neighbors = word_neighbors
    #initialized$model$neiPath = neiPath

    fitted <- keyATM_fit(initialized)
    used_iter <- get_used_iter(fitted, resume)
  }

  if ("resume" %in% names(options)) {
    fitted_save(options$resume, fitted, model_name, used_iter)
  }

  # 0 iterations
  if (fitted$options$iterations == 0) {
    cli::cli_alert_info("`options$iterations` is 0. keyATM returns an initialized object.")
    return(fitted)
  }

  # Get output
  out <- keyATM_output(fitted, keep, used_iter)
  out$number_of_topics <- number_of_topics
  out$no_keyword_topics <- NULL
  out$keyword_k <- NULL

  return(out)
}

# Resume related functions
fitted_load <- function(filename) {
  fitted <- readRDS(filename)
  if (! "keyATM_resume" %in% class(fitted))
    cli::cli_abort("The file is not a keyATM object.")
  cli::cli_alert_success("The fitted model is loaded from {.file {filename}}.")
  return(fitted)
}

fitted_update_iterations <- function(fitted, new_options) {
  if (! "keyATM_resume" %in% class(fitted))
    cli::cli_abort("The file is not a keyATM object.")
  if (! "iterations" %in% names(new_options))
    cli::cli_abort("The `options` argument must contain `iterations`.")
  fitted$model$options$iterations <- fitted$model$options$iterations + new_options$iterations  # total after the fitting
  fitted$model$options$iter_new <- new_options$iterations  # iterations to add
  return(fitted)
}

fitted_save <- function(filename, fitted, model_name, used_iter) {
  saveobj <- list(
    model = fitted,
    model_name = model_name,
    used_iter = used_iter,
    rand_state = .GlobalEnv$.Random.seed
  )
  class(saveobj) <- c("keyATM_resume", class(saveobj))
  saveRDS(object = saveobj, file = filename)
  cli::cli_alert_success("The fitted model is saved in {.file {filename}}.")
}

# Get indexes of used iterations
get_used_iter <- function(fitted, resume, exists = NULL) {
  thinning <- fitted$options$thinning
  if (resume) {
    exists_max <- max(exists)
    total_iter <- (exists_max + 1):(fitted$options$iterations)
    total_iter <- total_iter[(total_iter %% thinning == 0) | (total_iter == 1) | total_iter == max(total_iter)]
    used_iter <- c(exists, total_iter)
  } else {
    total_iter <- 1:(fitted$options$iterations)
    used_iter <- total_iter[(total_iter %% thinning == 0) | (total_iter == 1) | total_iter == max(total_iter)]
  }
  return(used_iter)
}
