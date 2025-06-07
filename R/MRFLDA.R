#' MRF-LDA
#'
#' This function do the deconvolution with genes network.
#' @export
WRLDA = function(corpus, k, type='base', seed=0, iterations=1500, priors=NULL, lambda = NULL,
                beta_init=NULL, ldamodel=NULL, verbose = FALSE, neiPath=NULL, threshold = 1,
                Z = NULL){

  if (!"alpha" %in% names(priors)) {
    priors$alpha = rep(1/k,k)
  }
  if (!"beta" %in% names(priors)) {
    priors$beta = 0.01
  }
  if(is.null(lambda)){
    lambda = 1
  }
  if (!"lambda" %in% names(priors)) {
    priors$lambda = lambda
  }

  Weight_Mat = Weight_Term(corpus, k, type=type,
                           seed=seed, beta_init=beta_init, ldamodel=ldamodel)
  #NeighborNodes
  if(is.null(neiPath)){
    neiPath = GenNeiPath(corpus = corpus, threshold = threshold)
  }

  out = weightedLDA(docs = corpus, Weight_Mat = Weight_Mat, neiPath = neiPath,
                    model = "base", number_of_topics = k,
                    options = list(seed = seed, iterations = iterations, verbose = verbose,
                                   estimate_alpha = 0, weights_type = type, store_theta = T,
                                   store_phi = T),
                    priors = priors, Z = Z)

  if (ncol(corpus) > length(out$vocab)){
    word = setdiff(colnames(corpus), out$vocab)
    out$vocab = append(out$vocab, word)
    miss_word <- matrix(10^-40, nrow = nrow(out$phi), ncol = length(word),
                        dimnames = list(NULL, word))
    out$phi <- cbind(out$phi, miss_word)
  }

  out$perplexity = Perplexity_wlda(out, corpus)
  return(out)
}

Weight_Term = function(corpus, k, type=NULL, seed=0, beta_init=NULL, ldamodel=NULL){
  if (is.null(type)) {
    type <- "base"
  } else {
    if (!type %in% c('base','infor_bdc', 'inverse_bdc', 'bdc', 'infor',
                                 'inverse', 'mutual_point', 'tf_idf'))
    {
      cli::cli_abort("An invalid value in `type`")
    }
  }

  if(type == 'base'){
    Weight_Mat= matrix(1, nrow = nrow(corpus), ncol = ncol(corpus), byrow = TRUE)
  }

  if(grepl("infor", type)){
    weight_log = -log2((colSums(corpus)+1)/sum(corpus))
    weight_log_norm = (weight_log-min(weight_log))/(max(weight_log)-min(weight_log))
  }
  if(grepl("inverse", type)){
    weight_inv = sum(corpus)/(colSums(corpus)+1)
    weight_inv_norm = (weight_inv-min(weight_inv))/(max(weight_inv)-min(weight_inv))
  }
  if(grepl("bdc", type)){

    if(is.null(beta_init)){
      if(!is.null(ldamodel)){
        beta = exp(ldamodel@beta)
      }else{
        ldamodel = topicmodels::LDA(x = slam::as.simple_triplet_matrix(corpus),
                                    k = k,
                                    control = list(seed=seed, verbose=0, keep=1, estimate.alpha=TRUE))
        beta = exp(ldamodel@beta)
      }
    }else{
      beta = beta_init
    }

    beta_col_sums <- colSums(beta)
    beta_normalized_mat <- sweep(beta, 2, beta_col_sums, FUN = "/")
    bdc <- apply(beta_normalized_mat, 2, function(col) sum(col * log(col)))
    weight_bdc = 1+bdc/log(nrow(beta))
    weight_bdc_norm = (weight_bdc-min(weight_bdc))/(max(weight_bdc)-min(weight_bdc))

  }

  if(type == 'infor'){
    Weight_Mat= matrix(weight_log_norm, nrow = nrow(corpus), ncol = ncol(corpus), byrow = TRUE)
  }
  if(type == 'inverse'){
    Weight_Mat= matrix(weight_inv_norm, nrow = nrow(corpus), ncol = ncol(corpus), byrow = TRUE)
  }
  if(type == 'bdc'){
    Weight_Mat= matrix(weight_bdc_norm, nrow = nrow(corpus), ncol = ncol(corpus), byrow = TRUE)
  }
  if(type == 'infor_bdc'){
    weight_log_bdc = weight_log_norm*weight_bdc_norm
    weight_log_bdc_norm = (weight_log_bdc-min(weight_log_bdc))/(max(weight_log_bdc)-min(weight_log_bdc))
    Weight_Mat= matrix(weight_log_bdc_norm, nrow = nrow(corpus), ncol = ncol(corpus), byrow = TRUE)
  }
  if(type == 'inverse_bdc'){
    weight_inv_bdc = weight_inv_norm*weight_bdc_norm
    weight_inv_bdc_norm = (weight_inv_bdc-min(weight_inv_bdc))/(max(weight_inv_bdc)-min(weight_inv_bdc))
    Weight_Mat= matrix(weight_inv_bdc_norm, nrow = nrow(corpus), ncol = ncol(corpus), byrow = TRUE)
  }

  if(type == 'mutual_point'){
    Weight_Mat = -log2(sweep(corpus+1, 2, colSums(corpus)+1, FUN = '/'))
  }

  if(type == 'tf_idf'){
    tf = -log2(sweep(corpus+1, 2, colSums(corpus)+1, FUN = '/'))
    idf = log(nrow(corpus)/(colSums(corpus>0)+1))
    Weight_Mat = sweep(tf, 2, idf, FUN = '*')
  }

  rownames(Weight_Mat) <- rownames(corpus)
  colnames(Weight_Mat) <- colnames(corpus)

  return(Weight_Mat)
}


Perplexity_wlda = function(model, corpus){
  beta = model$phi[, colnames(corpus)]
  theta = model$theta[rownames(corpus), ]
  x = slam::as.simple_triplet_matrix(corpus)
  perplexity = exp(-sum(log(colSums(beta[,x$j]*t(theta)[,x$i]))*x$v)/sum(x$v))
  # inner = log(theta%*%beta)
  # perplexity_check = exp(-sum(corpus*inner)/sum(corpus))
  return(perplexity)
}

GenNeiPath = function(corpus, threshold = 0.8){
  correlation = cor(corpus)

  neiPath <- lapply(1:nrow(correlation), function(i) {
    # Get the indices of entries > 0.5 for the i-th row
    neighbors <- which(correlation[i, ] > threshold)
    # Get the column names (genes) of these entries
    colnames(correlation)[neighbors]
  })

  # Set the names of neiPath to the rownames of the matrix
  names(neiPath) <- rownames(correlation)

  return(neiPath)
}

