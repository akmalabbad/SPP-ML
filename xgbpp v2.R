# =============================================================
# BAGIAN 1: DEFINISI FUNGSI-FUNGSI OBJECTIVE & METRIC
# =============================================================

# --- Objective & Metric untuk 'poisson' dan 'weighted_poisson' ---
make_poisson_obj <- function(vol_vec) {
  function(preds, dtrain) {
    labels <- getinfo(dtrain, "label"); weights <- getinfo(dtrain, "weight")
    if (is.null(weights)) weights <- rep(1, length(labels))
    dummy <- which(labels < 0); event <- which(labels > 0)
    norm <- if (length(dummy) > 0) sum(weights[dummy] * vol_vec[dummy]) else 1
    if (abs(norm) < 1e-9) norm <- 1
    weights <- weights / norm
    grad <- numeric(length(labels)); hess <- numeric(length(labels))
    grad[event] <- -weights[event]; hess[event] <- 1e-6
    if (length(dummy) > 0) {
      mu <- exp(preds[dummy]); gh <- weights[dummy] * mu * vol_vec[dummy]
      grad[dummy] <- gh; hess[dummy] <- pmax(gh, 1e-6)
    }
    return(list(grad = grad, hess = hess))
  }
}
make_weighted_poisson_obj <- function(vol_vec, F_prime_val, hess_min = 1e-6) {
  force(vol_vec); force(F_prime_val); force(hess_min)
  function(preds, dtrain) {
    label  <- getinfo(dtrain, "label")
    weight <- 1 / (1 + exp(preds) * F_prime_val)
    weight <- weight / sum(weight[label < 0] * vol_vec[label < 0])
    pos <- label > 0; neg <- !pos
    grad <- numeric(length(label)); hess <- numeric(length(label))
    grad[pos] <- -1; grad[neg] <- exp(preds[neg]) * vol_vec[neg]; grad <- weight * grad
    hess[pos] <- 0; hess[neg] <- exp(preds[neg]) * vol_vec[neg]; hess <- pmax(weight * hess, hess_min)
    list(grad = grad, hess = hess)
  }
}
make_poisson_metric <- function(vol_vec) {
  .metric <- function(preds, dtrain) {
    labels <- getinfo(dtrain, "label"); event <- which(labels > 0); dummy <- which(labels < 0)
    err <- numeric(length(labels))
    if (length(event) > 0) err[event] <- -preds[event]
    if (length(dummy) > 0) err[dummy] <- exp(preds[dummy]) * vol_vec[dummy]
    list(name = "poisson_err", value = sum(err) / 1e6)
  }
  return(.metric)
}

# --- Objective & Metric untuk 'logistic' ---
# Menggunakan objective yang telah diperbaiki dan metrik persis seperti yang Anda kirim.
make_logistic_obj <- function(vol_vec, F_prime_val, hess_min = 1e-6) {
  force(vol_vec); force(F_prime_val); force(hess_min)
  function(preds, dtrain) {
    label  <- getinfo(dtrain, "label")
    weight <- 1 / (1 + exp(preds) * F_prime_val)
    weight <- weight / sum(weight[label < 0] * vol_vec[label < 0])
    pos <- label > 0; neg <- !pos
    grad <- numeric(length(label)); hess <- numeric(length(label))
    delta_pos <- 1; delta_neg <- 1 / vol_vec[neg]
    exp_preds_pos <- exp(preds[pos]); exp_preds_neg <- exp(preds[neg])
    grad[pos] <- (-weight[pos] * delta_pos) / (delta_pos + exp_preds_pos)
    grad[neg] <- (weight[neg] * exp_preds_neg) / (exp_preds_neg + delta_neg)
    hess[pos] <- (weight[pos] * delta_pos * exp_preds_pos) / ((delta_pos + exp_preds_pos)^2)
    hess[neg] <- (weight[neg] * delta_neg * exp_preds_neg) / ((delta_neg + exp_preds_neg)^2)
    hess <- pmax(hess, hess_min)
    list(grad = grad, hess = hess)
  }
}

#' Membuat evaluation metric logistik (Sesuai Snippet Asli Anda)
make_logistic_metric <- function(vol_vec) {
  .metric <- function(preds, dtrain) {
    labels <- getinfo(dtrain, "label")
    event <- which(labels > 0)
    dummy <- which(labels < 0)
    err <- numeric(length(labels))
    if (length(event) > 0) err[event] <- -log(exp(preds[event])/(vol_vec[event]+exp(preds[event])))
    if (length(dummy) > 0) err[dummy] <- -log(exp(vol_vec[dummy])/(exp(preds[dummy])+vol_vec[dummy]))
    
    # Menggunakan nama "logistic_err" untuk konsistensi
    return(list(name = "logistic_err", value = sum(err) / 1e6))
  }
  return(.metric)
}


# =============================================================
# BAGIAN 2: FUNGSI WRAPPER UTAMA "xgbpp"
# =============================================================

xgbpp <- function(vol, params, ..., loss = "poisson", F_prime = 1) {
  
  # --- 1. Pilih Objective Function ---
  objective_func <- switch(
    loss,
    "poisson"          = make_poisson_obj(vol),
    "weighted_poisson" = make_weighted_poisson_obj(vol, F_prime),
    "logistic"         = make_logistic_obj(vol, F_prime),
    stop("Loss '", loss, "' tidak dikenali.")
  )
  
  # --- 2. Pilih Evaluation Metric ---
  eval_func <- switch(
    loss,
    "poisson"          = make_poisson_metric(vol),
    "weighted_poisson" = make_poisson_metric(vol),
    "logistic"         = make_logistic_metric(vol),
    stop("Loss '", loss, "' tidak dikenali.")
  )
  
  # --- 3. Penetapan Parameter & Pemanggilan xgb.train ---
  params$objective <- objective_func
  
  model <- xgb.train(
    params = params,
    feval = eval_func,
    maximize = FALSE,
    ...
  )
  
  return(model)
}