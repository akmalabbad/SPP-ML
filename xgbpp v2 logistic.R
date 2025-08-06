#' Wrapper xgb.train dengan dua pilihan custom objective Poisson
#'
#' Fungsi ini menyederhanakan training XGBoost dengan memilih antara objective standar
#' atau objective dengan pembobotan dinamis. Formatnya meniru lgbpp.
#'
#' @param vol Vektor numerik untuk custom loss.
#' @param params List parameter untuk XGBoost (gunakan nama parameter XGBoost seperti 'eta').
#' @param dynamic_weighted Jika TRUE, gunakan objective dengan pembobotan dinamis.
#' @param F_prime Nilai yang hanya digunakan jika dynamic_weighted = TRUE.
#' @param ... Argumen lain yang diteruskan ke `xgb.train` (misalnya, data, nrounds).
#'
#' @return Objek `xgb.Booster` yang sudah dilatih.

xgbpp <- function(vol, params, ..., dynamic_weighted = FALSE, F_prime = 1) {
  
  # 1. --- Pilih Objective Function Berdasarkan Argumen ---
  
  if (dynamic_weighted) {
    # --- Mode: Dynamic Weighted ---
    cat("--- Menggunakan objective function dengan Dynamic Weights (F_prime) ---\n")
    
    make_dw_logistic_obj <- function(vol_vec, F_prime_val, hess_min = 1e-6) {
      force(vol_vec); force(F_prime_val); force(hess_min)
      
      function(preds, dtrain) {
        label  <- getinfo(dtrain, "label") # getinfo() juga standar untuk XGBoost
        
        weight <- 1 / (1 + exp(preds) * F_prime_val)
        weight <- weight / sum(weight[label < 0] * vol_vec[label < 0])
        
        pos <- label > 0
        neg <- !pos
        
        grad <- numeric(length(label))
        hess <- numeric(length(label))
        
        delta_pos = 1/vol_vec[pos]
        delta_neg = 1/vol_vec[neg]
        
        grad[pos] <- (-weight*delta_pos)/(delta_pos+exp(preds[pos]))
        grad[neg] <- (weight*preds[neg])/(exp(preds[neg])+delta_neg)
        
        hess[pos] <- (weight*delta_pos*exp(preds[pos])) / ((delta_pos+exp(preds[pos]))**2)
        hess[neg] <- (weight*delta_neg*exp(preds[neg])) / ((delta_neg+exp(preds[neg]))**2)
        hess <- pmax(hess, hess_min)
        
        list(grad = grad, hess = hess)
      }
    }
    logistic_obj <- make_dw_logistic_obj(vol, F_prime)
    
  } 

  
  # --- Evaluasi Metrik (dengan format return yang disesuaikan untuk XGBoost) ---
  make_logistic_metric <- function(vol_vec) {
    .metric <- function(preds, dtrain) {
      labels <- getinfo(dtrain, "label")
      event <- which(labels > 0)
      dummy <- which(labels < 0)
      err <- numeric(length(labels))
      if (length(event) > 0) err[event] <- -log(exp(preds[event])/(vol_vec[event]+exp(preds[event])))
      if (length(dummy) > 0) err[dummy] <- -log(exp(vol_vec[dummy])/(exp(preds[dummy])+vol_vec[dummy]))
      
      # PERBEDAAN KUNCI: Format return value untuk feval di XGBoost
      return(list(name = "Error",value = sum(err) / 1e6, higher_better = FALSE))
    }
    return(.metric)
  }
  logistic_eval <- make_logistic_metric(vol)
  
  # 2. --- Penetapan Parameter & Pemanggilan xgb.train ---
  params$objective <- logistic_obj
  
  # PERBEDAAN KUNCI: Memanggil xgb.train dan menggunakan argumen 'feval'
  model <- xgb.train(
    params = params,
    feval = logistic_eval,
    maximize = FALSE,
    ...
  )
  
  return(model)
}
