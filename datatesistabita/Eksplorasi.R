library(spatstat)
library(xgboost)
library(RandomFields)
library(pracma)
library(caret)
library(openxlsx)

depth = load("C:/Users/Akmal/Documents/ITS/TA/Dataset/depth_tesis.Rda")
dip = load("C:/Users/Akmal/Documents/ITS/TA/Dataset/dip_tesis.Rda")
gempapp = load("C:/Users/Akmal/Documents/ITS/TA/Dataset/gempa_pp_tesis.Rda")
gempaqs = load("C:/Users/Akmal/Documents/ITS/TA/Dataset/gempa_qs_tesis.Rda")
ses = load("C:/Users/Akmal/Documents/ITS/TA/Dataset/ses_distfun.Rda")
sesar = load("C:/Users/Akmal/Documents/ITS/TA/Dataset/sesar_tesis.Rda")
strike = load("C:/Users/Akmal/Documents/ITS/TA/Dataset/strike_tesis.Rda")
subdist = load("C:/Users/Akmal/Documents/ITS/TA/Dataset/sub_distfun.Rda")
sub = load("C:/Users/Akmal/Documents/ITS/TA/Dataset/sub_tesis.Rda")
volcano = load("C:/Users/Akmal/Documents/ITS/TA/Dataset/volcano_tesis.Rda")
vuldist = load("C:/Users/Akmal/Documents/ITS/TA/Dataset/vul_distfun.Rda")

Poisson_objective <- function(preds, dtrain) {
  labels <- getinfo(dtrain, "label")
  weights <- getinfo(dtrain, "weight")
  preds1 <- preds[labels>0]
  preds2 <- preds[labels<0]
  weights <- weights/sum(weights[labels<0]*vol)
  grad1 <- rep(-1, length(preds1)) 
  grad2 <- exp(preds2)*vol
  grad <- c(grad1, grad2) 
  hess1 <- rep(0, length(preds1)) 
  hess2 <- exp(preds2)*vol
  hess <- c(hess1, hess2) 
  grad <- weights*grad
  hess <- weights*hess
  return(list(grad = grad, hess = hess))
}

DWPoisson_objective <- function(preds, dtrain) {
  labels <- getinfo(dtrain, "label")
  preds1 <- preds[labels>0]
  preds2 <- preds[labels<0]
  weights <- 1/(1+exp(preds)*F_prime)
  weights <- weights/sum(weights[labels<0]*vol)
  grad1 <- rep(-1, length(preds1)) 
  grad2 <- exp(preds2)*vol
  grad <- c(grad1, grad2) 
  hess1 <- rep(0, length(preds1)) 
  hess2 <- exp(preds2)*vol
  hess <- c(hess1, hess2)  
  grad <- weights*grad
  hess <- weights*hess
  return(list(grad = grad, hess = hess))
}

Poisson_error <- function(preds, dtrain) {
  labels <- getinfo(dtrain, "label")
  preds1 <- preds[labels>0]
  preds2 <- preds[labels<0]
  err1 <- -preds1
  err2 <- exp(preds2)*vol
  err <- c(err1, err2)
  return(list(metric = "Error", value = sum(err/10e6)))
}

data = data.frame()

for (i in 1:bei$n){
  x_index = floor((bei$x[i]+2.5)/5)+1
  y_index = floor((bei$y[i]+2.5)/5)+1
  convi = bci.covars$convi[y_index, x_index]
  Cu = bci.covars$Cu[y_index, x_index]
  dem = bci.covars$dem[y_index, x_index]
  devmean = bci.covars$devmean[y_index, x_index]
  difmean = bci.covars$difmean[y_index, x_index]
  grad = bci.covars$grad[y_index, x_index]
  elev = bei.extra$elev[y_index, x_index]
  K = bci.covars$K[y_index, x_index]
  mrvbf = bci.covars$mrvbf[y_index, x_index]
  Nmin = bci.covars$Nmin[y_index, x_index]
  P = bci.covars$P[y_index, x_index]
  pH = bci.covars$pH[y_index, x_index]
  solar = bci.covars$solar[y_index, x_index]
  twi = bci.covars$twi[y_index, x_index]
  data = rbind(data, c(bei$x[i], bei$y[i], convi, Cu, dem, devmean, difmean, grad, elev, K, mrvbf, Nmin, P, pH, solar, twi))
}
colnames(data) <- c('x','y','convi','Cu','dem','devmean','difmean','grad','elev','K','mrvbf','Nmin','P','pH','solar','twi')

data_pp <- ppp(bei$x, bei$y, window = bei$window)
plot(data_pp, main='bei data')

dummys = data.frame()
for (i in 1:201){
  for (j in 1:101){
    convi = bci.covars$convi[j, i]
    Cu = bci.covars$Cu[j, i]
    dem = bci.covars$dem[j, i]
    devmean = bci.covars$devmean[j, i]
    difmean = bci.covars$difmean[j, i]
    grad = bci.covars$grad[j, i]
    elev = bei.extra$elev[j, i]
    K = bci.covars$K[j, i]
    mrvbf = bci.covars$mrvbf[j, i]
    Nmin = bci.covars$Nmin[j, i]
    P = bci.covars$P[j, i]
    pH = bci.covars$pH[j, i]
    solar = bci.covars$solar[j, i]
    twi = bci.covars$twi[j, i]
    dummys <- rbind(dummys, c((i-1)*5, (j-1)*5, convi, Cu, dem, devmean, difmean, grad, elev, K, mrvbf, Nmin, P, pH, solar, twi))
  }
}
colnames(dummys) <- c('x','y','convi','Cu','dem','devmean','difmean','grad','elev','K','mrvbf','Nmin','P','pH','solar','twi')

#scale to unit window (2*1)
data$x = data$x/500
data$y = data$y/500
dummys$x = dummys$x/500
dummys$y = dummys$y/500

data$label <- 1
dummys$label <- -1
data_train <- rbind(data, dummys)
data_train$ids <- seq_len(nrow(data_train))
vol=1/10000

#original model
original_weights <- rep(1,dim(data_train)[1])
dtrain = xgb.DMatrix(data.matrix(data_train[,c('Cu','grad','elev','Nmin','P','pH','solar','twi')]), label=data.matrix(data_train$label), weight=original_weights)
score1 <- 10e6
for (eta in c(0.05)){
  for (lambda in c(10,30,50)){
    err_vector <- rep(0,2000)
    for (i in 1:3){
      print(i)
      set.seed(i)
      cv <- xgb.cv(data = dtrain, nrounds = 2000, obj = Poisson_objective, feval = Poisson_error, nfold = 2,
                   verbose = 0, num_parallel_tree = 3, subsample = 0.5, colsample_bynode=1/3, min_child_weight = 0.5, eta = eta, lambda = 0, alpha=lambda)
      err_vector <- err_vector + cv$evaluation_log$test_Error_mean
    }
    if (min(err_vector) < score1){
      eta1 = eta
      lambda1 = lambda
      best_round1 = which.min(err_vector)
      score1 <- min(err_vector)
    }
  }
}
model1 <- xgboost(data = dtrain, nrounds = best_round1, objective = Poisson_objective, eval_metric = Poisson_error,
                  verbose = 1, num_parallel_tree = 10, subsample = 0.5, colsample_bynode=1/3, min_child_weight = 0.5, eta = eta1, lambda=0, alpha=lambda1)
test_data = dummys[,c('Cu','grad','elev','Nmin','P','pH','solar','twi')]
pred <- predict(model1, data.matrix(test_data))
pred_model1 <- exp(pred)
pred_model1 <- Reshape(pred_model1, 101, 201)/250000 #if not divided by 250000, it is the intensity over the window 2*1
plot.im(as.im(pred_model1), main='bei; estimated intensities')

#importance1 <- xgb.importance(model = model1)
#xgb.plot.importance(importance_matrix = importance1)

preds_event = exp(predict(model1, data.matrix(data[,c('Cu','grad','elev','Nmin','P','pH','solar','twi')])))
pp = ppp(bei$x/500, bei$y/500, window = as.owin(c(c(0,2),c(0,1))))
sigma = 0.16
K_prime <- Kinhom(pp, preds_event, r=seq(0,sigma,0.01), correction = 'translation')
F_prime <- K_prime$trans[7] - K_prime$theo[7]


dtrain = xgb.DMatrix(data.matrix(data_train[,c('Cu','grad','elev','Nmin','P','pH','solar','twi')]), label=data.matrix(data_train$label))
score8 <- 10e6
for (eta in c(0.05)){
  for (lambda in c(10,30,50)){
    err_vector <- rep(0,2000)
    for (i in 1:3){
      print(i)
      set.seed(i)
      cv <- xgb.cv(data = dtrain, nrounds = 2000, obj = DWPoisson_objective, feval = Poisson_error, nfold = 2,
                   verbose = 0, num_parallel_tree = 3, subsample = 0.5, colsample_bynode=1/3, min_child_weight = 0.5, eta = eta, lambda = 0, alpha=lambda)
      err_vector <- err_vector + cv$evaluation_log$test_Error_mean
    }
    if (min(err_vector) < score8){
      eta8 = eta
      lambda8 = lambda
      best_round8 = which.min(err_vector)
      score8 <- min(err_vector)
    }
  }
}
model8 <- xgboost(data = dtrain, nrounds = best_round8, objective = DWPoisson_objective, eval_metric = Poisson_error,
                  verbose = 1, num_parallel_tree = 10, subsample = 0.5, colsample_bynode=1/3, min_child_weight = 0.5, eta = eta8, lambda=0, alpha=lambda8)
test_data = dummys[,c('Cu','grad','elev','Nmin','P','pH','solar','twi')]
pred <- predict(model8, data.matrix(test_data))
pred_model8 <- exp(pred)
pred_model8 <- Reshape(pred_model8, 101, 201)/250000
plot.im(as.im(pred_model8), main='bei data; estimated intensities')


importance <- xgb.importance(feature_names = c('Cu','Grad','Elev','Nmin','P','pH','Solar','Twi'), model = model8)
ggplot(importance, aes(x = Feature, y = Gain)) +
  geom_bar(stat = "identity", fill = "blue") +
  xlab("Covariates") +
  ylab("Importance") +
  theme(panel.background = element_blank(),  # Remove panel background
        plot.background = element_blank(),   # Remove overall plot background
        panel.grid = element_blank(), 
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle = 60, hjust = 1), )


wb <- createWorkbook()
addWorksheet(wb, "result")
result <- pred_model8
writeData(wb, "result", result, startRow = 1, startCol = 1)
saveWorkbook(wb, file = paste("result_bei.xlsx", sep=""))



# 4-fold cross validation
set.seed(1)
data$ids <- seq_len(nrow(data))
flds <- createFolds(data$ids, k = 4)
ll1s <- c()
ll8s <- c()

for (fold in 1:4){
  print(fold)
  data_infold <- data[-which(data$ids %in% flds[[fold]]),]
  test_infold <- data[which(data$ids %in% flds[[fold]]),]
  data_train_infold <- rbind(data_infold[,1:17], dummys)
  
  original_weights <- rep(1,dim(data_train_infold)[1])
  dtrain_infold = xgb.DMatrix(data.matrix(data_train_infold[,c('Cu','grad','elev','Nmin','P','pH','solar','twi')]), label=data.matrix(data_train_infold$label), weight=original_weights)
  score1 <- 10e6
  for (eta in c(0.05)){
    for (lambda in c(10,30,50)){
      err_vector <- rep(0,3000)
      for (i in 1:3){
        print(i)
        set.seed(i)
        cv <- xgb.cv(data = dtrain_infold, nrounds = 3000, obj = Poisson_objective, feval = Poisson_error, nfold = 2,
                     verbose = 0, num_parallel_tree = 3, subsample = 0.5, colsample_bynode=1/3, min_child_weight = 0.5, eta = eta, lambda = 0, alpha=lambda)
        err_vector <- err_vector + cv$evaluation_log$test_Error_mean
      }
      if (min(err_vector) < score1){
        eta1 = eta
        lambda1 = lambda 
        best_round1 = which.min(err_vector)
        score1 <- min(err_vector)
      }
    }
  }
  model1_infold <- xgboost(data = dtrain_infold, nrounds = best_round1, objective = Poisson_objective, eval_metric = Poisson_error,
                           verbose = 1, num_parallel_tree = 10, subsample = 0.5, colsample_bynode=1/3, min_child_weight = 0.5, eta = eta1, lambda=0, alpha=lambda1)
  test_data = dummys[,c('Cu','grad','elev','Nmin','P','pH','solar','twi')]
  pred <- predict(model1_infold, data.matrix(test_data))
  pred <- exp(pred)
  pred <- Reshape(pred, 101, 201)/250000
  plot.im(as.im(pred), main='bei; estimated intensities')
  
  ll1 = sum(log(exp(predict(model1_infold, data.matrix(test_infold[,c('Cu','grad','elev','Nmin','P','pH','solar','twi')])))/250000/3)) -
    sum(exp(predict(model1_infold, data.matrix(dummys[,c('Cu','grad','elev','Nmin','P','pH','solar','twi')])))/250000/3*25)
  ll1s = c(ll1s, ll1)
  
  preds_infold = exp(predict(model1_infold, data.matrix(data_infold[,c('Cu','grad','elev','Nmin','P','pH','solar','twi')])))
  pp = ppp(data_infold$x, data_infold$y, window = as.owin(c(c(0,2),c(0,1))))
  sigma = 0.16
  K_prime <- Kinhom(pp, preds_infold, r=seq(0,sigma,0.01), correction = 'translation')
  F_prime <- K_prime$trans[7] - K_prime$theo[7]
  
  dtrain_infold = xgb.DMatrix(data.matrix(data_train_infold[,c('Cu','grad','elev','Nmin','P','pH','solar','twi')]), label=data.matrix(data_train_infold$label))
  score8 <- 10e6
  for (eta in c(0.05)){
    for (lambda in c(10,30,50)){
      err_vector <- rep(0,3000)
      for (i in 1:3){
        print(i)
        set.seed(i)
        cv <- xgb.cv(data = dtrain_infold, nrounds = 3000, obj = DWPoisson_objective, feval = Poisson_error, nfold = 2,
                     verbose = 0, num_parallel_tree = 3, subsample = 0.5, colsample_bynode=1/3, min_child_weight = 0.5, eta = eta, lambda = 0, alpha=lambda)
        err_vector <- err_vector + cv$evaluation_log$test_Error_mean
      }
      if (min(err_vector) < score8){
        eta8 = eta
        lambda8 = lambda
        best_round8 = which.min(err_vector)
        score8 <- min(err_vector)
      }
    }
  }
  model8_infold <- xgboost(data = dtrain_infold, nrounds = best_round8, objective = DWPoisson_objective, eval_metric = Poisson_error,
                           verbose = 1, num_parallel_tree = 10, subsample = 0.5, colsample_bynode=1/3, min_child_weight = 0.5, eta = eta8, lambda=0, alpha=lambda8)
  test_data = dummys[,c('Cu','grad','elev','Nmin','P','pH','solar','twi')]
  pred <- predict(model8_infold, data.matrix(test_data))
  pred <- exp(pred)
  pred <- Reshape(pred, 101, 201)/250000
  plot.im(as.im(pred), main='bei data; estimated intensities')
  
  ll8 = sum(log(exp(predict(model8_infold, data.matrix(test_infold[,c('Cu','grad','elev','Nmin','P','pH','solar','twi')])))/250000/3)) -
    sum(exp(predict(model8_infold, data.matrix(dummys[,c('Cu','grad','elev','Nmin','P','pH','solar','twi')])))/250000/3*25)
  ll8s = c(ll8s, ll8)
}





