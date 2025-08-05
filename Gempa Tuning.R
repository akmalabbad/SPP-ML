library(spatstat)
library(xgboost)
#library(RandomFields)
library(pracma)
library(caret)
library(openxlsx)
library(raster)
library(spatstat.geom)
library(RColorBrewer)
library(akima)
library(viridis)
library(lightgbm)
library(tictoc)
library(spatstat.geom)
library(dplyr)

setwd('C:/Users/Akmal/Documents/ITS/TA')

# Load your data
depth = load("datatesistabita/depth_tesis.Rda")
dip = load("datatesistabita/dip_tesis.Rda")
gempapp = load("datatesistabita/gempa_pp_tesis.Rda")
gempaqs = load("datatesistabita/gempa_qs_tesis.Rda")
load("datatesistabita/ses_distfun.rda")
sesar = load("datatesistabita/sesar_tesis.Rda")
strike = load("datatesistabita/strike_tesis.Rda")
load("datatesistabita/sub_distfun.rda")
subdata = load("datatesistabita/sub_tesis.Rda")
volcano = load("datatesistabita/volcano_tesis.Rda")
load("datatesistabita/vul_distfun.rda")

megathrust_im <- distmap(megathrustSM)
fault_im <- distmap(faultSM)

x_data <- gempaSM1$data$x
y_data <- gempaSM1$data$y
x_dummy <- gempaSM1$dummy$x
y_dummy <- gempaSM1$dummy$y

# Coordinates as data.frame
coordinates <- data.frame(x = x_data, y = y_data)
dummycoor <- data.frame(x = x_dummy, y = y_dummy)

# Convert im to raster
depth_raster <- raster(as.im(depthSM))
dip_raster <- raster(as.im(dipSM))
strike_raster <- raster(as.im(strikeSM))

# Extract covariate values
depth_values <- raster::extract(depth_raster, coordinates)
dip_values <- raster::extract(dip_raster, coordinates)
strike_values <- raster::extract(strike_raster, coordinates)

depth_values_dummy <- raster::extract(depth_raster, dummycoor)
dip_values_dummy <- raster::extract(dip_raster, dummycoor)
strike_values_dummy <- raster::extract(strike_raster, dummycoor)

# Compute sesar, volcano, sub values
sesar_values <- ses(coordinates$x, coordinates$y)
volcano_values <- vul(coordinates$x, coordinates$y)
sub_values <- sub(coordinates$x, coordinates$y)

sesar_values_dummy <- ses(dummycoor$x, dummycoor$y)
volcano_values_dummy <- vul(dummycoor$x, dummycoor$y)
sub_values_dummy <- sub(dummycoor$x, dummycoor$y)


earthquakes <- read.csv('query.csv', sep = ';');head(earthquakes)
win_left  <- spatstat.geom::shift(gempaSM$window, vec = c(-11.3, 0))
pp_gempa <- ppp(x = earthquakes$longitude,
                y = earthquakes$latitude,
                window = win_left,
                marks = earthquakes$mag)
plot(pp_gempa)

# pilih hanya titik yang berada di dalam window
ok          <- inside.owin(pp_gempa$x, pp_gempa$y, pp_gempa$window)
pp_gempa_in <- pp_gempa[ok]          # subset method untuk ppp

plot(pp_gempa_in)
gempaSM <- pp_gempa_in
# make the plot window square so 1° in x equals 1° in y
plot(gempaSM$x, gempaSM$y,
     xlab = "long", ylab = "lat",
     col  = adjustcolor("black", alpha.f = 0.4),
     asp  = TRUE, axes = FALSE,
     cex  = 1.5 * (gempaSM$marks) / max(gempaSM$marks))

maps::map("world", add = TRUE, col = "grey10",
          xlim = gempaSM$window$xrange,
          ylim = gempaSM$window$yrange)

axis(1); axis(2)
polygon(gempaSM$window$bdry[[1]],
        border = 2)


#------------MAKE THE DUMMY POINTS----------------------
qd <- spatstat.geom::quadscheme(unmark(pp_gempa_in))
xg <- spatstat.geom::x.quad(qd)
yg <- spatstat.geom::y.quad(qd)
wg <- spatstat.geom::w.quad(qd)
plot(qd, cex = 0.6,
     main = "Grid Discretization")
sum(qd$w)
length(qd$w)
area.owin(gempaSM$window)
plot(qd$data)
plot(qd$dummy)
#------------EXTRACT THE COORDINATES--------------------
x_data <- pp_gempa_in$x
y_data <- pp_gempa_in$y
x_dummy <- xg
y_dummy <- yg

# Coordinates as data.frame
coordinates <- data.frame(x = x_data, y = y_data)
dummycoor <- data.frame(x = x_dummy, y = y_dummy)

# returns TRUE for every observation that has appeared before
dup_flag <- duplicated(cbind(coordinates$x, coordinates$y))

which(dup_flag)          # row indices of the second (and later) copies
sum(dup_flag)            # how many duplicated records

#depth_raster_left <- raster::shift(depth_raster, dx = -10, dy = 0);depth_raster_left
#pp_gempa_in
# Convert im to raster
depth_raster <- raster::shift(raster(as.im(depthSM)), dx = -11.3, dy = 0)
dip_raster <- raster::shift(raster(as.im(dipSM)), dx = -11.3, dy = 0)
strike_raster <- raster::shift(raster(as.im(strikeSM)), dx = -11.3, dy = 0)

# Extract covariate values
depth_values <- raster::extract(depth_raster, coordinates)
dip_values <- raster::extract(dip_raster, coordinates)
strike_values <- raster::extract(strike_raster, coordinates)

depth_values_dummy <- raster::extract(depth_raster, dummycoor)
dip_values_dummy <- raster::extract(dip_raster, dummycoor)
strike_values_dummy <- raster::extract(strike_raster, dummycoor)

# Compute sesar, volcano, sub values
sesar_values <- ses(coordinates$x, coordinates$y)
volcano_values <- vul(coordinates$x, coordinates$y)
sub_values <- sub(coordinates$x, coordinates$y)

sesar_values_dummy <- ses(dummycoor$x, dummycoor$y)
volcano_values_dummy <- vul(dummycoor$x, dummycoor$y)
sub_values_dummy <- sub(dummycoor$x, dummycoor$y)

#-----------------------DATA EXTRACTION----------------------
# Combine true data with covariates and label
truedata <- cbind(
  coordinates,
  depth_values, dip_values, sesar_values,
  strike_values, sub_values, volcano_values,
  label = 1,
  vol = 0.01
)
colnames(truedata) <- c("x", "y",
                        "depth", "dip", "sesar",
                        "strike", "megathrust", "volcano",
                        "label",
                        "vol")

# Combine dummy data with covariates and label
dummydata <- cbind(
  dummycoor,
  depth_values_dummy, dip_values_dummy, sesar_values_dummy,
  strike_values_dummy, sub_values_dummy, volcano_values_dummy,
  label = -1,
  vol = wg
)
colnames(dummydata) <- c("x", "y",
                         "depth", "dip", "sesar",
                         "strike", "megathrust", "volcano",
                         "label", "vol")
dummydata <- na.omit(dummydata)
datagempa <- rbind(truedata, dummydata); head(datagempa)
vol <- wg; length(vol)

#--------------------EDA--------------------------------
## 1. Fault Distance
dummy_ppp <- ppp(x = dummydata$x, 
                 y = dummydata$y, 
                 window = win_left)

# Add prediction as marks
marks(dummy_ppp) <- dummydata$sesar
plot(dummy_ppp)
# Smoothing (Intensity Estimation) from prediction
intensity_map <- Smooth(dummy_ppp)

par(pty = "s")                                   # square plotting frame

plot(intensity_map,                              # first layer
     col  = viridis::viridis(100),               # any palette
     ribbon = TRUE,                              # colour-bar
     asp  = 1,
     main = "Fault Image")                                   # enforces 1:1

maps::map("world",
          add  = TRUE,
          col  = "grey10",
          lwd  = 1,
          xlim = intensity_map$xr,               # 94.534 … 111.18
          ylim = intensity_map$yr)               # −6.85 … 7.35
points(gempaSM$x, gempaSM$y,
       col  = adjustcolor("black", alpha.f = 0.4),
       asp  = TRUE, axes = FALSE,
       cex  = 1 * (gempaSM$marks + 0.1) / max(gempaSM$marks))
polygon(gempaSM$window$bdry[[1]], border = 2)    # study window, optional
axis(1); axis(2)

## 2. Subduction Distance
dummy_ppp <- ppp(x = dummydata$x, 
                 y = dummydata$y, 
                 window = win_left)

marks(dummy_ppp) <- dummydata$megathrust
intensity_map <- Smooth(dummy_ppp)

par(pty = "s")                                   # square plotting frame

plot(intensity_map,                              # first layer
     col  = viridis::viridis(100),               # any palette
     ribbon = TRUE,                              # colour-bar
     asp  = 1,
     main = "Subduction Image")                                   # enforces 1:1

maps::map("world",
          add  = TRUE,
          col  = "grey10",
          lwd  = 1,
          xlim = intensity_map$xr,               # 94.534 … 111.18
          ylim = intensity_map$yr)               # −6.85 … 7.35
points(gempaSM$x, gempaSM$y,
       col  = adjustcolor("black", alpha.f = 0.4),
       asp  = TRUE, axes = FALSE,
       cex  = 1 * (gempaSM$marks + 0.1) / max(gempaSM$marks))
polygon(gempaSM$window$bdry[[1]], border = 2)    # study window, optional
axis(1); axis(2)

## 3. Volcanoes Distance
dummy_ppp <- ppp(x = dummydata$x, 
                 y = dummydata$y, 
                 window = win_left)


marks(dummy_ppp) <- dummydata$volcano
intensity_map <- Smooth(dummy_ppp)

par(pty = "s")                                   # square plotting frame

plot(intensity_map,                              # first layer
     col  = viridis::viridis(100),               # any palette
     ribbon = TRUE,                              # colour-bar
     asp  = 1,
     main = "Volcano Image")                                   # enforces 1:1

maps::map("world",
          add  = TRUE,
          col  = "grey10",
          lwd  = 1,
          xlim = intensity_map$xr,               # 94.534 … 111.18
          ylim = intensity_map$yr)               # −6.85 … 7.35
points(gempaSM$x, gempaSM$y,
       col  = adjustcolor("black", alpha.f = 0.4),
       asp  = TRUE, axes = FALSE,
       cex  = 1 * (gempaSM$marks + 0.1) / max(gempaSM$marks))
polygon(gempaSM$window$bdry[[1]], border = 2)    # study window, optional
axis(1); axis(2)

## 4. Dip
dummy_ppp <- ppp(x = dummydata$x, 
                 y = dummydata$y, 
                 window = win_left)

sum(is.na(dummydata$dip))
marks(dummy_ppp) <- dummydata$dip
intensity_map <- Smooth(dummy_ppp)

par(pty = "s")                                   # square plotting frame

plot(intensity_map,                              # first layer
     col  = viridis::viridis(100),               # any palette
     ribbon = TRUE,                              # colour-bar
     asp  = 1,
     main = "Dip Image")                                   # enforces 1:1

maps::map("world",
          add  = TRUE,
          col  = "grey10",
          lwd  = 1,
          xlim = intensity_map$xr,               # 94.534 … 111.18
          ylim = intensity_map$yr)               # −6.85 … 7.35
points(gempaSM$x, gempaSM$y,
       col  = adjustcolor("black", alpha.f = 0.4),
       asp  = TRUE, axes = FALSE,
       cex  = 1 * (gempaSM$marks + 0.1) / max(gempaSM$marks))
polygon(gempaSM$window$bdry[[1]], border = 2)    # study window, optional
axis(1); axis(2)

## 5. Depth
dummy_ppp <- ppp(x = dummydata$x, 
                 y = dummydata$y, 
                 window = win_left)


marks(dummy_ppp) <- dummydata$depth
intensity_map <- Smooth(dummy_ppp)

par(pty = "s")                                   # square plotting frame

plot(intensity_map,                              # first layer
     col  = viridis::viridis(100),               # any palette
     ribbon = TRUE,                              # colour-bar
     asp  = 1,
     main = "Depth Image")                                   # enforces 1:1

maps::map("world",
          add  = TRUE,
          col  = "grey10",
          lwd  = 1,
          xlim = intensity_map$xr,               # 94.534 … 111.18
          ylim = intensity_map$yr)               # −6.85 … 7.35
points(gempaSM$x, gempaSM$y,
       col  = adjustcolor("black", alpha.f = 0.4),
       asp  = TRUE, axes = FALSE,
       cex  = 1 * (gempaSM$marks + 0.1) / max(gempaSM$marks))
polygon(gempaSM$window$bdry[[1]], border = 2)    # study window, optional
axis(1); axis(2)

## 6. Strike
dummy_ppp <- ppp(x = dummydata$x, 
                 y = dummydata$y, 
                 window = win_left)


marks(dummy_ppp) <- dummydata$strike
intensity_map <- Smooth(dummy_ppp)

par(pty = "s")                                   # square plotting frame

plot(intensity_map,                              # first layer
     col  = viridis::viridis(100),               # any palette
     ribbon = TRUE,                              # colour-bar
     asp  = 1,
     main = "Strike Image")                                   # enforces 1:1

maps::map("world",
          add  = TRUE,
          col  = "grey10",
          lwd  = 1,
          xlim = intensity_map$xr,               # 94.534 … 111.18
          ylim = intensity_map$yr)               # −6.85 … 7.35
points(gempaSM$x, gempaSM$y,
       col  = adjustcolor("black", alpha.f = 0.4),
       asp  = TRUE, axes = FALSE,
       cex  = 1 * (gempaSM$marks + 0.1) / max(gempaSM$marks))
polygon(gempaSM$window$bdry[[1]], border = 2)    # study window, optional
axis(1); axis(2)

#write.csv(datagempa,'datagempa.csv')

datagempa <- datagempa %>% 
  slice(-head(which(datagempa$label == -1), sum(datagempa$label == 1)))
head(datagempa)
vol <- datagempa$vol

dummydata <- dummydata[datagempa$label==-1,]; head(dummydata)
dummy_ppp <- ppp(x = dummydata$x, 
                 y = dummydata$y, 
                 window = win_left)
plot(dummy_ppp)

features <- c("depth", "dip", "sesar", "strike", "megathrust", "volcano")
train_set <- lgb.Dataset(data = data.matrix(datagempa[, features]), label = datagempa$label)

# =============================================================
# BAGIAN 1: DEFINISI FUNGSI WRAPPER "lgbpp" (Tidak berubah)
# =============================================================
source('lgbpp v2.R')

# =============================================================
# BAGIAN 2: PROSES TUNING DENGAN OUTPUT DETAIL
# =============================================================

# --- Parameter Dasar ---
base_params <- list(
  # GOSS
  boosting_type = "goss",
  top_rate = 0.2,
  other_rate = 0.2,
  
  # EFB
  is_enable_bundle = TRUE,
  max_conflict_rate = 0.1,
  
  # Histogram
  max_bin = 255,
  bin_construct_sample_cnt = 200000,
  min_data_in_bin = 3,
  
  # Other
  min_data_in_leaf = 20,
  min_sum_hessian_in_leaf = 1e-3, 
  feature_fraction = 1/3,
  deterministic = TRUE,
  verbose = -1,
  num_threads = parallel::detectCores(),
  learning_rate = 0.021161616161616163,
  lambda_l1 = 0.7071067811865476,
  lambda_l2 = 0.005524271728019903,
  num_leaves = 3
)

model_final <- lgbpp(
  vol = vol,
  params = base_params,
  data = train_set,
  nrounds = 5000,
  valids = list(train = train_set),
  early_stopping_rounds = 50
)

# 1. Ekstrak data evaluasi dari objek model
eval_results <- c()
for(i in 1:model_final$best_iter){
  eval_results <- rbind(eval_results,model_final$record_evals$train$poisson_err$eval[[i]])
}

# 2. Buat data frame yang rapi untuk plotting
plot_df <- data.frame(
  iterasi = 1:length(eval_results),
  poisson_error = eval_results
)
# 3. Buat visualisasi menggunakan ggplot2
ggplot(plot_df, aes(x = iterasi, y = poisson_error)) +
  geom_line(color = "dodgerblue", size = 1) +  # Membuat plot garis
  
  # Tambahkan garis vertikal untuk menandai iterasi terbaik
  geom_vline(
    xintercept = model_final$best_iteration, 
    color = "red",
    linetype = "dashed",
    linewidth = 1
  ) +
  
  # Beri judul dan label yang jelas
  labs(
    title = "Perkembangan Poisson Error Selama Training",
    subtitle = "Metrik dievaluasi pada data training",
    x = "Jumlah Pohon (Iterasi)",
    y = "Poisson Error"
  ) +
  theme_minimal()

dummy_mask <- datagempa$label < 0
event_mask <- datagempa$label > 0
vol_dummy <- vol[dummy_mask]

y_pred_dummy <- predict(model_final, data.matrix(datagempa[dummy_mask, features]), raw = TRUE)
y_pred_event <- predict(model_final, data.matrix(datagempa[event_mask, features]), raw = TRUE)

cat("Number of true events : ", length(y_pred_event))
cat("Predicted number of events : ",sum(exp(y_pred_dummy)*vol_dummy))

# Add prediction as marks
marks(dummy_ppp) <- y_pred_dummy
plot(dummy_ppp)
# Smoothing (Intensity Estimation) from prediction
intensity_map <- Smooth(dummy_ppp)

par(pty = "s")                                   # square plotting frame

plot(intensity_map,                              # first layer
     col  = viridis::viridis(100, option = 'plasma'),               # any palette
     ribbon = TRUE,                              # colour-bar
     asp  = 1,                                   # enforces 1:1
     main = "Predicted Log-Intensity (LightGBMPP)")

maps::map("world",
          add  = TRUE,
          col  = "grey10",
          lwd  = 1,
          xlim = intensity_map$xr,               # 94.534 … 111.18
          ylim = intensity_map$yr)               # −6.85 … 7.35
points(gempaSM$x, gempaSM$y,
     col  = adjustcolor("black", alpha.f = 0.4),
     asp  = TRUE, axes = FALSE,
     cex  = 1 * (gempaSM$marks + 0.1) / max(gempaSM$marks))
polygon(gempaSM$window$bdry[[1]], border = 2)    # study window, optional
axis(1); axis(2)

# =============================================================
# BAGIAN 1: DEFINISI FUNGSI WRAPPER "lgbpp" (Tidak berubah)
# =============================================================
source('xgbpp v2.R')

# =============================================================
# BAGIAN 2: PROSES TUNING DENGAN OUTPUT DETAIL
# =============================================================
dtrain <- xgb.DMatrix(
  data = data.matrix(datagempa[, features]),
  label = datagempa$label
)

# --- Parameter Dasar ---
base_params_xgb <- list(
  booster = "gbtree",
  subsample = 0.8,
  colsample_bytree = 1/3, # Padanan dari feature_fraction
  nthread = parallel::detectCores(),
  tree_method = 'hist',
  eta = 0.021161616161616163,
  alpha = 0.7071067811865476,
  lambda = 0.005524271728019903,
  max_depth = 2
)

# 2. Latih model final dengan semua data training menggunakan xgbpp
model_final_xgb <- xgbpp(
  vol = vol,
  params = base_params_xgb,
  data = dtrain,
  nrounds = 5000,
  watchlist = list(train = dtrain), # Argumen XGBoost adalah 'watchlist'
  early_stopping_rounds = 50,
  verbose = 0 # Menonaktifkan log per ronde
)

# 1. Ekstrak data evaluasi dari objek model
eval_results_xgb <- c()
for(i in 1:model_final_xgb$best_iteration){
  eval_results_xgb <- rbind(eval_results_xgb,model_final_xgb$evaluation_log$train[i])
}

# 2. Buat data frame yang rapi untuk plotting
plot_df_xgb <- data.frame(
  iterasi = 1:length(eval_results_xgb),
  poisson_error = eval_results_xgb
)
# 3. Buat visualisasi menggunakan ggplot2
ggplot(plot_df_xgb, aes(x = iterasi, y = poisson_error)) +
  geom_line(color = "dodgerblue", size = 1) +  # Membuat plot garis
  
  # Tambahkan garis vertikal untuk menandai iterasi terbaik
  geom_vline(
    xintercept = model_final$best_iteration, 
    color = "red",
    linetype = "dashed",
    linewidth = 1
  ) +
  
  # Beri judul dan label yang jelas
  labs(
    title = "Perkembangan Poisson Error Selama Training",
    subtitle = "Metrik dievaluasi pada data training",
    x = "Jumlah Pohon (Iterasi)",
    y = "Poisson Error"
  ) +
  theme_minimal()

dummy_mask <- datagempa$label < 0
event_mask <- datagempa$label > 0
vol_dummy <- vol[dummy_mask]

y_pred_dummy_xgb <- predict(model_final_xgb, data.matrix(datagempa[dummy_mask, features]), raw = TRUE)
y_pred_event_xgb <- predict(model_final_xgb, data.matrix(datagempa[event_mask, features]), raw = TRUE)

cat("Number of true events : ", length(y_pred_event))
cat("Predicted number of events : ",sum(exp(y_pred_dummy_xgb)*vol_dummy))

# Add prediction as marks
marks(dummy_ppp) <- y_pred_dummy_xgb
plot(dummy_ppp)
# Smoothing (Intensity Estimation) from prediction
intensity_map <- Smooth(dummy_ppp)

par(pty = "s")                                   # square plotting frame

plot(intensity_map,                              # first layer
     col  = viridis::viridis(100, option = 'plasma'),               # any palette
     ribbon = TRUE,                              # colour-bar
     asp  = 1,                                   # enforces 1:1
     main = "Predicted Log-Intensity (XGBoostPP)")

maps::map("world",
          add  = TRUE,
          col  = "grey10",
          lwd  = 1,
          xlim = intensity_map$xr,               # 94.534 … 111.18
          ylim = intensity_map$yr)               # −6.85 … 7.35
points(gempaSM$x, gempaSM$y,
       col  = adjustcolor("black", alpha.f = 0.4),
       asp  = TRUE, axes = FALSE,
       cex  = 1 * (gempaSM$marks + 0.1) / max(gempaSM$marks))
polygon(gempaSM$window$bdry[[1]], border = 2)    # study window, optional
axis(1); axis(2)
