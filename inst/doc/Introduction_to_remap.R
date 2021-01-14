## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(magrittr) # For pipe %>% functionality
library(tibble)   # For light data wrangling
library(dplyr)    # For light data wrangling
library(ggplot2)  # For plots
library(maps)     # For a polygon of the state of Utah
library(sf)       # For spatial data manipulation
library(mgcv)     # For GAM modeling

library(remap)
data(utsnow)
data(utws)

# Reset CRS in case user has old version of GDAL 
sf::st_crs(utsnow) <- 4326
sf::st_crs(utws) <- 4326

## ----initial_map, fig.width = 6, fig.height = 6, fig.align='center'-----------
utstate <- maps::map("state", plot = FALSE, fill = TRUE) %>%
  sf::st_as_sf() %>%
  dplyr::filter(ID == "utah") %>%
  sf::st_transform(crs = 4326)

ggplot(utws, aes(fill = HUC2)) +
  geom_sf(alpha = 0.5) +
  geom_sf(data = utstate, fill = "NA", size = 1) +
  geom_sf(data = utsnow) +
  ggtitle("Modeling Data and Regions",
          "HUC2 regions are made up of smaller HUC4 regions.") +
  theme_void()

## ----utsnz--------------------------------------------------------------------
utsnz <- utsnow %>% dplyr::filter(WESD > 0)

## ----wesd_elev, fig.width = 5, fig.height = 2, fig.align='center'-------------
ggplot(utsnz, aes(x = ELEVATION, y = WESD)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_y_log10() +
  theme_minimal()

## ----lm-----------------------------------------------------------------------
lm_global <- lm(log(WESD) ~ ELEVATION, data = utsnz)

lm_global_mse <- mean((utsnz$WESD - exp(predict(lm_global, utsnz)))^2)
lm_global_mse

## ----wesd_elev2, fig.width = 5, fig.height = 7, fig.align='center'------------
ggplot(utsnz, aes(x = ELEVATION, y = WESD)) +
  facet_grid(HUC2 ~ .) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_y_log10() +
  theme_minimal()

## ----lm_huc2, message=FALSE---------------------------------------------------
t1 <- Sys.time()
lm_huc2 <- remap::remap(
  utsnz, utws, region_id = HUC2, 
  buffer = 20, min_n = 10,
  model_function = lm, 
  formula = log(WESD) ~ ELEVATION
)
  
lm_huc2_mse <- mean((utsnz$WESD - exp(predict(lm_huc2, utsnz, smooth = 10)))^2)
t2 <- Sys.time()

# mse
lm_huc2_mse

# runtime
round(difftime(t2, t1), 1)

## ----lm_huc2_models, message=FALSE--------------------------------------------
sapply(lm_huc2$models, function(x) x$coefficients)

## ----polygons, fig.width = 6, fig.height = 4, fig.align='center'--------------
utws_simp <- utws %>% sf::st_simplify(dTolerance = 0.05)

rbind(
  utws %>% dplyr::mutate(TYPE = "Original Watershed Polygons"),
  utws_simp %>% dplyr::mutate(TYPE = "Simplified Watershed Polygons")
) %>%
ggplot() +
  facet_grid(.~TYPE) +
  geom_sf() +
  theme_void()

## ----redist, message=FALSE----------------------------------------------------
huc2_dist_nz <- remap::redist(utsnz, utws_simp, region_id = HUC2)
head(huc2_dist_nz)

## ----lm_huc2_simp, message=FALSE----------------------------------------------
t1 <- Sys.time()
lm_huc2 <- remap(
  utsnz, utws_simp, region_id = HUC2, 
  buffer = 20, min_n = 10,
  model_function = lm, 
  formula = log(WESD) ~ ELEVATION,
  distances = huc2_dist_nz
)
  
lm_huc2_mse <- mean(
  (utsnz$WESD - 
     exp(predict(lm_huc2, utsnz, smooth = 10,
                 distances = huc2_dist_nz)))^2
)
t2 <- Sys.time()

# mse
lm_huc2_mse

# runtime
round(difftime(t2, t1), 1)

## ----gam----------------------------------------------------------------------
gam_limit <- function(data, fml) {
  g_model <- mgcv::gam(fml, data = data)
  upper_lim <- max(data$WESD)
  
  out <- list(g_model = g_model, upper_lim = upper_lim)
  class(out) <- "gam_limit"
  return(out)
}

predict.gam_limit <- function(object, newobs) {
  if (nrow(newobs) != 0) {
    preds <- predict(object$g_model, newobs)
    
    preds[preds < 0] <- 0
    preds[preds > object$upper_lim] <- object$upper_lim
    
    return(preds)
  }
  return(NULL)
}

## ----gam_global---------------------------------------------------------------
# Create vector for cross-validation
set.seed(42)
cv_k <- sample(1:10, nrow(utsnow), replace = TRUE)

# Initialize predictions
gam_global_preds <- rep(as.numeric(NA), nrow(utsnow))

# Formula for global GAM
global_fml <- WESD ~ s(ELEVATION, k = 5) + s(LATITUDE, LONGITUDE, bs = "sos", k = 50)

# Build and test models with 10 fold cross-validation
for (i in 1:10) {
  index <- cv_k == i
  gam_global <- gam_limit(utsnow[!index, ], fml = global_fml)
  gam_global_preds[index] <- predict(gam_global, utsnow[index, ])
}

# Calculate MSE
gam_global_mse <- mean((utsnow$WESD - gam_global_preds)^2)
gam_global_mse

## ----dist---------------------------------------------------------------------
huc2_dist <- remap::redist(utsnow, utws, region_id = HUC2)

## ----gam_huc2-----------------------------------------------------------------
# Initialize predictions
gam_huc2_preds <- rep(as.numeric(NA), nrow(utsnow))

# Formula for regional GAMs
gam_huc2_fml <- WESD ~ s(ELEVATION, k = 5) + s(LATITUDE, LONGITUDE, bs = "sos", k = 25)

# Build and test models with 10 fold cross-validation
for (i in 1:10) {
  index <- cv_k == i
  
  gam_huc2 <- remap::remap(
    utsnow[!index, ], utws, region_id = HUC2,
    model_function = gam_limit, 
    buffer = 20, min_n = 35,
    distances = huc2_dist[!index, ],
    fml = gam_huc2_fml
  )
  
  gam_huc2_preds[index] <- predict(
    gam_huc2, utsnow[index, ],
    smooth = 10, 
    distances = huc2_dist[index, ]
  )
}

# Calculate MSE
gam_huc2_mse <- mean((utsnow$WESD - gam_huc2_preds)^2)
gam_huc2_mse

## ----toy----------------------------------------------------------------------
# Make regions
toy_regions <- tibble::tribble(
  ~id, ~geometry,
  "a", sf::st_polygon(list(matrix(c(0, 0, 2, 0, 6, 3, 4, 10, 0, 10, 0, 0)*.1, ncol = 2, byrow = TRUE))),
  "b", sf::st_polygon(list(matrix(c(2, 0, 10, 0, 10, 4, 6, 3, 2, 0)*.1, ncol = 2, byrow = TRUE))),
  "c", sf::st_polygon(list(matrix(c(4, 10, 6, 3, 10, 4, 10, 10, 4, 10)*.1, ncol = 2, byrow = TRUE)))
) %>%
  sf::st_as_sf(crs = 4326)

# Manually make a toy remap model
make_toy <- function(x) {
  class(x) <- "toy_model"
  return(x)
}
remap_toy_model <- list(
  models = list("a" = make_toy("a"), 
                "b" = make_toy("b"), 
                "c" = make_toy("c")),
  regions = toy_regions,
  region_id = "id"
)
class(remap_toy_model) <- "remap"

# Make a prediction method for toy_model
predict.toy_model <- function(object, data) {
  x <- sf::st_coordinates(data)[, "X"]
  y <- sf::st_coordinates(data)[, "Y"]
  if (object == "a") {
    y - x
  } else if (object == "b") {
    x - y - 0.4
  } else {
    y - x + 0.3
  }
}

# Make a grid over the regions for predictions
grd <- sf::st_make_grid(toy_regions, cellsize = .01, what = "corners") %>%
  sf::st_sf()

## ----toy_regions, fig.width = 4, fig.height = 4, fig.align='center'-----------
ggplot2::ggplot(toy_regions, aes(fill = id)) +
    geom_sf(color = "black", size = 1) +
    ggtitle("Toy Regions") +
    theme_bw()

## ----grid---------------------------------------------------------------------
grd_pred <- grd %>%
  dplyr::mutate(SHARP = predict(remap_toy_model, grd, smooth = 0),
                SMOOTH = predict(remap_toy_model, grd, smooth = 30),
                LON = sf::st_coordinates(.)[, "X"],
                LAT = sf::st_coordinates(.)[, "Y"])

## ----sharp, fig.width = 5, fig.height = 4, fig.align='center'-----------------
ggplot(toy_regions) +
  geom_sf() +
  geom_tile(data = grd_pred, aes(x = LON, y = LAT, fill = SHARP)) +
  scale_fill_viridis_c(limits = c(-0.3, 1)) +
  geom_hline(yintercept = 0.8) +
  ggtitle("Sharp Predictions", "Black line corresponds to x-axis of the next plot.") +
  xlab("") + ylab("") +
  theme_bw()

## ----sharp08, fig.width = 4, fig.height = 2.5, fig.align='center'-------------
ggplot(grd_pred %>% dplyr::filter(LAT == 0.8),
         aes(x = LON, y = SHARP)) +
  geom_line(size = 1) +
  ggtitle("Sharp Predictions at 0.8 degrees N") +
  theme_minimal()

## ----smooth, fig.width = 5, fig.height = 4, fig.align='center'----------------
ggplot(toy_regions) +
  geom_sf() +
  geom_tile(data = grd_pred, aes(x = LON, y = LAT, fill = SMOOTH)) +
  scale_fill_viridis_c(limits = c(-0.3, 1)) +
  geom_hline(yintercept = 0.8) +
  ggtitle("Smooth Predictions", "Black line corresponds to x-axis of the next plot.") +
  xlab("") + ylab("") +
  theme_bw()

## ----smooth08, fig.width = 4, fig.height = 2.5, fig.align='center'------------
ggplot(grd_pred %>% dplyr::filter(LAT == 0.8),
         aes(x = LON, y = SMOOTH)) +
  geom_line(size = 1) +
  ggtitle("Smooth Predictions at 0.8 degrees N") +
  theme_minimal()

