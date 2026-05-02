# ============================================================
# Integrated Species Distribution Model (iSDM) and MaxEnt
# Species: Austral Pygmy-Owl (Glaucidium nana)
#
# Required input files (available in data repository):
#   AustralPygmyOwl_PA.csv  — eBird presence-absence data
#   AustralPygmyOwl_PO.csv  — Presence-only (GBIF + external) data
#
# Column descriptions:
#   site           = grid cell ID (links PA and PO datasets)
#   dummy          = sampling process indicator for presence-only data
#   occ            = presence (1) / background (0) — PO file
#   occ1–occ13     = eBird detection per breeding season (PA file)
#                    NA = grid cell not surveyed that season
#   dur.1–dur.13   = eBird checklist duration per season (minutes)
#   jday.1–jday.13 = Julian day per season
#   time.1–time.13 = time of day per season (nearest hour)
#   nobs.1–nobs.13 = number of observers per season
#   temp           = mean annual temperature
#   precp          = precipitation of warmest quarter
#   barren         = % cover barren lands
#   crop           = % cover cropland
#   decBroad       = % cover deciduous broadleaf forest
#   evBroad        = % cover evergreen broadleaf forest
#   herb           = % cover herbaceous vegetation
#   shrub          = % cover shrubland
#   elev           = elevation (m)
#   tpi            = topographic position index
#   tri            = terrain roughness index
#   X, Y           = grid cell centre coordinates (Albers Equal Area, km)
#                    Projection: +proj=aea +lat_0=-32 +lon_0=-60
#                                +lat_1=-5 +lat_2=-42
#                                +ellps=aust_SA +units=km
#   species        = species name
#   resolution     = spatial resolution of grid cell (km)
# ============================================================

library(spOccupancy)
library(maxnet)
library(ROCR)
library(raster)
library(redlistr)
library(downscale)
library(sf)

# ---- Load data ----
pa_data <- read.csv("AustralPygmyOwl_PA.csv")
po_data <- read.csv("AustralPygmyOwl_PO.csv")

# ============================================================
# DATA PREPARATION
# ============================================================

# ---- Build combined coordinate grid from union of PA and PO sites ----
# coords must cover ALL unique sites from both datasets combined
# because spIntPGOcc uses row positions in coords as the site index
coords_all <- rbind(
  pa_data[, c("site", "X", "Y")],
  po_data[!po_data$site %in% pa_data$site, c("site", "X", "Y")]
)
coords_all  <- coords_all[order(coords_all$site), ]
coords      <- as.matrix(coords_all[, c("X", "Y")])

# ---- Occupancy covariates — one row per site in coords_all ----
occ_all <- rbind(
  pa_data[, c("site", "elev", "tpi", "shrub", "temp", "precp", "distRoad",
              "barren", "tri", "herb", "crop", "decBroad", "evBroad")],
  po_data[!po_data$site %in% pa_data$site,
          c("site", "elev", "tpi", "shrub", "temp", "precp", "distRoad",
            "barren", "tri", "herb", "crop", "decBroad", "evBroad")]
)
occ_all  <- occ_all[order(occ_all$site), ]

occ_covs <- data.frame(
  elev     = occ_all$elev,
  tpi      = occ_all$tpi,
  shrub    = occ_all$shrub,
  temp     = occ_all$temp,
  precp   = occ_all$precp,
  barren   = occ_all$barren,
  tri      = occ_all$tri,
  herb     = occ_all$herb,
  crop     = occ_all$crop,
  decBroad = occ_all$decBroad,
  evBroad  = occ_all$evBroad
)

# ---- Site indices — row positions in coords_all ----
# NOT the raw site IDs — spIntPGOcc expects row positions
sites_int <- list(
  extGB = match(po_data$site, coords_all$site),
  ebird = match(pa_data$site, coords_all$site)
)

# ---- Presence-absence (eBird) ----
# y matrix: rows = grid cells, columns = breeding seasons (1-13)
# NA = grid cell not surveyed that season
y_pa <- as.matrix(pa_data[, paste0("occ", 1:13)])

det_pa <- list(
  dur  = as.matrix(pa_data[, paste0("dur.",  1:13)]),
  jday = as.matrix(pa_data[, paste0("jday.", 1:13)]),
  time = as.matrix(pa_data[, paste0("time.", 1:13)]),
  nobs = as.matrix(pa_data[, paste0("nobs.", 1:13)]),
  site = as.matrix(pa_data[, paste0("site")]),  
  distRoad = as.matrix(pa_data[, paste0("distRoad")]),
  year = as.matrix(pa_data[, paste0("year.", 1:13)])  
)

# ---- Presence-only (GBIF + external) ----
y_po <- matrix(po_data$occ, ncol = 1)

det_po <- list(
  dummy = matrix(po_data$dummy, ncol = 1),
  distRoad = matrix(po_data$distRoad, ncol = 1)
)

# ---- Integrated data list ----
data_int <- list(
  y        = list(extGB = y_po,   ebird = y_pa),
  occ.covs = occ_covs,
  det.covs = list(extGB = det_po, ebird = det_pa),
  sites    = sites_int,
  coords   = coords
)

# ============================================================
# iSDM — Integrated Spatial Occupancy Model
# ============================================================

occ_formula <- ~ scale(elev) + I(scale(elev)^2) +
  scale(tpi) + scale(shrub) +
  scale(temp) + scale(precp) +
  scale(barren) + scale(tri) +
  scale(herb) + scale(crop) +
  scale(decBroad) + scale(evBroad)

det_formula <- list(
  extGB = ~ dummy + scale(distRoad),
  ebird = ~ scale(dur) + scale(jday) + scale(time) +
    scale(nobs) + scale(distRoad) + (1|site) + (1|year)
)

# Priors
prior_list <- list(
  beta.normal  = list(mean = rep(0, 13),
                      var  = c(1, rep(2.72, 12))),
  alpha.normal = list(mean = list(0, 0),
                      var  = list(c(5, 5, 2.72),
                                  c(5, 2.72, 2.72, 2.72, 2.72, 2.72))),
  sigma.sq.ig  = c(20, 10),
  phi.unif     = c(0.8, 3.75)
)

# Initial values
J <- nrow(coords)
inits_list <- list(
  alpha    = list(0, 0),
  beta     = rep(0, 13),
  z        = rep(1, J),
  sigma.sq = 0.5,
  phi      = 2,
  w        = rep(0, J)
)

# MCMC settings
# Total iterations per chain = n.batch * batch.length = 375,000
# Effective posterior samples per chain = (375,000 - n.burn) / n.thin = 2,000
n_batch      <- 15000
batch_length <- 25
n_burn       <- 50000
n_thin       <- 50

iSDM <- spIntPGOcc(
  occ.formula   = occ_formula,
  det.formula   = det_formula,
  data          = data_int,
  inits         = inits_list,
  priors        = prior_list,
  tuning        = list(phi = 0.2),
  cov.model     = "exponential",
  NNGP          = TRUE,
  n.neighbors   = 5,
  n.batch       = n_batch,
  batch.length  = batch_length,
  n.burn        = n_burn,
  n.thin        = n_thin,
  n.chains      = 3,
  n.omp.threads = 4,
  n.report      = 200
)

# ============================================================
# MaxEnt
# ============================================================
# MaxEnt uses presence-background data combining:
#   1. GBIF + external presences (from po_data)
#   2. eBird presences — any grid cell with at least one detection
#      across any of the 13 breeding seasons
#   Background points are all non-presence grid cells in the extent

# ---- Identify eBird presence cells ----
# A cell is a presence if detected at least once in any season
ebird_pres_sites <- pa_data$site[apply(pa_data[, paste0("occ", 1:13)], 1,
                                        function(x) any(x == 1, na.rm = TRUE))]

# Build eBird presence rows with covariates from pa_data
# Include ALL ebird presence sites not already marked as presence in po_data
# (some may exist as background points in po_data and need upgrading)
ebird_pres_df <- pa_data[pa_data$site %in% ebird_pres_sites,
                          c("site", "temp", "barren", "crop", "decBroad", "elev",
                            "evBroad", "herb", "precp", "shrub", "tpi", "tri", "X", "Y")]
ebird_pres_df$occ <- 1

# Sites already marked as presence in po_data — no need to add again
already_pres <- po_data$site[po_data$occ == 1]
ebird_new    <- ebird_pres_df[!ebird_pres_df$site %in% already_pres, ]

# Combine: start with all po_data, add new ebird presences,
# then for any site appearing twice keep the presence (occ=1)
maxent_combined <- bind_rows(
  po_data[, c("site", "occ", "temp", "barren", "crop", "decBroad", "elev",
               "evBroad", "herb", "precp", "shrub", "tpi", "tri", "X", "Y")],
  ebird_new
) %>%
  group_by(site) %>%
  slice(which.max(occ)) %>%   # keeps occ=1 where site appears as both background and presence
  ungroup() %>%
  as.data.frame()

cat("Total presences for MaxEnt:", sum(maxent_combined$occ == 1), "\n")
cat("Total background for MaxEnt:", sum(maxent_combined$occ == 0), "\n")
cat("Total cells:", nrow(maxent_combined), "\n")

# ---- Fit MaxEnt ----
# Column order matches original: temp, barren, crop, decBroad, elev,
# evBroad, herb, precp, shrub, tpi, tri
maxent_data <- maxent_combined[, c("occ", "temp", "barren", "crop", "decBroad",
                                    "elev", "evBroad", "herb", "precp", "shrub",
                                    "tpi", "tri")]

SDM <- maxnet(
  p       = maxent_data$occ,
  data    = as.data.frame(maxent_data[, -1]),
  f       = maxnet.formula(p       = maxent_data$occ,
                           data    = as.data.frame(maxent_data[, -1]),
                           classes = "lq"),
  regmult = 1
)

# ============================================================
# MAXENT PREDICTIONS AND EOO/AOO
# ============================================================
# Predict on ALL grid cells in the study extent using the fitted model.
# Background cells (occ=0) are not true absences — they characterise
# the available environmental space for the species.

pred_covs <- as.data.frame(maxent_combined[, c("temp", "barren", "crop", "decBroad",
                                                 "elev", "evBroad", "herb", "precp",
                                                 "shrub", "tpi", "tri")])
pred_vals <- predict(SDM, pred_covs, type = "cloglog", clamp = TRUE)

# Build prediction raster from grid coordinates
pred_df     <- data.frame(x = maxent_combined$X, y = maxent_combined$Y, pred = pred_vals)
pred_raster <- rasterFromXYZ(pred_df, res = 25,
                              crs = "+proj=aea +lat_0=-32 +lon_0=-60 +lat_1=-5 +lat_2=-42 +x_0=0 +y_0=0 +ellps=aust_SA +units=km +no_defs")

# ---- Threshold 1: maxSSS (maximum sum of sensitivity + specificity) ----
pred_rocr     <- prediction(pred_vals, maxent_combined$occ)
ss            <- performance(pred_rocr, "sens", "spec")
thresh_idx    <- which.max(ss@x.values[[1]] + ss@y.values[[1]])
thresh_maxSSS <- ss@alpha.values[[1]][thresh_idx]
cat("MaxEnt maxSSS threshold:", round(thresh_maxSSS, 3), "\n")

# ---- Threshold 2: mean predicted value across ALL grid cells ----
# meanPred = global mean of the prediction surface
# equivalent to terra::global(pred.mn, fun = "mean", na.rm = TRUE)
thresh_meanPred <- mean(pred_vals, na.rm = TRUE)
cat("MaxEnt meanPred threshold:", round(thresh_meanPred, 3), "\n")

# ============================================================
# MaxEnt EOO and AOO — maxSSS threshold
# ============================================================

# Reclassify: cells >= threshold = 1 (suitable), < threshold = 0 (unsuitable)
rcmat_maxSSS <- matrix(c(0, thresh_maxSSS, 0,
                          thresh_maxSSS, 1, 1), ncol = 3, byrow = TRUE)
rcm_maxSSS   <- reclassify(pred_raster, rcmat_maxSSS)

# AOO — downscale from 25km cells to 4km² (2x2km) using ensemble of models
occupancy_maxSSS <- upgrain(atlas.data = rcm_maxSSS, cell.width = 25,
                             scales = 3, method = "All_Sampled", plot = FALSE,
                             return.rasters = FALSE)

ens_maxSSS <- ensemble.downscale(occupancies = occupancy_maxSSS,
                                  new.areas = c(4),
                                  models = c("Nachman", "PL", "Logis", "GNB",
                                             "FNB", "Hui", "Poisson", "INB", "NB"),
                                  verbose = FALSE, plot = FALSE)
AOO_maxSSS <- ens_maxSSS$AOO$Means
cat("MaxEnt AOO (maxSSS):", AOO_maxSSS, "km²\n")

# EOO — minimum convex polygon around suitable cells
rcmat_eoo      <- matrix(c(0, 0, NA, 1, 1, 1), ncol = 3, byrow = TRUE)
rcm_maxSSS_eoo <- reclassify(rcm_maxSSS, rcmat_eoo, include.lowest = TRUE)
EOO_poly_maxSSS  <- makeEOO(rcm_maxSSS_eoo)
EOO_area_maxSSS  <- getAreaEOO(EOO_poly_maxSSS)
cat("MaxEnt EOO (maxSSS):", EOO_area_maxSSS, "km²\n")

# Save EOO polygon as sf object
E_maxSSS <- st_as_sf(EOO_poly_maxSSS)

# ============================================================
# MaxEnt EOO and AOO — meanPred threshold
# ============================================================

rcmat_meanPred <- matrix(c(0, thresh_meanPred, 0,
                            thresh_meanPred, 1, 1), ncol = 3, byrow = TRUE)
rcm_meanPred   <- reclassify(pred_raster, rcmat_meanPred)

occupancy_meanPred <- upgrain(atlas.data = rcm_meanPred, cell.width = 25,
                               scales = 3, method = "All_Sampled", plot = FALSE,
                               return.rasters = FALSE)

ens_meanPred <- ensemble.downscale(occupancies = occupancy_meanPred,
                                    new.areas = c(4),
                                    models = c("Nachman", "PL", "Logis", "GNB",
                                               "FNB", "Hui", "Poisson", "INB", "NB"),
                                    verbose = FALSE, plot = FALSE)
AOO_meanPred <- ens_meanPred$AOO$Means
cat("MaxEnt AOO (meanPred):", AOO_meanPred, "km²\n")

rcm_meanPred_eoo <- reclassify(rcm_meanPred, rcmat_eoo, include.lowest = TRUE)
EOO_poly_meanPred  <- makeEOO(rcm_meanPred_eoo)
EOO_area_meanPred  <- getAreaEOO(EOO_poly_meanPred)
cat("MaxEnt EOO (meanPred):", EOO_area_meanPred, "km²\n")

E_meanPred <- st_as_sf(EOO_poly_meanPred)

# ============================================================
# PRINT SUMMARY
# ============================================================
cat("\n=== MaxEnt Results ===\n")
cat("Threshold maxSSS:   ", round(thresh_maxSSS,   3), "\n")
cat("Threshold meanPred: ", round(thresh_meanPred, 3), "\n")
cat("AOO (maxSSS):       ", AOO_maxSSS,   "km²\n")
cat("AOO (meanPred):     ", AOO_meanPred, "km²\n")
cat("EOO (maxSSS):       ", EOO_area_maxSSS,   "km²\n")
cat("EOO (meanPred):     ", EOO_area_meanPred, "km²\n")

# ============================================================
# iSDM EOO AND AOO
# ============================================================
# Uses posterior samples from the fitted iSDM to calculate EOO and AOO
# across 2000 posterior samples, then summarises with mean and 90% CI
#
# Each posterior sample gives a continuous occurrence probability surface
# which is converted to a binary presence/absence raster before calculating
# EOO (minimum convex polygon) and AOO (downscaled to 4km² cells)

# Convert iSDM coordinates from km to m (rasterFromXYZ expects metres)
xy      <- as.data.frame(iSDM$coords)
colnames(xy) <- c("x", "y")
xy$x    <- xy$x * 1000
xy$y    <- xy$y * 1000

# Reference raster built from model coordinates — no external file needed
# coords are in km so convert to m for rasterFromXYZ
r_ref <- rasterFromXYZ(
  data.frame(x = coords[, 1] * 1000,
             y = coords[, 2] * 1000,
             z = 1),
  crs = "+proj=aea +lat_0=-32 +lon_0=-60 +lat_1=-5 +lat_2=-42 +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs"
)

# Sample indices — 2000 evenly spaced posterior samples across all chains
# Total samples = n.batch * batch.length - n.burn = 325,000 / n.thin = 2,000 per chain x 3 chains
k <- round(seq(from = 1, to = 19500, by = 9.75), 0)

eoo_samples <- aoo_samples <- list()

for (i in 1:length(k)) {

  # Extract one posterior sample of occupancy probability for all sites
  p1   <- as.data.frame(iSDM$z.samples[k[[i]], 1:3708])
  data <- cbind(xy, p1)

  # Build raster from this posterior sample
  rast_i <- rasterFromXYZ(data, res = res(r_ref), crs = crs(r_ref))

  # ---- AOO: downscale occupied 25km cells to 4km² ----
  occupancy_i <- upgrain(atlas.data = rast_i, cell.width = 25,
                          scales = 3, method = "All_Sampled", plot = FALSE,
                          return.rasters = FALSE)

  ens_i <- ensemble.downscale(occupancies = occupancy_i,
                               new.areas   = c(4),
                               models      = c("Nachman", "PL", "Logis", "GNB",
                                               "FNB", "Hui", "Poisson", "INB", "NB"),
                               verbose = FALSE, plot = FALSE)

  aoo_samples[[i]] <- ens_i$AOO$Means

  # ---- EOO: reclassify to binary then compute MCP ----
  rcmat_i <- matrix(c(0, 0, NA, 1, 1, 1), ncol = 3, byrow = TRUE)
  rast_i  <- reclassify(rast_i, rcmat_i, include.lowest = TRUE)

  eoo_poly_i       <- makeEOO(rast_i)
  eoo_samples[[i]] <- getAreaEOO(eoo_poly_i)

}

# ---- Summarise across posterior samples ----
eoo_vec <- unlist(eoo_samples)
aoo_vec <- unlist(aoo_samples)

cat("\n=== iSDM Results ===\n")
cat("EOO mean:",   round(mean(eoo_vec), 1), "km²\n")
cat("EOO 90% CI:", round(quantile(eoo_vec, 0.05), 1), "–",
                   round(quantile(eoo_vec, 0.95), 1), "km²\n")
cat("AOO mean:",   round(mean(aoo_vec), 1), "km²\n")
cat("AOO 90% CI:", round(quantile(aoo_vec, 0.05), 1), "–",
                   round(quantile(aoo_vec, 0.95), 1), "km²\n")
