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

# PO data: occ == 1 are presences, occ == 0 are background
maxent_data <- po_data[, c("occ", "elev", "tpi", "shrub", "temp", "precp",
                            "barren", "tri", "herb", "crop", "decBroad", "evBroad")]

SDM <- maxnet(
  p       = maxent_data$occ,
  data    = maxent_data[, -1],
  f       = maxnet.formula(p       = maxent_data$occ,
                           data    = maxent_data[, -1],
                           classes = "lq"),
  regmult = 1
)
