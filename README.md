# Using-Integrated-Species-Distribution-Models
There are 10 datasets that represent five raptors species used in our integrated species distribution modeling efforts. Each species has two datasets, the ones labeled with "PA" are the presences-absence datasets, while the "PO" are the presence-only datasets. Please see the metadata for detailed information about both types of datasets. 

File names:
1.	AustralPygmyOwl_PA.csv
2.	AustralPygmyOwl_PO.csv
3.	RufousLeggedOwl_PA.csv
4.	RufousLeggedOwl_PO.csv
5.	RufousTailedHawk_PA.csv
6.	RufousTailedHawk_PO.csv
7.	StriatedCaracara_PA.csv
8.	StriatedCaracara_PO.csv
9.	WhiteThroatedCaracara_PA.csv
10.	WhiteThroatedCaracara_PO.csv

Column Names:
_PA files are the eBird presence absence data. PA data has 13 years therefore all detection data will have 13 numbered columns. Not all cells have an eBird survey, those cell/years will contain a NA instead of 0 or 1.

_PO files are the presence only data and only have one column for detection covariates. 

“site” = grid cell number (not comparable between species ONLY within species datasets)

“dummy” = sampling process for presence-only data 

“occ” = the presence (1) and background or absences depending on dataset (0)

“dur” = total duration of minutes for eBird data

“time” = time of day to the nearest hour

“nobs” = number of observers

“distRoad” = distance to nearest road

“distCoast” = distance to coastline for occupancy of Striated Caracara only

“jday” = julian day

“year” = year (1-13)

“temp” = either annual mean temperature or max temperature warmest month depending on the species

“precip” = either precipitation of warmest quarter or precipitation seasonality

“barren” = % cover of barren lands

“crop” = % cover of cropland

“decBroad” = % cover deciduous broadleaf forest

“mixed” = % cover mixed forest

“evBroad” = % cover evergreen broadleaf forest

“herb” = % cover herbaceous vegetation

“shrub” = % cover shrubland

“elev” = elevation

“tpi” = topographical position index

“tri” = terrain roughness index

“X” and “Y” = coordinates of grid cell center. Spatial data were reprojected to Albers Equal Area projection centered at 60°W and 32°S, with standard parallels at 5°S and 42°S. The projection uses the Australian National & South American 1969 ellipsoid with units defined in kilometers

“species” = the species the data file represents

“resolution” = the spatial resolution of the grid cell in kilometers



