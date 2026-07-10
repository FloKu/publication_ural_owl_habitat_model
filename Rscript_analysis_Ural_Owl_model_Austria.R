# *************************************************************
# Spatial analyses of ural owl
# *************************************************************
# authors: Florian Kunz; Matthias Amon
# last changed: 06.2026 by FK, added chapter six parts on territories
# script to the publication

{library("dplyr")
  library("readxl")
  library("lubridate")
  library("terra")
  library("gdistance")
  library("tidysdm")
  library("sf")
  library("kuenm")
  library("spatstat")}

# load path variable
path <- "YOUR_PATH"

# load helper function
source(paste0(path, "/Rfunction_extract.by.mask.R")) # see other Github repository


### Chapter One - Creating presence data ###
############################################

# 0) Read in data from excel database
# **************************
database <- read_excel(paste0(path, "presence_points/Datenbank_Habichtskauz.xlsx"), 
                       col_types = c("text", "numeric", "numeric", "text", "text", "date", "numeric", "text")) %>% 
  mutate(source = as.factor(source)) %>% 
  mutate(type = as.factor(type)) %>% 
  mutate(year = year(date)) %>% 
  rownames_to_column(var="ID")

str(database)


# 1) Data checks
# **************************
# check number of NA values
for (i in colnames(database)) {
  print(paste0("column: ", i, " / Missing values: ", length(which(is.na(database[[i]])))))
}

# look for unusual coordinates
head(sort(database$longitude, decreasing=T))
tail(sort(database$longitude, decreasing=T))

head(sort(database$latitude, decreasing=T))
tail(sort(database$latitude, decreasing=T))

# year
head(sort(database$year, decreasing=T))
tail(sort(database$year, decreasing=T))

# identify and correct/exclude
database[database$longitude==45702.00000,]
database[database$longitude==30286,]
database[database$longitude==45.70200,]
database[database$year==2025,]
database[database$year==1946,]


# 2) Filtering
# **************************
str(database)
levels(database$source)

# number of obsverations per source
for (i in levels(database$source)) {
  print(paste0(nrow(database[database$source == i, ]),"   " ,i))
}

# number of obsverations per year
for (i in levels(as.factor(database$year))) {
  print(paste0(nrow(database[database$year == i, ]),"   " ,i))
}

# create spatvector from presence points
presence_shp <- vect(st_as_sf(database, coords=c("longitude", "latitude"), crs=4326))
# 839 records 


# 2.1) Filter wrong records 
# **************************
database <- database[database$longitude != 45702,]
database <- database[database$longitude != 45.70200,]
database <- database[database$longitude != 30286,]


# 2.2) Filter for study area 
# **************************
# exclude all points outside of the study area
study_area <- rast(paste0(path,"modelling/variables_asc/elevation.asc"))
values(study_area)[!is.na(values(study_area))] <- 1
study_area_proj <- project(as.polygons(study_area), crs(presence_shp))

# exclude all points outside of the study area
presence_shp <- terra::mask(presence_shp, study_area_proj)
# 797 records remaining

levels(presence_shp$source)
# 8 data sources present. Data from 5 sources has been reviewed and filtered by Richard Zink manually a priori during data gathering: 
#   App Wildtiere, Beobachtungen_HK, Naturbeobachter, Website Stadtwildtiere, Website Wilde Nachbarn
# 3 data sources need assessment and filtering: gbif_hk, Wildnissgebiet_hk, Birdlife


# 2.3) Filter the Austrian Ornithologist Centre and citizen science data  
# **************************
# was done a priori by R Zink and M Amon, no shape files needed


# 2.4) Filter GBIF data
# **************************
# Include only category of HUMAN_OBSERVATIONS 
# Also inspect points manually, as GBIF data is know to host unreliable data

# gbif types before year filtering
gbif <- presence_shp[presence_shp$source == "gbif_hk",]
levels(droplevels(gbif$type))

# filter data so that only HUMAN_OBSERVATIONS are included
gbif <- gbif %>% 
  subset(presence_shp$type != "MATERIAL_CITATION") %>% 
  subset(presence_shp$type != "PRESERVED_SPECIMEN")
levels(droplevels(gbif$type))

# save as shape and inspect manually
writeVector(gbif, filename=paste0(path, "presence_points/shapes/FILTERING_gbif_for_inspection.shp"), overwrite=TRUE)

# After inspection, we chose to exclude gbif points generally
# There are only 34 points, and about half of them seem unreliable (due to points being located in the midst of agricuztlural areas)

# exclude all gbif records
presence_shp <- subset(presence_shp, presence_shp$source != "gbif_hk")
# 763 records remaining


# 2.4) Filter BirdLife data
# **************************
# Filtering to precision: Matthias Amon only included datapoints from BirdLife data into our database that were noted as "Exakte Lokalisierung" / "Exact location"
# As for gbif, birdlife points needed manual inspection

# birdlife data into shape file and inspect manually
birdlife <- presence_shp[presence_shp$source == "birdlife",]
writeVector(birdlife, filename=paste0(path, "presence_points/shapes/FILTERING_birdlife_for_inspection.shp"), overwrite=TRUE)

# After inspection: one point excluded where commentary indicated faulty coordinates
presence_shp <- presence_shp[presence_shp$ID != 594,] # Birdlife datapoint that has a imprecise localization
# 762 records remaining


# 2.5) Filter years 
# **************************
# we chose only individuals that we habituated after reintroduction (so as year, we chose the year of the first successful reproduction)
presence_shp <- subset(presence_shp, presence_shp$year >= 2009)  
# 762 records remaining - no record excluded as old records were excluded in steps above already


# 2.6) Spatial filtering by nearest neighbor distance, depending on data source
# **************************
# The fopllowing filter steps will include randomization, as points will be exlcuded based on nearest neighbour distance
# Steps: 
#     1) filter citicen science data to resolve spatial clustering, using nearest meighbour distance
#     1) filter  Wildnissgebiet data (which is GPS locations per individual) using nearest neighbour distance to reduce individual bias (overrepresenation)
#     2) filter all data based on nearest neighbour distance to reduce spatial bias

# As these three steps involve random subsampling, we repeated all folowing filter steps 5 times and create 5 final datasets
# Models were calculated with all 5 datasets and then averaged to the final model

# how many individuals with gps
length(subset(presence_shp, presence_shp$source == "Wildnisgebiet_HK"))

# how many samples per gps logged individuum
table(subset(presence_shp$comment, presence_shp$source == "Wildnisgebiet_HK"))

# read raster on forest and shrub cover
forest <- rast(paste0(path, "presence_points/land_use/forest_shrubs.tif"))

# number of obsverations per source
presence_shp$source <- droplevels(presence_shp$source)
for (i in levels(presence_shp$source)) {
  print(paste0(nrow(presence_shp[presence_shp$source == i, ]),"   " ,i))
}

# looping filter steps 5 times to create 5 datasets
for (i in 1:5) {
  print(paste0("Iteration number ", i))
  
  # step 1: filter citizen science data
  # split data
  cs <- subset(presence_shp, presence_shp$source %in% c("App Wildtiere", "Birdlife", "Naturbeobachter", "Website Stadtwildtiere", "Website Wilde Nachbarn"))
  no_cs <- subset(presence_shp, !(presence_shp$source %in% c("App Wildtiere", "Birdlife", "Naturbeobachter", "Website Stadtwildtiere", "Website Wilde Nachbarn")))
  
  # NOT USED IN FINAL RUNS: exclude points not in forest area
  #pp_forest <- terra::extract(forest, cs)
  #cs$isForest <-pp_forest[,2]  
  
  #cs_thinned <- subset(cs, cs$isForest == 1)
  cs_thinned <- cs

  # filter based on nearest neighbor distance of 500m
  cs_thinned_sf <- st_as_sf(cs_thinned, coords = c("longitude", "latitude"))
  cs_thinned_sf_thinned <- tidysdm::thin_by_dist(cs_thinned_sf, dist_min=km2m(0.5))
  cs_thinned <- vect(cs_thinned_sf_thinned)
  print(paste0("Filtering citizen science: from ", nrow(cs), " to ", nrow(cs_thinned)))
  
  # add citizen science data back to data
  presence_shp_1 <- rbind(no_cs, cs_thinned)

  # step 2: filter gps data
  # divide data into two sets
  wg <- subset(presence_shp_1, presence_shp_1$source == "Wildnisgebiet_HK")
  no_wg <- subset(presence_shp_1, presence_shp_1$source != "Wildnisgebiet_HK")
  
  # create empty spatvector to be filled
  wg_thinned <- wg[-(1:nrow(wg)),]
  
  # apply spatial subsampling per individual, loop overall individuals and randomly exclude records within 500m
  for (j in 1:length(unique(wg$comment))) {
    IND <- unique(wg$comment)[j]
    wg_ind <- subset(wg, wg$comment == IND)
    wg_ind_sf <- sf::st_as_sf(wg_ind, coords = c("longitude", "latitude"))
    wg_ind_sf_thinned <- tidysdm::thin_by_dist(wg_ind_sf, dist_min=km2m(0.50))
    wg_ind_thinned <- vect(wg_ind_sf_thinned)
    print(paste0("Individuum ", IND, ": ", nrow(wg_ind), " to ", nrow(wg_ind_sf_thinned)))
    wg_thinned <- rbind(wg_thinned, wg_ind_thinned)
  }
  
  # add subsampled records of "Wildnisgebiet_HK" back to data
  presence_shp_fin <- rbind(no_wg, wg_thinned)
  print(paste0("Filtering GPS data: from ", nrow(wg), " to ", nrow(presence_shp_fin)))
 
  # step 3: apply spatial subsampling on all records, 100m
  presence_shp_fin_sf <- st_as_sf(presence_shp_fin, coords = c("longitude", "latitude"))
  presence_shp_fin_sf_thinned <- tidysdm::thin_by_dist(presence_shp_fin_sf, dist_min=km2m(0.1))
  presence_shp_fin_thinned <- vect(presence_shp_fin_sf_thinned)
  print(paste0("Retained records after spatial subsampling: ", nrow(presence_shp_fin_thinned)))
  
  # project to CS
  presence_shp_fin_thinned <- terra::project(presence_shp_fin_thinned, "epsg:32633")
  
  # export as shapefile for biasfile creattion
  file <- paste0(path, "presence_points/final_data/presence_points_thinned_iteration-",i,".shp") 
  writeVector(presence_shp_fin_thinned, filename=file, overwrite=TRUE)
  
  # export as csv for maxent
  maxent_input <- as.data.frame(presence_shp_fin_thinned) %>% 
    mutate(species="ural_owl") %>% 
    mutate(longitude=geom(presence_shp_fin_thinned)[,3]) %>% 
    mutate(latitude=geom(presence_shp_fin_thinned)[,4]) %>%
    dplyr::select("species", "longitude", "latitude")
  file <- paste0(path, "presence_points/final_data/presence_points_thinned_iteration-",i,".csv") 
  write.csv(maxent_input, file, row.names = F, quote=F)
}


### Chapter Two - prepare biasfile ###
############################################
# Due to RAM issues, biasfiles using Gaussian Kernel Density were created in ArcGIS SDM Toolbox
# Then, multiplied by 1000 and round to integer

# read in files
ref <- rast(paste0(path, "modelling/variables_asc/elevation.asc"))
files <- list.files(path=paste0(path, "modelling/biasfiles_SDMtoolbox"), pattern = "^biasfile_[0-9]+\\.asc$")

# project, multiply and round
for (i in 1:length(files)) {
  bias <- rast(paste0(path, "modelling/biasfiles_SDMtoolbox/", files[i]))
  
  # define projection
  crs(bias) <- "EPSG:32633"
  
  # mulitply and round
  bias <- round(bias*1000)
  
  # cleaning step: making sure extent, CS and raster correspond
  bias <- Rfunction.extract.by.mask(bias, ref, F)
  
  # save
  filepath <- paste0(path, "modelling/biasfiles_clean/biasfile_clean_it-", i, "_.asc")
  writeRaster(bias, filepath, NAflag=-9999, overwrite=TRUE)
}

# Duplicate biasfile of first iteration for the kuenm run
bias <- rast(paste0(path, "modelling/biasfiles_SDMtoolbox/", files[1]))
crs(bias) <- "EPSG:32633"
bias <- round(bias*1000)
bias <- Rfunction.extract.by.mask(bias, ref, F)
writeRaster(bias, paste0(path, "modelling/2_step_kuenm/bias.asc"), NAflag=-9999, overwrite=TRUE) # note that kuenm seems to need the biasfile be named exactly bias.asc, nothing else works


### Chapter Three - running kuenm ###
############################################

rm(list=setdiff(ls(), "path"))


# 1) write csv for kuenm run
# **************************
setwd(paste0(path, "modelling/2_step_kuenm"))

occ_kauz <- read.csv(paste0(path, "presence_points/final_data/presence_points_thinned_iteration-1.csv"), header=TRUE)
split <- kuenm_occsplit(occ_kauz, train.proportion = 0.75, method = "random", save = TRUE, name = "occ")


# 2) run kuenm
# **************************
setwd(paste0(path, "modelling/2_step_kuenm"))

# preparing sets of variables
help("kuenm_varcomb")

kuenm_varcomb(var.dir="Variables_asc", out.dir="M_variables_2404",
              min.number=5, in.format="ascii", out.format="ascii")

# create candidate models
help("kuenm_cal")

oj <- "occ_joint.csv"
otr <- "occ_train.csv"
mvars <- "M_variables_2404"
candir <- "Candidate_Models_2404"
bcal <- "Candidate_Models_2404"
regm <- c(1,2,3)
maxpath <- "C:/Users/floriankunz/Desktop/maxent"
max.memory <- 1300
fclas <- "lqh"

kuenm_cal(occ.joint = oj, occ.tra = otr, M.var.dir = mvars, batch = bcal,
          out.dir = candir, reg.mult = regm, f.clas = fclas, 
          args = "biasfile=C:/Users/floriankunz/seadrive_root/Florian_1/Für mich freigegeben/Habichtskauz - Resubmission/modelling/2_step_kuenm/bias.asc biastype=3 maximumiterations=5000",
          max.memory = max.memory, maxent.path = maxpath, wait = TRUE, run = TRUE)

# evaluate candidate models
help(kuenm_ceval)

ote <- "occ_test.csv"
path <- "Candidate_Models_2404"
cresdir <- "Calibration_results_2404"
treshold <- 5

Final_evaluation <- kuenm_ceval(path=path, occ.joint=oj, occ.tra=otr, occ.test=ote,
                                batch=bcal, out.eval=cresdir, threshold=treshold, rand.percent=50, iterations=5000,
                                kept=TRUE, selection="OR_AICc")

# create final models
help(kuenm_mod)

bfmod <- "Final_Models_2404"
moddir <- "Final_Models_2404"

kuenm_mod(occ.joint = oj, M.var.dir = mvars, out.eval = cresdir, maxent.path = maxpath, 
          out.dir = moddir, batch = bfmod, rep.n = 20, rep.type = "Bootstrap", 
          jackknife = TRUE, max.memory = 1300, out.format = "cloglog", project = FALSE,
          write.mess = FALSE, write.clamp = FALSE, wait = TRUE, run = TRUE,
          args = "biasfile=C:/Users/floriankunz/seadrive_root/Florian_1/Für mich freigegeben/Habichtskauz - Resubmission/modelling/2_step_kuenm/bias.asc biastype=3 maximumiterations=5000")


### Chapter Four - averaging models ### 
############################################

# read in model files
models_vec <- c()

for (i in 1:5) {
  model <- rast(paste0(path, "modelling/3_step_final_models/final_model_", i, "/ural_owl_avg.asc"))
  
  # define projection
  crs(model) <- "EPSG:32633"

  # add to vector
  models_vec <- append(models_vec, model, after=length(models_vec))
  
  # save as tif for visualization
  filepath <- paste0(path, "modelling/3_step_final_models/final_model_avg/final_model_", i, ".tif")
  writeRaster(model, filepath, overwrite=TRUE)
}
  

# average
model_fin <- mean(models_vec)

# write final model
writeRaster(model_fin, paste0(path, "modelling/3_step_final_models/final_model_avg/final_model.asc"), NAflag=-9999, overwrite=TRUE)
writeRaster(model_fin, paste0(path, "modelling/3_step_final_models/final_model_avg/final_model.tif"), overwrite=TRUE)

### Chapter Five - summarize results ### 
############################################

# initialize vectors
{testAUC <- c()
  PI_decid <- c()
  PI_timber <- c()
  PI_slope <- c()
  PI_edge <- c()
  PI_aspect <- c()
  PC_decid <- c()
  PC_timber <- c()
  PC_slope <- c()
  PC_edge <- c()
  PC_aspect <- c()
  binary <- c()}

# read in maxent results
for (i in 1:5) {
  maxent <- read.csv(paste0("E:/Habichtskauz_final/3_step_final_models/final_model_", i, "/maxentResults.csv"))
  
  # append to vector
  testAUC <- append(testAUC, maxent$Test.AUC[1:20], after=length(testAUC))
  
  PI_decid <- append(PI_decid, maxent$proportion_of_deciduous_forest.permutation.importance[1:20], after=length(PI_decid))
  PI_timber <- append(PI_timber, maxent$proportion_of_timber_wood.permutation.importance[1:20], after=length(PI_timber))
  PI_slope <- append(PI_slope, maxent$slope.permutation.importance[1:20], after=length(PI_slope))
  PI_edge <- append(PI_edge, maxent$proportion_of_forest_edge.permutation.importance[1:20], after=length(PI_edge))
  PI_aspect <- append(PI_aspect, maxent$cat_aspect.permutation.importance[1:20], after=length(PI_aspect))
  
  PC_decid <- append(PC_decid, maxent$proportion_of_deciduous_forest.contribution[1:20], after=length(PC_decid))
  PC_timber <- append(PC_timber, maxent$proportion_of_timber_wood.contribution[1:20], after=length(PC_timber))
  PC_slope <- append(PC_slope, maxent$slope.contribution[1:20], after=length(PC_slope))
  PC_edge <- append(PC_edge, maxent$proportion_of_forest_edge.contribution[1:20], after=length(PC_edge))
  PC_aspect <- append(PC_aspect, maxent$cat_aspect.contribution[1:20], after=length(PC_aspect))
  
  binary <- append(binary, maxent$Maximum.test.sensitivity.plus.specificity.Cloglog.threshold[1:20], after=length(binary))
}

# averages 
mean(testAUC)

mean(PI_decid)
mean(PI_timber)
mean(PI_slope)
mean(PI_aspect)
mean(PI_edge)

mean(PC_decid)
mean(PC_timber)
mean(PC_slope)
mean(PC_aspect)
mean(PC_edge)

mean(binary)


### Chapter Six - Zonal statitics, post modelling  ### 
############################################

path <- "C:/Users/mett_meister/Documents/Habichtskauz_lokal/ArcGIS/"
path <- "C:/Users/floriankunz/seadrive_root/Florian_1/Für mich freigegeben/Habichtskauz - Resubmission/modelling/3_step_final_models/final_model_avg/"

# read in conservation area files

biosphere_reserve <- vect(paste0(path,"Biospharenpark_Wien_Kernzonen/RRU_BPWW_KERNZONENPolygon_UTM33N.shp"))
protected_landscape_vi <- vect(paste0(path, "Naturschutzgebiete_Wien/NATURSCHUTZGEBOGDPolygon_UTM33N.shp"))
protected_landscape_noe <- vect(paste0(path, "Naturschutzgebiete_NOE/RNA_NSGEBPolygon_UTM33N.shp"))
natura_2000_vi <- vect(paste0(path, "Natura2000_Wien/NATURA2TOGDPolygon_UTM33N.shp"))
natura_2000_noe <- vect(paste0(path, "Natura2000_NOE/RNA_N2K_FFHPolygon_UTM33N.shp"))

# Merge all protected area shapefiles into one vector
all_protected_areas <- rbind(
  biosphere_reserve,
  protected_landscape_vi,
  protected_landscape_noe,
  natura_2000_vi,
  natura_2000_noe
)

# Load final model (tif)
final_model_avg <- rast(paste0(path, "final_model.tif"))

# 1. Caluate how much of the area is classified as high habitat potential 
##########################################################################

# truncate the final model with a threshold to create a binary map of only high habitat potential (threshold of maximum test sensitivity plus specificity)

final_model_avg[final_model_avg < 0.490] <- NA

# Calculate percentage of values ≥ 0.490 
original_values <- values(final_model_avg)

valid_values <- original_values[!is.na(original_values)]

suitable_count <- sum(valid_values >= 0.490)

valid_total <- length(valid_values)

# Calculate the percentage
percent_suitable <- (suitable_count / valid_total) * 100

# Output the result
cat("Percentage of valid raster cells ≥ 0.490 :", round(percent_suitable, 2), "%\n")

# 2. calulate how much of the high habitat potential area is also currently under protection 
#############################################################################################

# Crop and mask the raster with the merged protected areas
masked_raster <- mask(crop(final_model_avg, all_protected_areas), all_protected_areas)

# Count non-NA cells in the entire raster
total_area_cells <- sum(!is.na(values(final_model_avg)))

# Count non-NA cells in the overlapping (masked) area
overlap_area_cells <- sum(!is.na(values(masked_raster)))

# Calculate percentage overlap
percent_overlap <- (overlap_area_cells / total_area_cells) * 100

# Print result
cat("Percentage of raster overlapping with protected areas:", round(percent_overlap, 2), "%\n")

# 3. calulate the number of territories supported 
#############################################################################################
path <- "C:/Users/floriankunz/seadrive_root/Florian_1/Für mich freigegeben/Habichtskauz - Resubmission/modelling/3_step_final_models/"

m1 <- rast(paste0(path, "final_model_1/ural_owl_avg.asc"))
m2 <- rast(paste0(path, "final_model_2/ural_owl_avg.asc"))
m3 <- rast(paste0(path, "final_model_3/ural_owl_avg.asc"))
m4 <- rast(paste0(path, "final_model_4/ural_owl_avg.asc"))
m5 <- rast(paste0(path, "final_model_5/ural_owl_avg.asc"))
avg <- rast(paste0(path, "final_model_avg/final_model.tif"))

models <- list(m1, m2, m3, m4, m5)

# total suitable area
avg_bin <- avg
avg_bin[avg_bin < 0.49] <- NA
avg_bin[avg_bin >= 0.49] <- 1

# total suitable area size in nectare: 
terra::global(avg_bin, fun="notNA")

# divided through minimum territory size: 
terra::global(avg_bin, fun="notNA")/260
# The total area supports 1942 territories 

# function to calculate number of territories
n_territories <- function (model, th, terr_size) { # terr_size in hectares
  
  # create binary map per iteration
  binary <- model
  binary[binary < th] <- 0
  binary[binary >= th] <- 1
  
  # check if patches occur
  clusters <- patches(binary, directions = 8, zeroAsNA=TRUE)  # 8 directions for connectivity (queen's case)
  cluster_sizes <- freq(clusters)  # Frequency table of cluster IDs and sizes
  
  # remove all patches that dont meet the criteria of a minimum habitat of [terr_size] km2
  if (nrow(cluster_sizes) != 1) {
    big_ids <- cluster_sizes$value[cluster_sizes$count >= terr_size]
    binary[!(clusters %in% big_ids)] <- NA
  }
  
  # calculate number of territories
  binary[binary == 0] <- NA
  n <- terra::global(binary, fun="notNA")/terr_size
  return(n)
}

# how many territories based on average run
model <- avg
th <- 0.49
terr_size <- 260

n_territories(model, th, terr_size)
# While the total area supports 1942 territories, when accounting for patch size, 1534.8 territories are supported

# run function on all 5 repeats
territories <- c()

for (i in 1:5) {
  
  print(i)
  
  # read in data
  hsi <- models[[i]]
  sum <- read.csv(paste0(path, "/maxentResults_", i, ".csv"))
  binary <- sum$Maximum.test.sensitivity.plus.specificity.Cloglog.threshold[21]
  
  # calculate territories
  n <- n_territories(hsi, binary, 260)
  print(n)
  
  territories <- c(territories, n$notNA)
}
  
# results
mean(territories)
sd(territories)
# mean: 1531 sd: 111.7

#########################
