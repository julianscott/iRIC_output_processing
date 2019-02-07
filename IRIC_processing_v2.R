##################################################
## Script purpose:  iRIC ouput processing for building an inundation map (i.e. dagwood) for ascending discharge levels for 
#                   species distribution modeling
## Date:  02/07/2019
## Author: Julian Scott
##################################################

# Requirements for running this code:

# 1. This code uses iric output using Calculation Result -> Export > 2D -> Format = CSV files
# 2. The sequential layers of inundating discharge surfaces has an order that is inherited when it is read in from the directory.
#    Therefore, the naming convention of the IRIC output csv files must be such that they are ordered correctly 
#    E.g., when viewing the csvs in File Explorer, when the folder is sorted by name, are they in the correct ascending order?
#    We achieve this using the following convention, IRICoutput_0000.01, IRICoutput_0001.5, IRICoutput_0010.0,IRICoutput_1010.0, etc
#    There are many other ways to do this, but this is how the code we have written works. 

# 3. change setwd to location where all of the downloaded test files are
# 4. When ready to use on your own dataset, there are multiple things that will have to be changed
#   - coordinate projection system (see variable projr)
#   - directories (see setwd)
#   - resolution (see res)
# 5. Other variables are inherited from the input, so they are dynamic and may not need to be changed.

# This is version 2 and is a work in progress. Please keep in touch with developments as you use this or other code for processing iric data.
# Updates from version 1 include:
# 1. Inherit extent and other raster information from one of the iric output csvs, rather than the DEM. 
# 2. Rather than using cover to combine layers of cells inundated at ascending sequential discharges, 
#     use raster::boundaries and raster::merge to combine layers of the boundary of inundation of asencding sequential discharges.
#     Then, the "white" space in between boundaries may be interpolated/triangulated by exporting cells as points and TINing. 
#     I used arc gis for TINing, then making a digital discharge map (DQM)

#  packages used
packages <- c("SDMTools","sp","raster","rgeos","rgdal","sf","spatstat","spdep","tidyverse","rasterVis","ggplot2")

#  Check to see if each is installed, and install if not.
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {    
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

# load the installed libraries in the packages list 
lapply(packages,library,character.only=TRUE)

# Set this to the IRIC_Processing_in_R_v1 folder to run this code
setwd("E:\\_DoD\\_Camp_Pendleton_Survey\\IRIC\\_Modeling_dir\\IRIC_Processing_in_R_v1")

# get names of iric output csv files in the working directory
iric_results <- list.files(pattern = ".csv")

# create an empty raster with the extent, resolution, and projection of the iric mesh.
forextent <- read.csv(iric_results[1],skip = 2,header = T) 
e <- extent(with(forextent,c(min(X),max(X),min(Y),max(Y))))
# The dem we input to iRIC was 0.5x0.5 m resolution. The iRIC mesh was set to match this.
# Note that, in the test case, sinuousity of valley distorts cell sizes, so iric coordinate output has roughly 
# the same resolution, but not exactly.  
res <- c(0.5,0.5)
projr <- CRS("+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
r <- raster(x = e,resolution = res,crs = projr)

# read in one of the result files to view head and QAQC
read.csv(iric_results[1],skip = 2,header = T) %>%
  head()

# select the value (column) that will define the values of the processed iric result grid
val = "WaterSurfaceElevation"
# val = "Elevation"

# Create empty list to contain rasters
SMR_Q_l <- list()
# Create empty vector to contain modeled Qs
modeled_q <- c()

# for testing, set i to one of the iric_results objects, e.g.:
# i = "Result_0000.67.csv"

# For every iric output file (csv) in iric_results, 
for (i in iric_results) {
  # get iric output i
  df_i <-  read.csv(i,skip = 2,header = T)
  # set the discharge of iric output i as q_i [character]
  q_i <- sub("Result_","",sub(".csv","",i))
  # Transfer values from iric output i to the cells of raster r, but only when the cell is inundated (!= 0)
  # If multiple points from the iric output are within a cell of r, the mean points is used for the cell value.
  # Note that because the points in the iric output do not form a perfectly regular grid, where each point is the centroid of a 
  # square cell in r, the rasterize funciton may result in empty cells. These cells are filled in using bilinear interpolation in a 
  # following step.
  r_i <- rasterize(x = df_i[,c("X","Y")],y = r, field = ifelse(df_i[,"Depth"] == 0,NA,df_i[,val]),fun = mean)
  # set projection of raster i
  proj4string(r_i) <- projr
  # Resample raster i using bilinear interpolation to fill in those cells in r_i that did not have a cell value due to no point overlap
  r_rsmpl_i <- raster::projectRaster(from=r_i,to=r_i,method = 'bilinear')
  # update the name of the resampled raster to include the units of the discharge measure
  names(r_rsmpl_i) <- paste0("cms",q_i)
  # add resampled raster to the list
  SMR_Q_l[[i]] <- r_rsmpl_i
  # Create vector modeled qs
  modeled_q <- c(modeled_q,as.numeric(q_i))
}

plot(SMR_Q_l[[1]])

# If desired, write each raster in the list of rasters to drive
#sapply(SMR_Q_l,FUN = function(x) writeRaster(x,paste0("WSE_",names(x),".tif"),overwrite=FALSE))

# For creating dagwood sandwich of aerial extent of all inundation levels, first use the following reclassification
# that will change non NA values to the value of Q
raster_list_bi <- sapply(SMR_Q_l,function(x) {
  # create numeric object q, equal to the discharge of iric ouput i. 
  q <- as.numeric(sub("cms","",names(x)))
  # THIS IS AN ADDED LINE. It creates a raster reflecting the boundary cells of the inundating area; all other cells are set to NA
  x <- boundaries(x,asNA = TRUE)
  # set non NA cells to discharge, q
  x[!is.na(x)] <- as.numeric(q)
  # update the name of the new raster to include the units of the discharge measure
  names(x) <- paste0("cms",q) 
  return(x)
}
)

plot(raster_list_bi[[1]])

# If desired, write each raster in the list of rasters to drive
#sapply(raster_list_bi,FUN = function(x) writeRaster(x,paste0("Inund_",names(x),".tif")))


# NOW USE RASTER::MERGE to merge the rasters of sequential ascending discharge levels. 
# Priority in the case of overlaps is given to argument order (which is why sequential ascending qs is critical)
# See critical note below on sequential ascending q levels


# CRITICAL: The sequential layers of inundating discharge surfaces has an order that is inherited when it is read in from the directory.
# CRTIICAL: Therefore, the naming convention of the IRIC output csv files must be such that they are ordered correctly 
# CRITICAL: E.g., when viewing the csvs in File Explorer, when the folder is sorted by name, are they in the correct ascending order?
# CRITICAL: We achieve this using the following convention, IRICoutput_0000.01, IRICoutput_0001.5, IRICoutput_0010.0,IRICoutput_1010.0, etc
# There are many other ways to do this, but this is how the code we have written works for now. 

# Merge() requires adding each raster as an individual argument separated by a comma,
# we will use do.call() to run the merge function. 

# copy raster_list_bi to a new name, for work.
raster_list_bi_Qdagwood <- raster_list_bi

# set names of list object to NULL (do.call requirement)
names(raster_list_bi_Qdagwood) <- NULL
# add a object to the list for specifying each argument required in the raster::merge function, starting with filesname
# (see ?raster::merge)
raster_list_bi_Qdagwood$filename <- "Inundation_boundary_test.tif"
# argument overwrite set to false/true
raster_list_bi_Qdagwood$overwrite <- TRUE
# Now run the cover function on the raster_list_bi_docall object, which contains each dagwood layer and the required arguments
Q_Dagwood <- do.call(raster::merge,raster_list_bi_Qdagwood)

plot(Q_Dagwood)
zoom(Q_Dagwood)

# create points from Q_Dagwood that can be used for creating a TIN surface
QD_pt <- rasterToPoints(Q_Dagwood,spatial = TRUE)
st_write(st_as_sf(QD_pt),"SAL_complete_inun_map5_pt.shp")
plot(QD_pt)






