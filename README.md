# iRIC_output_processing 

This repository contains the R code developed for processing iric model output. 
It also contains example data for running the r script. Download and place all of 
the files into a directory titled "IRIC_Processing_in_R_v1" and update the setwd command
to run the example. 

V2 added on 2/7/19. This version creates a raster that reflects the combined boundaries of inundation of the ascending sequential discharges. TINing can then be used to interpolate the cells in between the edges of inundation to create a continuous surface of inundating discharge. 
