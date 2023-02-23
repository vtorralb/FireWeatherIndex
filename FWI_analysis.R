#---------------------------------------------------------
# Script to assess the fwi in the seasonal forecasts
# to be run in the servers (server74)
# Veronica Torralba
# February 2023
#---------------------------------------------------------

library(ncdf4)
library(RColorBrewer)
library(s2dverification)
library(multiApply)
library(easyVerification)
library(ClimProjDiags)
library(s2dv)
library(lubridate)
library(startR)
library(CSTools)
library(qmap)

years <- 1993:2016

path_obs <- '/work3/veronicatorralba/ERA5_fwi/Europe/'
path_exp <- '/work3/veronicatorralba/CMCC_fwi/'

# data array structure
fwi_seasonal_cmcc <- array(dim=c(75,45,92,40,24))
fwi_seasonal_obs <- array(dim=c(301,177,92,24))

# Loading the data (CMCC seasonal forecasts and ERA5) in JJA 
for (y in 1:length(years)){
  ncexp <- nc_open(paste0(path_exp,'cmcc-35_fwi_',years[y],'0501.nc'))
  exp <- ncvar_get(ncexp,'fwi')
  names(dim(exp)) <- c('lon','lat','time','member')
  fwi_seasonal_cmcc[,,,,y] <- Subset(exp,along='time',indices= 32:123)

  obs_sea <- NULL
  for (mon in 6:8){
    ncobs <- nc_open(paste0(path_obs,'ERA5_FWI_',years[y],sprintf('%02d',mon),'.nc'))
    obs <- ncvar_get(ncobs,'fwi')
    names(dim(obs)) <- c('lon','lat','time')
    obs_sea <- startR:::.MergeArrays(obs_sea,obs,along ='time')
  }
  fwi_seasonal_obs[,,,y] <- obs_sea
}
# Original coordinates
lon_obs <- as.vector(ncobs$dim$longitude$vals)
lat_obs <- as.vector(ncobs$dim$latitude$vals)
lon_exp <- as.vector(ncexp$dim$lon$vals)   
lat_exp <- as.vector(ncexp$dim$lat$vals) 

# Seasonal mean 
fwi_seasonal_cmcc <- Mean1Dim(fwi_seasonal_cmcc,3)
fwi_seasonal_obs <- Mean1Dim(fwi_seasonal_obs,3)
names(dim(fwi_seasonal_cmcc)) <- c('lon','lat','member','sdate')
names(dim(fwi_seasonal_obs)) <- c('lon','lat','sdate')

# Downscaling (bilinear interpolation)
fwi_seasonal_down <- CDORemap(data_array = fwi_seasonal_cmcc,
                              lats=lat_exp,lons=lon_exp,
                              method = 'conservative',
                              grid = paste0(path_obs,'ERA5_FWI_199306.nc'))

# Selecting the observations in the common grid
#obs1 <- Subset(fwi_seasonal_obs, along = 'lat', indices=which(lat_obs %in% fwi_seasonal_down$lat))

fwi_seasonal_bc <- Calibration(exp=fwi_seasonal_down$data_array,obs=fwi_seasonal_obs,
                               cal.method = 'mse_min',
                               eval.method = 'in-sample')

data <- list(data_exp=fwi_seasonal_bc,
               data_obs=fwi_seasonal_obs,
             lon=fwi_seasonal_down$lon,
             lat=fwi_seasonal_down$lat)


for (y in 1:length(years)){
  fwi_bc <- Subset(data$data_exp,along='sdate',indices=1,drop=TRUE)
  fwi_bc <- aperm(fwi_bc,c(2,3,1))
  metadata <- list(fwi_bc = list(units = 'Unitless'))
  attr(fwi_bc, 'variables') <- metadata
  lon <- data$lon
  dim(lon) <- length(lon)
  metadata <- list(lon = list(units = 'degrees_east'))
  attr(lon, 'variables') <- metadata
  names(dim(lon)) <- 'lon'
  lat <- data$lat
  dim(lat) <- length(lat)
  metadata <- list(lat = list(units = 'degrees_north'))
  attr(lat, 'variables') <- metadata
  names(dim(lat)) <- 'lat'
  member <- 1:40
  dim(member) <- length(member)
  metadata <- list(member = list(units = ''))
  attr(member, 'variables') <- metadata
  names(dim(member)) <- 'member'
  ArrayToNetCDF(list(fwi_bc, lon, lat, member), paste0(path_exp,'fwi_down_bc/cmcc-35_fwi_bc_',years[y],'0501_JJA.nc'))
}
# 
# obs_low <- CDORemap(data_array = fwi_seasonal_obs,
#                               lats=lat_obs,lons=lon_obs,
#                               method = 'conservative',
#                               grid = paste0(path_exp,'cmcc-35_fwi_19930501.nc'))
# 
# # Basic tests
# #-------------------------
# clim_obs <- Mean1Dim(fwi_seasonal_obs,3)
# clim_exp_or <- MeanListDim(fwi_seasonal_cmcc,c(3,4))
# clim_exp_down <- MeanListDim(fwi_seasonal_down$data,c(3,4))
# clim_exp_down_bc <- MeanListDim(fwi_seasonal_bc,c(3,4))
# 
# 
# exp <- aperm(fwi_seasonal_bc,c(3,4,2,1))
# em <- MeanDims(exp,dims='member')
# cor_map <- s2dv::Corr(exp =em ,obs=fwi_seasonal_obs,time_dim='sdate',dat_dim = NULL)
# rpss <- s2dv::RPSS(exp=exp,obs=fwi_seasonal_obs,time_dim='sdate',memb_dim='member')
# 
# outputs <- list(clim_obs=clim_obs,clim_exp_or=clim_exp_or,
#                 clim_exp_down=clim_exp_down,
#                 clim_exp_down_bc=clim_exp_down_bc,
#                 corr_bc=cor_map,
#                 rpss_bc=rpss,
#                 lat_high=lat_obs,
#                 lon_high=lon_obs,
#                 lat_low=lat_exp,
#                 lon_low=lon_exp)
#         
# saveRDS(outputs,'/work3/veronicatorralba/CMCC_fwi/fwi_data_4_plots.RDS')


