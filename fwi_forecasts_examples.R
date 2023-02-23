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

binning <- function(x, thresholds) {
  n_samples <- length(x)
  n_bins <- length(thresholds) + 1
  
  thresholds <- c(thresholds, max(x))
  result <- 1:n_bins
  lower_threshold <- min(x) - 1
  for (i in 1:n_bins) {
    result[i] <- sum(x > lower_threshold & x <= thresholds[i]) / n_samples
    lower_threshold <- thresholds[i]
  }
  
  dim(result) <- c(bin = n_bins)
  result
}

prob_thresholds <- c(1/3,2/3)
fix_thresholds <- c(5.2,11.2,21.3,38,50)
data <- readRDS('/work3/veronicatorralba/CMCC_fwi/CMCC_downscaled_calibrated_sdate_05_JJA_1993_2016.RDS')
exp_probs <- Apply(data$data_exp,target_dims='member',binning,thresholds=fix_thresholds)

#x11(width=14,height=10)
  col_sets <- list( c('#cbc9e2', '#9e9ac8', '#756bb1', '#54278f'),
                    c("#6BAED6FF", "#4292C6FF", "#2171B5FF", "#08519CFF"),
                    c("#A1D99B", "#74C476", "#41AB5D", "#238B45"),
                    c("#FFEDA0FF", "#FED976FF", "#FEB24CFF", "#FD8D3CFF"),
                    c("#FC4E2AFF", "#E31A1CFF", "#BD0026FF", "#800026FF"),
                    c("#FCC5C0", "#FA9FB5", "#F768A1", "#DD3497"))
  
  
  years <- 1993:2016
  for (i in 1:length(years)){
    cairo_ps(paste0('/work3/veronicatorralba/CMCC_fwi/EU_fwi_seasonal_',years[i],'.ps'),width = 14,height = 10)
    PlotMostLikelyQuantileMap(exp_probs[[1]][,i,,],lat=data$lat,lon=data$lon,colNA='white',
                              col_unknown_cat = 'grey',cols=col_sets,
                              bar_titles = c('Very low','Low','Moderate',
                                             'High','Very high','Extreme'))
    dev.off()
    cairo_ps(paste0('/work3/veronicatorralba/CMCC_fwi/SI_fwi_seasonal_',years[i],'.ps'),width = 14,height = 10)
    dat <- SelBox(exp_probs[[1]][,i,,],lat=as.vector(data$lat),
                  lon=as.vector(data$lon),region=c(13.5,19,36,43))
    PlotMostLikelyQuantileMap(dat$data,dat$lon,dat$lat,colNA='white',intxlon=2,intylat=1,
                              col_unknown_cat = 'grey',cols=col_sets,
                              bar_titles = c('Very low','Low','Moderate',
                                             'High','Very high','Extreme'))
    dev.off()
  }





