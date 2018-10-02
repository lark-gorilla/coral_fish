# 02/10/2018

# Download 

rm(list=ls())
# rerdap
library(rerddap)
library(raster)
library(ncdf4)
library(sf)

setwd("M:/coral_fish")

# We use rerddap package to make direct calls to the erddap server and get 
# gridded data via the OPeNDAP hyperslab protocol

#product code from here:
#https://coastwatch.pfeg.noaa.gov/erddap/griddap/index.html?page=1&itemsPerPage=1000

# pull all grid datasets
g1<-ed_datasets('grid')

# get summary of all 'Coral' datasets
g1[grep('Coral',g1$Title),]$Summary

# NOAA_DHW_5km is the 5 km product suite
# https://coralreefwatch.noaa.gov/satellite/index.php

info('NOAA_DHW_5km') # see variables
# rerddap package not working well.. try manual download

# use erddap data form to make query
# https://coastwatch.pfeg.noaa.gov/erddap/griddap/NOAA_DHW.html

### OK now read in sample sites and buffer to get reasonably sized bounding box
rmi_sites<-read.csv('data/RMI/RMI_~Sitetraits.csv', h=T)

rmi_sp<-st_as_sf(rmi_sites, coords = c("Long","Lat"), crs="+proj=longlat +datum=WGS84") # as sf object
rmi_buf<-st_buffer(rmi_sp, dist=0.5)# buffer it. dist in decimal degrees ~ 50 km
rmi_buf<-st_union(rmi_buf)# dissolve buffers

# get list of bboxes for download parameters 
bbox_list<-rmi_buf %>% st_cast('POLYGON') %>% lapply(st_bbox) 

#do.call('rbind', bbox_list)

bbox_list
# select number..
i=bbox_list[[1]]

x1=i$xmin
x2=i$xmax
y1=i$ymin
y2=i$ymax


CRW_DHW[(2004-09-01T12:00:00Z):1:(2014-09-01T12:00:00Z)][(6.56696):1:(7.72105)][(170.52708):1:(171.88225)],CRW_HOTSPOT[(2004-09-01T12:00:00Z):1:(2014-09-01T12:00:00Z)][(6.56696):1:(7.72105)][(170.52708):1:(171.88225)],CRW_SST[(2004-09-01T12:00:00Z):1:(2014-09-01T12:00:00Z)][(6.56696):1:(7.72105)][(170.52708):1:(171.88225)],CRW_SSTANOMALY[(2004-09-01T12:00:00Z):1:(2014-09-01T12:00:00Z)][(6.56696):1:(7.72105)][(170.52708):1:(171.88225)]

dl_key<-paste0('https://coastwatch.pfeg.noaa.gov/erddap/griddap/NOAA_DHW.nc?CRW_DHW[(2004-09-01T12:00:00Z):1:(2014-09-01T12:00:00Z)][(',y1,'):1:(',y2,')][(',x1,'):1:(',x2,')],
CRW_HOTSPOT[(2004-09-01T12:00:00Z):1:(2014-09-01T12:00:00Z)][(',y1,'):1:(',y2,')][(,',x1,'):1:(',x2,')],
CRW_SST[(2004-09-01T12:00:00Z):1:(2014-09-01T12:00:00Z)][(',y1,'):1:(',y2,')][(',x1,'):1:(',x2,')],
CRW_SSTANOMALY[(2004-09-01T12:00:00Z):1:(2014-09-01T12:00:00Z)][(',y1,'):1:(',y2,')][(',x1,'):1:(',x2,')]')

## Issue: only goes back until 2013 on ERDDAP servers grr

## Use these instad and write manual code
# ftp://ftp.star.nesdis.noaa.gov/pub/sod/mecb/crw/data/5km/


download.file(url=dl_key,
              destfile='C:/ocean_data/NOAA_coral_bleaching/RMI_1.nc')


https://coastwatch.pfeg.noaa.gov/coastwatch/CWBrowserWW360.jsp?get=gridData&dataSet=CGCssta&timePeriod=monthly&centeredTime=0001-09-16T00:00:00&minLon=165.5&maxLon=172.5&minLat=6&maxLat=12.5&fileType=.nc

# chl use erdVH3chla8day or erdMH1chla8day 
# sst erdAGssta8day or ncdcOisst2Agg
# ssh nrlHycomGLBu008e911S
# wind stress and upwelling erdQMstress3day
# sst anomaly ncdcOisst2Agg or erdAGtanm8day

# should do a extract3d on the blended sst anomaly
"nrlHycomGLBu008e911S", "erdVH3chla8day", "erdMH1chla8day",
"erdAGssta8day", "ncdcOisst2Agg", "erdAGtanm8day"
"erdQMstress8day" 

# we're gonna save and export as a .nc file as waaaaaay smaller than csv file size (50 vs 600 Mb)

ed_search(query = 'erdMH1chla16day', which = "grid")$info

info("erdVH3chlamday")

(res <- griddap("erdMH1chlamday",
                time = c('2014-02-01', '2014-04-30'),
                latitude = c(-10, -42),
                longitude = c(140, 170))) 