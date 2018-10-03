# 02/10/2018

# Download 

rm(list=ls())
# rerdap
library(rerddap)
library(raster)
library(ncdf4)
library(sf)
library(RCurl)

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

dl_key<-paste0('https://coastwatch.pfeg.noaa.gov/erddap/griddap/NOAA_DHW.nc?CRW_DHW[(2004-09-01T12:00:00Z):1:(2014-09-01T12:00:00Z)][(',y1,'):1:(',y2,')][(',x1,'):1:(',x2,')],
CRW_HOTSPOT[(2004-09-01T12:00:00Z):1:(2014-09-01T12:00:00Z)][(',y1,'):1:(',y2,')][(,',x1,'):1:(',x2,')],
CRW_SST[(2004-09-01T12:00:00Z):1:(2014-09-01T12:00:00Z)][(',y1,'):1:(',y2,')][(',x1,'):1:(',x2,')],
CRW_SSTANOMALY[(2004-09-01T12:00:00Z):1:(2014-09-01T12:00:00Z)][(',y1,'):1:(',y2,')][(',x1,'):1:(',x2,')]')

download.file(url=dl_key,
              destfile='C:/ocean_data/NOAA_coral_bleaching/RMI_1.nc')

## Issue: only goes back until 2013 on ERDDAP servers - need further back in time.

## Solution: Use FTP from https://coralreefwatch.noaa.gov to gain access to monthly product:
# ftp://ftp.star.nesdis.noaa.gov/pub/sod/mecb/crw/data/5km/

#get url files for monthly products
u1<-getURL('ftp://ftp.star.nesdis.noaa.gov/pub/sod/mecb/crw/data/5km/v3.1/nc/v1.0/monthly/', dirlistonly=T)

# have alook at years
strsplit(u1, '\r\n')

#have a look at annual contents
strsplit(getURL('ftp://ftp.star.nesdis.noaa.gov/pub/sod/mecb/crw/data/5km/v3.1/nc/v1.0/monthly/1990/', dirlistonly=T),'\r\n')

# have alook at years
strsplit(u1, '\r\n')
#
# ok so there is MIN, MAX and MEAN products per month and
# sst - sea surface temp
# dha - degree heating weeks
# hs - heat stress
# ssta - sea surface temp anomoly

# remember these are global rasters so we'll grab from 09/2017 (Australia - most recent) back to 
# 10 years prior to 02/2010 (Australia - oldest), so 02/2000. This also covers RMI and Japan (2014, 15/16)

url1='ftp://ftp.star.nesdis.noaa.gov/pub/sod/mecb/crw/data/5km/v3.1/nc/v1.0/monthly/'

for (i in 2000:2017)
{
  for(j in c('01','02','03','04','05','06','07','08','09','10','11','12'))
    {
    dl=paste0(url1,i,'/ct5km_sst-mean_v3.1_',i, j, '.nc')
  
    download.file(url=dl,
                destfile=paste0('C:/ocean_data/NOAA_coral_bleaching/sst/',
                                'ct5km_sst-mean_v3.1_',i, j, '.nc'))
    }  
print(paste(i,j))  
}


