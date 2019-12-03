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
# ok so there is MIN, MAX and MEAN products per month for sst and ssta,
# only MAX availiable for dha, hs and baa
# sst - sea surface temp
# ssta - sea surface temp anomoly
# hs - bleaching hotspot
# dhw - degree heating weeks
# baa - bleaching alert area

# remember these are global rasters so we'll grab from 09/2017 (Australia - most recent) back to 
# 10 years prior to 02/2010 (Australia - oldest), so 02/2000. This also covers RMI and Japan (2014, 15/16)

url1='ftp://ftp.star.nesdis.noaa.gov/pub/sod/mecb/crw/data/5km/v3.1/nc/v1.0/monthly/'

## NOTE: ncdf4 cannot read downloaded .nc unless mode='wb' !!!
# https://stackoverflow.com/questions/50048989/downloading-netcdf-files-with-r-manually-works-download-file-produces-error

for (i in 2000:2017)
{
  for(j in c('01','02','03','04','05','06','07','08','09','10','11','12'))
    {
    dl=paste0(url1,i,'/ct5km_sst-mean_v3.1_',i, j, '.nc')
  
    download.file(url=dl,
                destfile=paste0('C:/ocean_data/NOAA_coral_bleaching/sst/',
                               'ct5km_sst-mean_v3.1_',i, j, '.nc'), mode='wb')
    
    dl2=paste0(url1,i,'/ct5km_ssta-mean_v3.1_',i, j, '.nc')
    
    download.file(url=dl2,
                  destfile=paste0('C:/ocean_data/NOAA_coral_bleaching/ssta/',
                                  'ct5km_ssta-mean_v3.1_',i, j, '.nc'), mode='wb')
    
    dl3=paste0(url1,i,'/ct5km_dhw-max_v3.1_',i, j, '.nc')
    
    download.file(url=dl3,
                  destfile=paste0('C:/ocean_data/NOAA_coral_bleaching/dhw/',
                                  'ct5km_dhw-max_v3.1_',i, j, '.nc'), mode='wb')
    print(paste(i,j)) 
    }  
 
}

# slightly different loop for daily data

# read in site data

jpn_sites<-read.csv('M:/coral_fish/data/Japan/JP2015_16_waypoints.csv', h=T)
rmi_sites<-read.csv('M:/coral_fish/data/RMI/RMI_~Sitetraits.csv', h=T)
aus_sites<-read.csv('M:/coral_fish/data/Australia/Australia_SitesMar2010toAug2017.csv', h=T)

sites<-rbind(data.frame(site=jpn_sites$Site, region='JPN', date=jpn_sites$Date, lat=jpn_sites$lat, long=jpn_sites$lon),
             data.frame(site=rmi_sites$SiteID, region='RMI', date=rmi_sites$Date, lat=rmi_sites$Lat, long=rmi_sites$Long),
             data.frame(site=aus_sites$Site.name, region='AUS', date=aus_sites$Date, lat=aus_sites$Lat, long=aus_sites$Long))

sites[108,]$lat= -28.61087
sites[108,]$long=	153.62809

sites<-na.omit(sites) # remove NA rows

sites<-sites[-120,] # remove blank

sites<-sites[-which(duplicated(sites$site)),] # remove site duplicates


#setup master dataframe
r1<-raster('C:/ocean_data/template_sst.nc')
ext2<-extract(r1, cbind(sites$long, sites$lat), buffer=50000, cellnumbers=T)

out<-NULL
for(i in 1:nrow(sites))
{
  df<-data.frame(sites[i,], coordinates(r1)[ext2[[i]][,1],])
  out<-rbind(out, df)
  print(i)
}


#https://stackoverflow.com/questions/52182635/r-reading-geotiff-data-straight-from-web-url-httrget-raw-content

for (i in 2000:2017)
{
  # for sst
  url1<-paste0('ftp://ftp.star.nesdis.noaa.gov/pub/sod/mecb/crw/data/5km/v3.1/nc/v0.1/daily/sst/', i, '/')
  
  filez<-unlist(strsplit(getURL(url1, dirlistonly=T),'\r\n'))
  filez<-filez[order(substr(filez, 15,22))]
  
  yr_out<-out
  for(j in filez)
  {
    
    download.file(url=paste0(url1,j),
                  destfile='C:/ocean_data/temp_sst.nc', mode='wb')
    
    
    r1<-raster('C:/ocean_data/temp_sst.nc')
    
    ext1<-raster::extract(r1, cbind(sites$long, sites$lat), buffer=50000)
    
    yr_out<-cbind(yr_out, unlist(ext1))
    names(yr_out)[names(yr_out)=='unlist(ext1)']<-paste0('s', substr(j, 15,22))
    
    file.remove('C:/ocean_data/temp_sst.nc')
  
    
    # delete temp file
    
 
    
    #dl2=paste0(url1,i,'/ct5km_ssta-mean_v3.1_',i, j, '.nc')
    
    #download.file(url=dl2,
    #              destfile=paste0('C:/ocean_data/NOAA_coral_bleaching/ssta/',
    #                              'ct5km_ssta-mean_v3.1_',i, j, '.nc'), mode='wb')
    
    #dl3=paste0(url1,i,'/ct5km_dhw-max_v3.1_',i, j, '.nc')
    
    #download.file(url=dl3,
    #              destfile=paste0('C:/ocean_data/NOAA_coral_bleaching/dhw/',
    #                              'ct5km_dhw-max_v3.1_',i, j, '.nc'), mode='wb')
    print(paste(i,j)) 
    
    
  }
  write.csv(yr_out, paste0('C:/ocean_data/NOAA_coral_bleaching/daily_sst/sst_', i, '.csv'), quote=F, row.names=F)
  
}



## Now extract chl-a data from ERDDAP servers. rerrdap not working

library(rerrdap)

info('erdMH1chla8day')

sites_sp<-st_as_sf(sites, coords = c("long","lat"), crs="+proj=longlat +datum=WGS84") # as sf object
sites_buf<-st_buffer(sites_sp, dist=0.5)# buffer it. dist in decimal degrees ~ 50 km

for(i in 1: nrow(sites_buf))
{
  beeb<-st_bbox(sites_buf[i,])
  
  # MODIS chlorophyll
  #urly<-paste0('https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMH1chla8day.nc?chlorophyll[(2003-01-29):1:(2017-12-31)][(',
  #             beeb$ymax,'):1:(',beeb$ymin,')][(',beeb$xmin,'):1:(',beeb$xmax,')]')
  
  #download.file(url=urly,
  #              destfile=paste0('C:/ocean_data/AQUA_MODIS_chla/',sites_buf[i,]$site,'.nc'), mode='wb')
  
  # PML processed CHL and kd490 data 
  urly<-paste0('https://coastwatch.pfeg.noaa.gov/erddap/griddap/pmlEsaOCCCI31KD490Weekly.nc?chlor_a[(2000-01-01):1:(2017-12-31)][(',
               beeb$ymax,'):1:(',beeb$ymin,')][(',beeb$xmin,'):1:(',beeb$xmax,')]')
  
  
  err<-try(download.file(url=urly,
                destfile=paste0('C:/ocean_data/PML_chla/',sites_buf[i,]$site,'.nc'), mode='wb'))
  
  while(class(err) == "try-error"){
    err<-try(download.file(url=urly,
                           destfile=paste0('C:/ocean_data/PML_chla/',sites_buf[i,]$site,'.nc'), mode='wb'))
                            } 
  
  urly<-paste0('https://coastwatch.pfeg.noaa.gov/erddap/griddap/pmlEsaOCCCI31KD490Weekly.nc?kd_490[(2000-01-01):1:(2017-12-31)][(',
               beeb$ymax,'):1:(',beeb$ymin,')][(',beeb$xmin,'):1:(',beeb$xmax,')]')
  
  
  err<-try(download.file(url=urly,
                destfile=paste0('C:/ocean_data/PML_kd490/',sites_buf[i,]$site,'.nc'), mode='wb'))
  
  while(class(err) == "try-error"){
    err<-try(download.file(url=urly,
                          destfile=paste0('C:/ocean_data/PML_kd490/',sites_buf[i,]$site,'.nc'), mode='wb'))
                          } 
 print(i)          
}

# get list of bboxes for download parameters 
bbox_list<-sites_buf %>% st_cast('POLYGON') %>% lapply(st_bbox) 




r2<-raster('C:/ocean_data/AQUA_MODIS_chla/JP1.nc', band=2)

as.Date(as.POSIXlt(getZ(r2), origin="1970-01-01", "GMT"))


https://coastwatch.pfeg.noaa.gov/erddap/griddap/pmlEsaOCCCI31KD490Weekly.nc?chlor_a[(2000-01-01):1:(2017-12-31)][(24.84577202):1:(23.84577202)][(123.189956):1:(124.189956)],kd_490[(2000-01-01):1:(2017-12-31)][(24.84577202):1:(23.84577202)][(123.189956):1:(124.189956)]
https://coastwatch.pfeg.noaa.gov/erddap/griddap/pmlEsaOCCCI31KD490Weekly.nc?chlor_a[(2000-01-01):1:(2017-12-31)][(24.84577202):1:(23.84577202)][(123.189956):1:(124.189956)],kd_490[(2000-01-01):1:(2017-12-31)][(24.84577202):1:(23.84577202)][(123.189956):1:(124.189956)]


# Get doldrums data
v1<-expand.grid(2014:2016, 01:12, 01:31)

for (i in 1:nrow(v1))
{
  tstamp<-paste0(v1[i,]$Var1, sprintf("%02d", v1[i,]$Var2), sprintf("%02d", v1[i,]$Var3))
  
  urly=paste0('https://coralreefwatch.noaa.gov/satellite/doldrums_v2/data/hdf/archive/wind.doldrum.v3.field.25km.seawind.',
              tstamp, '.hdf')
  
  err<-try(download.file(url=urly,
                         destfile=paste0('C:/ocean_data/doldrums/',tstamp,'.hdf'), mode='wb'))
  print(i)
  
  if(class(err) == "try-error"){print(tstamp);next}
  
}

#reading in
#library(gdalUtils); library(raster)
#gdalinfo('C:/ocean_data/doldrums/20160303.hdf')
#https://stackoverflow.com/questions/48599875/converting-hdf-to-georeferenced-file-geotiff-shapefile
