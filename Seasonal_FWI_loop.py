#!/usr/bin/env python
# coding: utf-8

# In[1]:


import xarray as xr
import pandas as pd
from scipy import signal
import os
import datetime
import time
import shutil
import sys
import numpy as np
import numpy.ma as ma
from dateutil.relativedelta import relativedelta
import math
from glob import glob
from netCDF4 import num2date, date2num
import calendar
from netCDF4 import num2date, date2num, Dataset
from numpy import dtype
import xclim as xc
import netCDF4 as nc


# In[2]:


def psurface(temp,slp,elev):
    g=9.807
    elevation=elev/g
    R=8.314
    mair=0.029
    lapse=6.5/1000
    tsurf = temp + (lapse * elevation)
    psurf = slp * math.exp((-elevation * g * mair)/(tsurf * R))
    return(psurf)

Psurf=np.vectorize(psurface)


def RH_ecmwf(x,y,p):
    a1 = 611.21
    a3_aux=22.587
    a3 = 17.502
    Rdry = 287.0597
    Rvap = 461.5250
    a4 = 32.19
    a4_aux= -0.7
    T0 = 273.16
    Tice = 250.16
    R = Rdry/Rvap
    E = a1 * math.exp(a3*((y - T0)/(y-a4)))
    Esat1 = a1 * math.exp(a3*((x - T0)/(x-a4)))
    Esat2 = a1 * math.exp(a3_aux*((x - T0)/(x-a4_aux)))
    if (x<=Tice):
        alpha=0
    elif (Tice < x < T0):
        alpha=((x-Tice)/(T0-Tice))
    elif (x>= T0):
        alpha = 1
    Esat= (alpha* Esat1) + ((1-alpha) * Esat2)
    q = R * (E/(p -((1-R)*E)))
    RH = (p * (q*(1/R)))/(Esat * (1 +(q*((1/R)-1))))
    return(RH *100)

RH_ECMWF=np.vectorize(RH_ecmwf)

def load_eu(path1):
    min_lon = -15
    min_lat = 25
    max_lon = 60
    max_lat = 70
    ds = xr.load_dataset(path1, engine="cfgrib")
    ds.coords['longitude']=(ds.coords['longitude'] + 180) % 360 - 180
    ds = ds.sortby(ds.longitude)
    out = ds.sel(latitude=slice(max_lat,min_lat),longitude=slice(min_lon,max_lon))
    return out


# In[ ]:


# loading orography (needed for the computation of the surface pressure from the SLP)
orog = xr.load_dataset('/data/csp/vt17420/CLINT_proj/ERA5/ERA5_masks/oro_Europe_1deg.nc')
oro_eu = orog[['z']].to_array()[0,:,:]

path='/work/asc/dias/Datastore/Seasonal/'
sys='cmcc-35'
nmemb=4

#year=2008
for year in range(1993,1994):
    for i in range(0,nmemb):
        t2m_file='seasonal_original_single_levels-2m_temperature-'+sys+'-'+str(year)+'.grib'
        d2m_file='seasonal_original_single_levels-2m_dewpoint_temperature-'+sys+'-'+str(year)+'.grib'
        psl_file='seasonal_original_single_levels-mean_sea_level_pressure-'+sys+'-'+str(year)+'.grib'
        pr_file='seasonal_original_single_levels-total_precipitation-'+sys+'-'+str(year)+'-05-01.grib'
        uas_file='seasonal_original_single_levels-10m_u_component_of_wind-'+sys+'-'+str(year)+'-05-01.grib'
        vas_file='seasonal_original_single_levels-10m_v_component_of_wind-'+sys+'-'+str(year)+'-05-01.grib'
        t2m=load_eu(path+t2m_file).to_array()[0,i,0,:,:,:]
        d2m=load_eu(path+d2m_file).to_array()[0,i,0,:,:,:]
        psl=load_eu(path+psl_file).to_array()[0,i,0,:,:,:]
        uas=load_eu(path+uas_file).to_array()[0,i,:,:,:]
        vas=load_eu(path+vas_file).to_array()[0,i,:,:,:]
        pr=load_eu(path+pr_file).to_array()[0,i,:,:,:]

        ntime=psl.shape[0]
        nlat=psl.shape[1]
        nlon=psl.shape[2]
        elevat_aux=np.zeros([nlat,nlon,ntime])
        for n in range(ntime):
                elevat_aux[:,:,n]=oro_eu[0,0:nlat,0:nlon]

        # computation of the surface pressure from the SLP
        psl2ps = xr.apply_ufunc(
            Psurf,
            t2m,
            psl,
            elevat_aux,
            dask="parallelized",
            input_core_dims=[['step'], ['step'],['step']],
            output_core_dims=[['step']],
            output_dtypes=[float])

        # computation of the relative humidity (temeprature, dewpoint and surface pressure)
        RH = xr.apply_ufunc(
                RH_ECMWF,
                t2m,
                d2m,
                psl2ps,
                dask="parallelized",
                input_core_dims=[['step'], ['step'],['step']],
                output_core_dims=[['step']],
                output_dtypes=[float])


        # Daily values
        rh=RH.resample(step='1D').mean()
        sfcwind= 3.6 * np.sqrt(uas**2 + vas**2)
        sfcWind=sfcwind.resample(step='1D').mean()
        # units adjustment
        tmax=t2m.resample(step='1D').max()-273.15
        prec=pr*1000

        # Coordinates
        lats=t2m.coords['latitude']
        lons=t2m.coords['longitude']
        dates = pd.date_range(start=str(year)+"-05-01", end=str(year)+"-10-31", freq="D")

        # DATA ARRAY OBJECTS NEED TO BE CREATED
        tmax_ds = xr.DataArray(data=tmax.values,
         dims=["time","latitude","longitude"],
         coords=({'latitude':lats,
                 'longitude':lons,
                 'time':dates}),
         attrs=dict(units="C"))

        sfcWind_ds = xr.DataArray(data=sfcWind.values,
         dims=["time","latitude","longitude"],
         coords=({'latitude':lats,
                 'longitude':lons,
                 'time':dates}),
         attrs=dict(units="km/h"))

        rh_ds = xr.DataArray(data=rh.values,
         dims=["latitude","longitude","time"],
         coords=({'latitude':lats,
                 'longitude':lons,
                 'time':dates}),
         attrs=dict(units="pct"))

        prec_ds = xr.DataArray(data=prec.values,
         dims=["time","latitude","longitude"],
         coords=({'latitude':lats,
                 'longitude':lons,
                 'time':dates}),
         attrs=dict(units="mm/day"))


        # Computation of the Fire Weather Index
        fwi_system = xc.indices.cffwis_indices(tas=tmax_ds,
                                               pr=prec_ds,
                                               sfcWind=sfcWind_ds,
                                               hurs=rh_ds, 
                                               lat=lats)

        # OUTPUT VARIABLES
        varn_short = ["dc", "dmc", "ffmc", "isi", "bui", "fwi"] 
        varn_long = ["Drought Code",                                 
                     "Duff Moisture Code",
                     "Fine Fuel Moisture Code",
                     "Initial Spread Idex",
                     "Buildup Index",
                     "Fire Weather Index"]

        for ivar, var in enumerate(varn_short):
                        print("    {}...".format(varn_short[ivar]))
                        dirout = "/work/csp/vt17420/FWI_seasonal/"
                        fout = dirout+sys+'_'+var+'_'+str(year)+'0501'+'_memb_'+str(i)+'.nc'
                        fout_nc = nc.Dataset(fout, mode="w", format="NETCDF4_CLASSIC")

                        # Definición de las dimensiones
                        lat_dim = fout_nc.createDimension('lat', len(lats))
                        lon_dim = fout_nc.createDimension('lon', len(lons))
                        time_dim = fout_nc.createDimension('time', None)

                        # Definición de las variables
                        lat_def = fout_nc.createVariable('lat', np.float32, ('lat',))
                        lat_def.units = 'degrees_north'
                        lat_def.long_name = 'latitude'
                        lon_def = fout_nc.createVariable('lon', np.float32, ('lon',))
                        lon_def.units = 'degrees_east'
                        lon_def.long_name = 'longitude'
                        time_def = fout_nc.createVariable('time', np.float64, ('time',))
                        time_def.long_name = 'time'
                        time_def.calendar = '365_day'
                        var_def = fout_nc.createVariable(var,np.float64,('time','lat','lon'))
                        var_def.units = 'Unitless'
                        var_def.standard_name = varn_long[ivar]

                        # Almacenamiento de datos
                        lat_def[:] = lats
                        lon_def[:] = lons
                        time_def[:] = dates
                        # En el almacenamiento del FWI se intercambian los ejes del array: (lat,lon,time) ---> (time,lat,lon)
                        var_def[:] = np.moveaxis(fwi_system[ivar].values, 2, 0)
                        fout_nc.close()


# In[ ]:




