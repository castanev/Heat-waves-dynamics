# This code is for the collection of NCEP data and construction of netcdf files.

# =================================== NCEP ================================================================
# Originaly, the latitudes in NCEP files are in descending order, the outputs here are in ascending order.
 
#  - Surface temperature (at 1000 hPa)
#  - u component of wind at 300 hPa
#  - v component of wind at 300 hPa


import wget
from netCDF4 import Dataset
import pandas as pd
import numpy as np
import datetime as dt
import os


# INPUTS TO CHANGE ===========================================================================================================================
path_data_t = f'/scratch/brown/castanev/NCEP/NCEP_data_T/'
path_data_v = f'/scratch/brown/castanev/NCEP/NCEP_data_v/'
path_data_u = f'/scratch/brown/castanev/NCEP/NCEP_data_u/'

# ============================================================================================================================================

name = 'NCEP'
path_figures = f'{os.path.abspath(os.getcwd())}/{name}/Figures/'
path_outputs = f'{os.path.abspath(os.getcwd())}/{name}/Data/'


# T at the surface ============================================================================
# Saving T at approximately 1000 hPa (.nc) ====================================================
list_data = np.sort(os.listdir(f'{path_data_t}'))
years = pd.date_range(start=dt.date(1948, 1,1), end=dt.date(2023, 1, 1), freq='A')

# Downloading the data of surface temperature
for year in years.strftime('%Y'):
    print(f'ftp://ftp2.psl.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/surface/air.sig995.{year}.nc')
    filename = wget.download(f'ftp://ftp2.psl.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/surface/air.sig995.{year}.nc', out = f'{path_data_t}')


ncfile = Dataset(f'{path_data_t}air.sig995.1948.nc')
lats = np.array(ncfile['lat'])
lons = np.array(ncfile['lon'])

new_lats = lats[::-1]

# Concatenating the data into one array, considering only u in 300hPa
contt = 0
for i in list_data:
    if 'air.sig995.' in i:
        contt = contt + 1
        datat_i_nc = Dataset(f'{path_data_t}{i}')  #(t, lat, lon)
        datat_i = np.array(datat_i_nc['air'])[:, ::-1, :]  #t at surface
        timei_t = np.array(datat_i_nc['time']) #Daily
        if contt == 1: t = datat_i; dates_t = timei_t
        else: t = np.concatenate((t, datat_i), axis = 0); dates_t = np.concatenate((dates_t, timei_t))
        print(contt)
        print(f'{path_data_t}{i}')


nc_name = f't1000_{name}.nc'
ncfile = Dataset(f'{path_outputs}{nc_name}', 'w')

ncfile.createDimension('lat', len(new_lats))
ncfile.createDimension('lon', len(lons))
ncfile.createDimension('time', len(dates_t))

var_lats = ncfile.createVariable('lat', 'f', ('lat'))
var_lons = ncfile.createVariable('lon', 'f', ('lon'))
var_time = ncfile.createVariable('time', 'f', ('time'))

var_lats[:] = new_lats
var_lons[:] = lons
var_time[:] = dates_t

vwnd = ncfile.createVariable('t1000', 'f', ('time', 'lat', 'lon'))
vwnd[:, :, :] = t[:, :, :]
ncfile.close()




# u wind at 300 hPa ===========================================================================
#==============================================================================================
list_data = np.sort(os.listdir(f'{path_data_u}'))
nc_name = f'uwnd.1948.nc'
ncfile = Dataset(f'{path_data_u}/{nc_name}')
lev = np.array(ncfile['level'])  
lats = np.array(ncfile['lat'])  
lons = np.array(ncfile['lon'])  
pos_lev = np.where(abs(lev-300) == np.min(abs(lev-300)))

new_lats = lats[::-1]

variable = 'u_300'
contt = 0
for i in list_data:
    if 'uwnd' in i:
        contt = contt + 1
        datat_i_nc = Dataset(f'{path_data_u}{i}')  #(t, level, lat, lon)
        datat_i = np.array(datat_i_nc['uwnd'])[:, pos_lev[0][0], :, :]
        datat_i = datat_i[:,::-1,:]

        timei = np.array(datat_i_nc['time'])
        df = pd.DataFrame(np.reshape(datat_i, (datat_i.shape[0],datat_i.shape[1]*datat_i.shape[2])), index = pd.DatetimeIndex([dt.datetime(1800,1,1) + dt.timedelta(hours = int(timei[i])) for i in range(len(timei))]))

        if contt == 1: dataa = datat_i; dates = df.index
        else: dataa = np.concatenate((dataa, datat_i), axis = 0); dates = np.concatenate((dates, df.index))


# Saving the data
# Lats in ascending order 
nc_name = f'{variable}_{name}.nc'
ncfile = Dataset(f'{path_outputs}{nc_name}', 'w')

ncfile.createDimension('lat', len(new_lats))
ncfile.createDimension('lon', len(lons))
ncfile.createDimension('time', len(dates))

var_lats = ncfile.createVariable('lat', 'f', ('lat'))
var_lons = ncfile.createVariable('lon', 'f', ('lon'))
var_time = ncfile.createVariable('time', 'f', ('time'))

var_lats[:] = new_lats
var_lons[:] = lons
var_time[:] = dates

varr = ncfile.createVariable(variable, 'f', ('time', 'lat', 'lon'))
varr[:, :, :] = dataa[:, :, :]
ncfile.close()




# v wind at 300 hPa ===========================================================================
#==============================================================================================
list_data = np.sort(os.listdir(f'{path_data_v}'))
nc_name = f'vwnd.1948.nc'
ncfile = Dataset(f'{path_data_v}/{nc_name}')
lev = np.array(ncfile['level'])  
lats = np.array(ncfile['lat'])  
lons = np.array(ncfile['lon'])  
pos_lev = np.where(abs(lev-300) == np.min(abs(lev-300)))

new_lats = lats[::-1]

variable = 'v_300'
contt = 0
for i in list_data:
    if 'vwnd' in i:
        contt = contt + 1
        datat_i_nc = Dataset(f'{path_data_v}{i}')  #(t, level, lat, lon)
        datat_i = np.array(datat_i_nc['vwnd'])[:, pos_lev[0][0], :, :]
        datat_i = datat_i[:,::-1,:]

        timei = np.array(datat_i_nc['time'])
        df = pd.DataFrame(np.reshape(datat_i, (datat_i.shape[0],datat_i.shape[1]*datat_i.shape[2])), index = pd.DatetimeIndex([dt.datetime(1800,1,1) + dt.timedelta(hours = int(timei[i])) for i in range(len(timei))]))

        if contt == 1: dataa = datat_i; dates = df.index
        else: dataa = np.concatenate((dataa, datat_i), axis = 0); dates = np.concatenate((dates, df.index))
        print(i)



# Saving the data
# Lats in ascending order 
nc_name = f'{variable}_{name}.nc'
ncfile = Dataset(f'{path_outputs}{nc_name}', 'w')

ncfile.createDimension('lat', len(new_lats))
ncfile.createDimension('lon', len(lons))
ncfile.createDimension('time', len(dates))

var_lats = ncfile.createVariable('lat', 'f', ('lat'))
var_lons = ncfile.createVariable('lon', 'f', ('lon'))
var_time = ncfile.createVariable('time', 'f', ('time'))

var_lats[:] = new_lats
var_lons[:] = lons
var_time[:] = dates

varr = ncfile.createVariable(variable, 'f', ('time', 'lat', 'lon'))
varr[:, :, :] = dataa[:, :, :]
ncfile.close()