# This code works with NCEP data. It generates plots of: 
# - composites for heatwaves events considering: surface air temperature, streamfunction anomalies and the RWP envelope
# - Hovmoller diagram considering the same variables 
# NOTE: 
# The RWP envelope is calculated following Fragkoulidis, 2020. (Hilbert transform using a filter band). In this case, 
# we use a filter to consider only wavenumbers between 4 and 8

from Functions import *
import pandas as pd
import datetime as dt
from netCDF4 import Dataset
import scipy as scp
from dateutil.relativedelta import relativedelta
import matplotlib.pyplot as plt
import numpy as np
import scipy as scp
from cartopy import crs
import cartopy
import matplotlib.ticker as mticker
import matplotlib.ticker as ticker
import matplotlib.colors as mcolors
import os


# INPUTS TO CHANGE ===========================================================================================================================
# name: It can be NCEP or the name of the experiment (output of the model)
# var: the variable name in the nc file. 't1000' for NCEP, 'temp' for the dry core GCM 
# lat_minHW, lat_maxHW, lon_minHW, lon_minHW: to define the area where the events will be detected
# delete_days is the number of days to cut (only for output models). 

# name = 'NCEP'  
# nc_name_t = 't1000_NCEP.nc'
# var_t = 't1000'
# nc_name_v = 'v_300_NCEP.nc'
# var_v = 'v_300'
# nc_name_sf = 'sf_vp_300_1948-2022.nc' 
# var_sf = 'SF'
# lat_minHW = 25; lat_maxHW = 50; lon_minHW = 235; lon_maxHW = 290; midlat=45
# topography = True
# seasons = True 
# delete_days = 0
# path_data = f'/scratch/brown/castanev/{name}/'

name = 'exp1_Held_Suarez'  
nc_name_t = 't.atmos_daily.nc'
var_t = 'temp'
nc_name_v = 'v300.atmos_daily.nc'
var_v = 'v'
nc_name_sf = 'streamfunction_daily.nc' 
var_sf = 'SF'
lat_minHW = 30; lat_maxHW = 50; lon_minHW = 245; lon_maxHW = 275; midlat=45
topography = False
seasons = False 
delete_days = 1000
path_data = f'/scratch/brown/castanev/DryCore_Wu/output/{name}/post_processed/output/'

# ============================================================================================================================================
# path_data = f'{os.path.abspath(os.getcwd())}/{name}/Data/'
# path_figures = f'{os.path.abspath(os.getcwd())}/{name}/Figures/'
# path_outputs = f'{os.path.abspath(os.getcwd())}/{name}/'


path_figures = f'/home/castanev/Heat-waves-dynamics/{name}/Figures/'
path_outputs = f'/home/castanev/Heat-waves-dynamics/{name}/'


# Lats and lons 
ncfile = Dataset(f'{path_data}{nc_name_v}')
lats = np.array(ncfile['lat'])
lons = np.array(ncfile['lon'])  
timei = np.array(ncfile['time'])  

if   name == 'NCEP': time = np.array([dt.datetime(1800,1,1) + dt.timedelta(hours = int(timei[i])) for i in range(len(timei))])
elif name in ['CAM', 'LENS']:  
    timei = [str(int(timei[i]))[-4:] for i in range(len(timei))]
    time = []
    for i in range(147000):
        x = dt.datetime.strptime(timei[0], '%m%d') + i * dt.timedelta(days=1)    
        if x.month == 2 and x.day == 29: continue
        time.append(x)
    time = np.array(time[delete_days:])
else: time = timei[delete_days:]

pos_lats_NH = np.where((lats > 0))[0]
lats_NH = lats[pos_lats_NH]
pos_lats_hw = np.where((lats_NH >= lat_minHW) & (lats_NH <= lat_maxHW))[0]   
pos_lons_hw = np.where((lons >= lon_minHW) & (lons <= lon_maxHW))[0]   


# File with the positions of heat waves 
name_file_posHW = f'resume_positions_HWdays_{name}_Teng.csv'
pos_HWdays = pd.read_csv(f'{path_outputs}{name_file_posHW}', index_col=0)

# ==================================================== Temperature ==============================================================================
ncfile = Dataset(f'{path_data}/{nc_name_t}') 
t = np.array(ncfile[var_t])[:,pos_lats_NH,:]
t = t[delete_days:,:,:]
t = t[:len(time),:,:]

if seasons == True: 
    df_t = pd.DataFrame(index=time, data=np.reshape(t, [t.shape[0], t.shape[1] * t.shape[2]]))
    t_anom = anomalies_seasons(df_t)
elif seasons == False: 
    df_t = pd.DataFrame(data=np.reshape(t, [t.shape[0], t.shape[1] * t.shape[2]]))
    t_anom = anomalies_noseasons(df_t)

t_anom = t_anom.values.reshape(t_anom.shape[0], lats_NH.shape[0], lons.shape[0])

composites_matriz_t, composites_matrix_complete_t = calculate_composites2(pos_HWdays, t_anom)

# ==================================================== v' ============================================================================== 
ncfile = Dataset(f'{path_data}/{nc_name_v}') 
v_300 = np.array(ncfile[var_v])[:,pos_lats_NH,:]
v_300 = v_300[delete_days:]
v_300 = v_300[:len(time)]

if seasons == True: 
    df_v_300 = pd.DataFrame(index=time, data=np.reshape(v_300, [v_300.shape[0], v_300.shape[1] * v_300.shape[2]]))
    v300_anom = anomalies_seasons(df_v_300)
elif seasons == False: 
    df_v_300 = pd.DataFrame(data=np.reshape(v_300, [v_300.shape[0], v_300.shape[1] * v_300.shape[2]]))
    v300_anom = anomalies_noseasons(df_v_300)
v300_anom = v300_anom.values.reshape(v300_anom.shape[0], lats_NH.shape[0], lons.shape[0])

composites_matriz_v, composites_matrix_complete_v = calculate_composites2(pos_HWdays, v300_anom)


# ==================================================== Streamfunction ==============================================================================  
ncfile = Dataset(f'{path_data}/{nc_name_sf}')
sf = np.array(ncfile[var_sf])[:,pos_lats_NH,:]
sf = sf[delete_days:]
sf = sf[:len(time)]


# Calculating the 'subseasonal anomalies' for the streamfunction. Following (Teng, 2013), calculate 
# both: anomalies with respect to the climatological day and anomalies with respect to the actual season
if seasons == True: 
    df_sf = pd.DataFrame(index=time, data=np.reshape(sf, [sf.shape[0], sf.shape[1] * sf.shape[2]]))
    sf_anom_1 = anomalies_seasons(df_sf)
    sf_anom_sub = subseasonal_anomalies(sf_anom_1)
elif seasons == False: 
    df_sf = pd.DataFrame(data=np.reshape(sf, [sf.shape[0], sf.shape[1] * sf.shape[2]]))
    sf_anom_sub = anomalies_noseasons(df_sf)
sf_anom = sf_anom_sub.values.reshape(sf_anom_sub.shape[0], lats_NH.shape[0], lons.shape[0])

composites_matriz_sf, composites_matrix_complete_sf = calculate_composites2(pos_HWdays, sf_anom)
composites_matriz_sf = composites_matriz_sf/10000000
composites_matrix_complete_sf =composites_matrix_complete_sf/10000000


# =================================================== RWP Envelope ==========================================================
# Filtering between 
kmin = 8
kmax = 16
RWP_envelope = np.zeros_like(v300_anom)*np.nan
for t in range(v300_anom.shape[0]):
    for num, lat in enumerate(lats_NH):
        serie = v300_anom[t,num,:]
        
        analytic_signal_filt = hilbert_filtered(serie, kmin, kmax)
        analytic_signal_hilbert_filt = np.abs(analytic_signal_filt)
        RWP_envelope[t,num,:] = analytic_signal_hilbert_filt

composites_matriz_envelope, composites_matrix_complete_envelope = calculate_composites2(pos_HWdays, RWP_envelope)


# =================================================== COMPOSITES ==========================================================
# ================================================================================================================

# ================================================ HOVMOLLER DIAGRAM ==========================================================
# ================================================================================================================
pos_midlat = np.where(abs(lats_NH - midlat) == np.min(abs(lats_NH - midlat)))[0][0]
pos_middle_HW = np.where(abs(lons - (lon_minHW+lon_maxHW)/2) == np.min(abs(lons - (lon_minHW+lon_maxHW)/2)))[0][0]
lons_hovmoller = np.roll(lons, - (round(abs(pos_middle_HW-len(lons)/2))))
time_lags = np.arange(-20, 21, 1)


# Hovmoller diagram for one particular latitud in the midlatitudes
composites_matrix_complete_t_midlat = composites_matrix_complete_t[:,pos_midlat,:]
hovmoller_t = np.roll(composites_matrix_complete_t_midlat,  -(round(abs(pos_middle_HW-len(lons)/2))), axis = 1)

composites_matrix_complete_sf_midlat = composites_matrix_complete_sf[:,pos_midlat,:]
hovmoller_sf = np.roll(composites_matrix_complete_sf_midlat,  -(round(abs(pos_middle_HW-len(lons)/2))), axis = 1)

composites_matrix_complete_v_midlat = composites_matrix_complete_v[:,pos_midlat,:]
hovmoller_v = np.roll(composites_matrix_complete_v_midlat,  -(round(abs(pos_middle_HW-len(lons)/2))), axis = 1)

composites_matrix_complete_env_midlat = composites_matrix_complete_envelope[:,pos_midlat,:]
hovmoller_envelope = np.roll(composites_matrix_complete_env_midlat,  -(round(abs(pos_middle_HW-len(lons)/2))), axis = 1)




# Hovmoller diagram for a band of latitudes
pos_lats_hw = np.where((lats_NH >= 30) & (lats_NH <= 60))[0]
composites_matrix_complete_t_midlat_mean = np.mean(composites_matrix_complete_t[:,pos_lats_hw,:], 1)
hovmoller_t_mean = np.roll(composites_matrix_complete_t_midlat_mean,  -(round(abs(pos_middle_HW-len(lons)/2))), axis = 1)

composites_matrix_complete_sf_midlat_mean = np.mean(composites_matrix_complete_sf[:,pos_lats_hw,:], 1)
hovmoller_sf_mean = np.roll(composites_matrix_complete_sf_midlat_mean,  -(round(abs(pos_middle_HW-len(lons)/2))), axis = 1)

composites_matrix_complete_v_midlat_mean = np.mean(composites_matrix_complete_v[:,pos_lats_hw,:], 1)
hovmoller_v_mean = np.roll(composites_matrix_complete_v_midlat_mean,  -(round(abs(pos_middle_HW-len(lons)/2))), axis = 1)

composites_matrix_complete_env_midlat_mean = np.mean(composites_matrix_complete_envelope[:,pos_lats_hw,:], 1)
hovmoller_envelope_mean = np.roll(composites_matrix_complete_env_midlat_mean,  -(round(abs(pos_middle_HW-len(lons)/2))), axis = 1)

   

if name == 'NCEP':
    li, ls, intervalos, limite, color = -3, 3, 15, 1, 'RdYlBu_r'
    bounds = np.round(np.linspace(li, ls, intervalos), 3)
    colormap = center_white_anom(color, intervalos, bounds, limite)

    composites_coastlines(lats_NH, lons, composites_matriz_t, -3, 3 + np.abs(bounds[0] - bounds[1]), composites_matriz_sf, [-0.2,0.1,0.2], composites_matriz_envelope, [10, 12],lat_minHW, lat_maxHW, lon_minHW, lon_maxHW, colormap, f'{path_figures}composites_evolution_coastlines_v200_t_sf_E.png', var = 'T anomaly [°C]')

    # Particular event 
    pos_hw_event = 26828 
    year = 2021
    composites_matriz_t, composites_matrix_complete_t = calculate_composites_event(pos_hw_event, t_anom)
    composites_matriz_sf, composites_matrix_complete_sf = calculate_composites_event(pos_hw_event, sf_anom)
    composites_matriz_sf = composites_matriz_sf/10000000
    composites_matrix_complete_sf =composites_matrix_complete_sf/10000000
    composites_matriz_envelope, composites_matrix_complete_envelope = calculate_composites_event(pos_hw_event, RWP_envelope)

    li, ls, intervalos, limite, color = -10, 10, 15, 4, 'coolwarm'
    bounds = np.round(np.linspace(li, ls, intervalos), 3)
    colormap = center_white_anom(color, intervalos, bounds, limite)

    composites_coastlines(lats_NH, lons, composites_matriz_t, -10, 10 + np.abs(bounds[0] - bounds[1]), composites_matriz_sf, [-1.5, -0.5, 0.5, 1.5], composites_matriz_envelope, [17, 20,25,28],lat_minHW, lat_maxHW, lon_minHW, lon_maxHW, colormap, f'{path_figures}composites_evolution_coastlines_v200_t_sf_E_{year}.png', var = 'T anomaly [°C]', topography=topography)

    li, ls, intervalos, limite, color = -4, 4, 15, 1.2, 'RdYlBu_r'
    bounds_t = np.round(np.linspace(li, ls, intervalos), 3)
    colormap_t = center_white_anom(color, intervalos, bounds_t, limite)
    norm_t = mcolors.DivergingNorm(vmin=-4, vcenter=0, vmax = 4)

    li, ls, intervalos, limite, color = -1.1, 1.1, 15, 0.1, 'RdGy_r'
    bounds_sf = np.round(np.linspace(li, ls, intervalos), 3)
    colormap_sf = center_white_anom(color, intervalos, bounds_sf, limite)
    norm_sf = mcolors.DivergingNorm(vmin=-1.1, vcenter=0, vmax = 1.1)

    hovmoller(time_lags, np.linspace(-180,180,len(lons)), hovmoller_t, colormap_t, norm_t, -4, 4, hovmoller_sf, colormap_sf, norm_sf, -1.1, 1.1, hovmoller_envelope, [18, 21], path_figures + f'Hovmoller.png', var_1 = 'T anomalies [K]', var_2 = r'Streamfunction anomalies [10 x $10^{7}$ $m^{2}$/s]')

    li, ls, intervalos, limite, color = -2, 2, 15, 0.4, 'RdYlBu_r' 
    bounds_t = np.round(np.linspace(li, ls, intervalos), 3)
    colormap_t = center_white_anom(color, intervalos, bounds_t, limite)
    norm_t = mcolors.DivergingNorm(vmin=-2, vcenter=0, vmax = 2)

    li, ls, intervalos, limite, color = -0.34, 0.34, 15, 0.03, 'RdGy_r'
    bounds_sf = np.round(np.linspace(li, ls, intervalos), 3)
    colormap_sf = center_white_anom(color, intervalos, bounds_sf, limite)
    norm_sf = mcolors.DivergingNorm(vmin=-0.35, vcenter=0, vmax = 0.35)

    hovmoller(time_lags, np.linspace(-180,180,len(lons)), hovmoller_t_mean, colormap_t, norm_t, -2, 2, hovmoller_sf_mean, colormap_sf, norm_sf, -0.35, 0.35, hovmoller_envelope_mean, [8.5,8.8], path_figures + f'Hovmoller_mean_midlat.png', var_1 = 'T anomalies [K]', var_2 = r'Streamfunction anomalies [1x$10^{7}$ $m^{2}$/s]')


elif name in ['CAM', 'LENS']:
    li, ls, intervalos, limite, color = -6, 6, 15, 1, 'RdYlBu_r'
    bounds = np.round(np.linspace(li, ls, intervalos), 3)
    colormap = center_white_anom(color, intervalos, bounds, limite)
    composites_coastlines(lats_NH, lons, composites_matriz_t, -6, 6 + np.abs(bounds[0] - bounds[1]), composites_matriz_sf, [-1,-0.2,0.2,1], composites_matriz_envelope, [10, 15],lat_minHW, lat_maxHW, lon_minHW, lon_maxHW, colormap, f'{path_figures}composites_evolution_coastlines_v200_t_sf_E.png', var = 'T anomaly [°C]', topography=topography)

    li, ls, intervalos, limite, color = -4, 4, 15, 1.2, 'RdYlBu_r'
    bounds_t = np.round(np.linspace(li, ls, intervalos), 3)
    colormap_t = center_white_anom(color, intervalos, bounds_t, limite)
    norm_t = mcolors.DivergingNorm(vmin=-4, vcenter=0, vmax = 4)

    li, ls, intervalos, limite, color = -1.1, 1.1, 15, 0.1, 'RdGy_r'
    bounds_sf = np.round(np.linspace(li, ls, intervalos), 3)
    colormap_sf = center_white_anom(color, intervalos, bounds_sf, limite)
    norm_sf = mcolors.DivergingNorm(vmin=-1.1, vcenter=0, vmax = 1.1)

    hovmoller(time_lags, np.linspace(-180,180,len(lons)), hovmoller_t, colormap_t, norm_t, -4, 4, hovmoller_sf, colormap_sf, norm_sf, -1.1, 1.1, hovmoller_envelope, [18, 21], path_figures + f'Hovmoller.png', var_1 = 'T anomalies [K]', var_2 = r'Streamfunction anomalies [10 x $10^{7}$ $m^{2}$/s]')

    li, ls, intervalos, limite, color = -2, 2, 15, 0.4, 'RdYlBu_r' 
    bounds_t = np.round(np.linspace(li, ls, intervalos), 3)
    colormap_t = center_white_anom(color, intervalos, bounds_t, limite)
    norm_t = mcolors.DivergingNorm(vmin=-2, vcenter=0, vmax = 2)

    li, ls, intervalos, limite, color = -1.1, 1.1, 15, 0.1, 'RdGy_r'
    bounds_sf = np.round(np.linspace(li, ls, intervalos), 3)
    colormap_sf = center_white_anom(color, intervalos, bounds_sf, limite)
    norm_sf = mcolors.DivergingNorm(vmin=-1.1, vcenter=0, vmax = 1.1)

    hovmoller(time_lags, np.linspace(-180,180,len(lons)), hovmoller_t_mean, colormap_t, norm_t, -2, 2, hovmoller_sf_mean, colormap_sf, norm_sf, -0.35, 0.35, hovmoller_envelope_mean, [8.5,8.8], path_figures + f'Hovmoller_mean_midlat.png', var_1 = 'T anomalies [K]', var_2 = r'Streamfunction anomalies [1x$10^{7}$ $m^{2}$/s]')


elif name == 'exp1_Held_Suarez':
    li, ls, intervalos, limite, color = -5, 5, 15, 1, 'coolwarm'
    bounds = np.round(np.linspace(li, ls, intervalos), 3)
    colormap = center_white_anom(color, intervalos, bounds, limite)
    composites_coastlines(lats_NH, lons, composites_matriz_t, -5, 5 + np.abs(bounds[0] - bounds[1]), composites_matriz_sf, [-0.17,0.07,0.17], composites_matriz_envelope, [28, 30, 36, 39], lat_minHW, lat_maxHW, lon_minHW, lon_maxHW, colormap, f'{path_figures}composites_evolution_t_sf_E.png', var = 'T anomaly [°C]', topography=topography)

    li, ls, intervalos, limite, color = -5, 5, 15, 1.2, 'RdYlBu_r'
    bounds_t = np.round(np.linspace(li, ls, intervalos), 3)
    colormap_t = center_white_anom(color, intervalos, bounds_t, limite)
    norm_t = mcolors.DivergingNorm(vmin=-5, vcenter=0, vmax = 5)

    li, ls, intervalos, limite, color = -1.2, 1.2, 15, 0.3, 'RdGy_r'
    bounds_sf = np.round(np.linspace(li, ls, intervalos), 3)
    colormap_sf = center_white_anom(color, intervalos, bounds_sf, limite)
    norm_sf = mcolors.DivergingNorm(vmin=-1.2, vcenter=0, vmax = 1.2)

    hovmoller(time_lags, np.linspace(-180,180,len(lons)), hovmoller_t, colormap_t, norm_t, -5, 5, hovmoller_sf, colormap_sf, norm_sf, -1.2, 1.2 + np.abs(bounds_sf[0] - bounds_sf[1]), hovmoller_envelope, [27,30,33], path_figures + f'Hovmoller.png', var_1 = 'T anomalies [K]', var_2 = r'Streamfunction anomalies [10 x $10^{7}$ $m^{2}$/s]')


    li, ls, intervalos, limite, color = -4, 4, 15, 0.7, 'RdYlBu_r'
    bounds_t = np.round(np.linspace(li, ls, intervalos), 3)
    colormap_t = center_white_anom(color, intervalos, bounds_t, limite)
    norm_t = mcolors.DivergingNorm(vmin=-5, vcenter=0, vmax = 5)

    li, ls, intervalos, limite, color = -0.6, 0.6, 15, 0.09, 'RdGy_r'
    bounds_sf = np.round(np.linspace(li, ls, intervalos), 3)
    colormap_sf = center_white_anom(color, intervalos, bounds_sf, limite)
    norm_sf = mcolors.DivergingNorm(vmin=-0.2, vcenter=0, vmax = 0.18)

    hovmoller(time_lags, np.linspace(-180,180,len(lons)), hovmoller_t_mean, colormap_t, norm_t, -4, 4, hovmoller_sf_mean, colormap_sf, norm_sf, -0.6, 0.6 + np.abs(bounds_sf[0] - bounds_sf[1]), hovmoller_envelope_mean, [25, 27, 29], path_figures + f'Hovmoller_mean_midlat.png', var_1 = 'T anomalies [K]', var_2 = r'Streamfunction anomalies [10 x $10^{7}$ $m^{2}$/s]')



elif name == 'exp2_NCEPsymm_noSeason_noTop':
    li, ls, intervalos, limite, color = -2, 2, 15, 0.4, 'coolwarm'
    bounds = np.round(np.linspace(li, ls, intervalos), 3)
    colormap = center_white_anom(color, intervalos, bounds, limite)
    composites_coastlines(lats_NH, lons, composites_matriz_t, -2, 2 + np.abs(bounds[0] - bounds[1]), composites_matriz_sf, [-0.02,-0.01,0.01,0.02], composites_matriz_envelope, [8, 10], lat_minHW, lat_maxHW, lon_minHW, lon_maxHW, colormap, f'{path_figures}composites_evolution_t_sf_E.png', var = 'T anomaly [°C]', topography=topography)

    li, ls, intervalos, limite, color = -2, 2, 15, 0.4, 'RdYlBu_r'
    bounds_t = np.round(np.linspace(li, ls, intervalos), 3)
    colormap_t = center_white_anom(color, intervalos, bounds_t, limite)
    norm_t = mcolors.DivergingNorm(vmin=-2, vcenter=0, vmax = 2)

    li, ls, intervalos, limite, color = -0.3, 0.3, 15, 0.09, 'RdGy_r'
    bounds_sf = np.round(np.linspace(li, ls, intervalos), 3)
    colormap_sf = center_white_anom(color, intervalos, bounds_sf, limite)
    norm_sf = mcolors.DivergingNorm(vmin=-0.3, vcenter=0, vmax = 0.3)

    hovmoller(time_lags, np.linspace(-180,180,len(lons)), hovmoller_t, colormap_t, norm_t, -2, 2, hovmoller_sf, colormap_sf, norm_sf, -0.3, 0.3 + np.abs(bounds_sf[0] - bounds_sf[1]), hovmoller_envelope, [8.8,9.5], path_figures + f'Hovmoller.png', var_1 = 'T anomalies [K]', var_2 = r'Streamfunction anomalies [10 x $10^{7}$ $m^{2}$/s]')

    li, ls, intervalos, limite, color = -1.4, 1.4, 15, 0.25, 'RdYlBu_r'
    bounds_t = np.round(np.linspace(li, ls, intervalos), 3)
    colormap_t = center_white_anom(color, intervalos, bounds_t, limite)
    norm_t = mcolors.DivergingNorm(vmin=-1.4, vcenter=0, vmax = 1.4)

    li, ls, intervalos, limite, color = -0.18, 0.18, 15, 0.06, 'RdGy_r'
    bounds_sf = np.round(np.linspace(li, ls, intervalos), 3)
    colormap_sf = center_white_anom(color, intervalos, bounds_sf, limite)
    norm_sf = mcolors.DivergingNorm(vmin=-0.18, vcenter=0, vmax = 0.18)

    hovmoller(time_lags, np.linspace(-180,180,len(lons)), hovmoller_t, colormap_t, norm_t, -1.4, 1.4, hovmoller_sf, colormap_sf, norm_sf, -0.18, 0.18 + np.abs(bounds_sf[0] - bounds_sf[1]), hovmoller_envelope, [5.3, 5.5], path_figures + f'Hovmoller_mean_midlat.png', var_1 = 'T anomalies [K]', var_2 = r'Streamfunction anomalies [10 x $10^{7}$ $m^{2}$/s]')


elif name == 'exp3_NCEPasymm_noSeason_noTop':
    li, ls, intervalos, limite, color = -2, 2, 15, 0.4, 'coolwarm'
    bounds = np.round(np.linspace(li, ls, intervalos), 3)
    colormap = center_white_anom(color, intervalos, bounds, limite)
    composites_coastlines(lats_NH, lons, composites_matriz_t, -2, 2 + np.abs(bounds[0] - bounds[1]), composites_matriz_sf, [-0.07,-0.025,0.025,0.07], composites_matriz_envelope, [7,8.5], lat_minHW, lat_maxHW, lon_minHW, lon_maxHW, colormap, f'{path_figures}composites_evolution_t_sf_E.png', var = 'T anomaly [°C]', topography=topography)

    li, ls, intervalos, limite, color = -2, 2, 15, 0.4, 'RdYlBu_r'
    bounds_t = np.round(np.linspace(li, ls, intervalos), 3)
    colormap_t = center_white_anom(color, intervalos, bounds_t, limite)
    norm_t = mcolors.DivergingNorm(vmin=-2, vcenter=0, vmax = 2)

    li, ls, intervalos, limite, color = -0.3, 0.3, 15, 0.09, 'RdGy_r'
    bounds_sf = np.round(np.linspace(li, ls, intervalos), 3)
    colormap_sf = center_white_anom(color, intervalos, bounds_sf, limite)
    norm_sf = mcolors.DivergingNorm(vmin=-0.3, vcenter=0, vmax = 0.3)

    hovmoller(time_lags, np.linspace(-180,180,len(lons)), hovmoller_t, colormap_t, norm_t, -2, 2, hovmoller_sf, colormap_sf, norm_sf, -0.3, 0.3 + np.abs(bounds_sf[0] - bounds_sf[1]), hovmoller_envelope, [8.8,9.5], path_figures + f'Hovmoller.png', var_1 = 'T anomalies [K]', var_2 = r'Streamfunction anomalies [10 x $10^{7}$ $m^{2}$/s]')

    li, ls, intervalos, limite, color = -1.4, 1.4, 15, 0.25, 'RdYlBu_r'
    bounds_t = np.round(np.linspace(li, ls, intervalos), 3)
    colormap_t = center_white_anom(color, intervalos, bounds_t, limite)
    norm_t = mcolors.DivergingNorm(vmin=-1.4, vcenter=0, vmax = 1.4)

    li, ls, intervalos, limite, color = -0.18, 0.18, 15, 0.06, 'RdGy_r'
    bounds_sf = np.round(np.linspace(li, ls, intervalos), 3)
    colormap_sf = center_white_anom(color, intervalos, bounds_sf, limite)
    norm_sf = mcolors.DivergingNorm(vmin=-0.18, vcenter=0, vmax = 0.18)

    hovmoller(time_lags, np.linspace(-180,180,len(lons)), hovmoller_t, colormap_t, norm_t, -1.4, 1.4, hovmoller_sf, colormap_sf, norm_sf, -0.18, 0.18 + np.abs(bounds_sf[0] - bounds_sf[1]), hovmoller_envelope, [5.3, 5.5], path_figures + f'Hovmoller_mean_midlat.png', var_1 = 'T anomalies [K]', var_2 = r'Streamfunction anomalies [10 x $10^{7}$ $m^{2}$/s]')
