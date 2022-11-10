# This code is for the detection of heatwaves. It's based on the methodology used by Teng, 2013
# Here, a heat wave event is defined as at least 5 consecutive days following:
# i). More than 5% of the domain (US) has daily averaged SAT exceeding the threshold value
# ii). Centre of these warm points does not move faster than 5◦ latitude or longitude per day

# Threshold: 97.5 percentile for historical t within a 15-day window centred on the day (for each day and each gridpoint)
# To avoid contamination, we use only events that have no heat wave days in the preceding 20 day

# NOTE: 
# For the detection of heat waves in the Dry Core GCM, the threshold is defined as the 97.5 percentile of all the data (the model has no seasons)

from Functions import *
from netCDF4 import Dataset
import pandas as pd
import numpy as np
import pandas as pd
import datetime as dt
from dateutil.relativedelta import relativedelta
import matplotlib.colors
import os

import matplotlib.pyplot as plt
import cartopy
from cartopy import crs
import matplotlib.ticker as mticker
import matplotlib.ticker as ticker


# INPUTS TO CHANGE ===========================================================================================================================
# name: It can be NCEP or the name of the experiment (output of the model)
# nc_name: name of the file containing the surface air temperature. 't.atmos_daily.nc' for the dry core GCM or 't1000_NCEP.nc' for NCEP
# var: the variable name in the nc file. 't1000' for NCEP, 'temp' for the dry core GCM 
# min_duration: minimun consecutive days for a heat wave event
# vel: vel is a vector to limit the speed of the center of the warm points [max. degrees, number of days]
# threshold_value: percentile 
# lat_minHW, lat_maxHW, lon_minHW, lon_minHW: to define the area where the events will be detected
# delete_days is the number of days to cut (only for output models). 

name = 'NCEP'  
nc_name = 't1000_NCEP.nc'
var = 't1000'
vel = [5, 1]
min_duration = 5
threshold_value = 97.5
lat_minHW = 30; lat_maxHW = 50; lon_minHW = 245; lon_maxHW = 290; midlat=45
seasons = True 
topography = True
delete_days = ''

# ============================================================================================================================================
path_data = f'{os.path.abspath(os.getcwd())}/{name}/Data/'
path_figures = f'{os.path.abspath(os.getcwd())}/{name}/Figures/'
path_outputs = f'{os.path.abspath(os.getcwd())}/{name}/'

methodology = 'Teng'
vel_str = f'vel{str(vel[0])}{str(vel[1])}'

ncfile = Dataset(f'{path_data}{nc_name}')
t_k = np.array(ncfile[var])  #(t, lat, lon) °K
t = pd.DataFrame(data=np.reshape(t_k, [t_k.shape[0], t_k.shape[1] * t_k.shape[2]]))
t = t - 273.15 # °C
t = t.values.reshape(t_k.shape[0], t_k.shape[1], t_k.shape[2])
lats = np.array(ncfile['lat'])
lons = np.array(ncfile['lon'])
time = np.array(ncfile['time'])


# Spatial cut: United States
pos_lats = np.where((lats >= lat_minHW) & (lats <= lat_maxHW))
pos_lons = np.where((lons >= lon_minHW) & (lons <= lon_maxHW))
lats_US = lats[pos_lats]
lons_US = lons[pos_lons]
min_grid_5 = round(lats_US.shape[0] * lons_US.shape[0] * 0.05)
t_US = t[:, pos_lats[0], :]
t_US = t_US[:, :, pos_lons[0]]

min_grid_5 = round(lats_US.shape[0] * lons_US.shape[0] * 0.05)
print(f"Total de días analizados = {t_US.shape[0]}")
print(f'5% of the total grids corresponds to = {min_grid_5}')  


if name == 'NCEP':
    dates_d = np.array([dt.datetime(1800,1,1) + dt.timedelta(hours = int(time[i])) for i in range(len(time))])
    Month = np.array([ii.month for ii in dates_d])
    df_t = pd.DataFrame(index=dates_d, data=np.reshape(t_US, [t_US.shape[0], t_US.shape[1] * t_US.shape[2]]))
    days_summer = np.array([dt.datetime(2021, 6, 1) + relativedelta(days=int(xx)) for xx in range(92)])  # random year
    
    threshold = np.zeros([len(days_summer), len(lats_US), len(lons_US)])
    for i, d in enumerate(days_summer):
        t_pos1 = df_t.index.get_indexer_for((df_t.loc[(df_t.index.month == d.month) & (df_t.index.day == d.day)].index))
        # Positions of the 15-day window centred on the day of the year of the potential heat wave day:
        t_pos = [np.concatenate((t_pos1, t_pos1 + iii)) for iii in range(-7, 8)]
        t_pos = np.unique(t_pos)
        t_pos = t_pos[t_pos < t_US.shape[0]] # For the last days of the dataset
        data = t_US[t_pos, :, :]
        threshold_i = np.percentile(data, threshold_value, axis=0)
        threshold[i, :, :] = threshold_i

    df_threshold = pd.DataFrame(index=days_summer,
                                data=np.reshape(threshold, [days_summer.shape[0], t_US.shape[1] * t_US.shape[2]]))

    pos_summer = np.where([ii in [6, 7, 8] for ii in Month])[0]
    dates_summer = dates_d[pos_summer]
    t_summer = t_US[pos_summer, :, :]

    pos_heat_wavesi = []
    for i, date in enumerate(dates_summer):

        threshold_pos = df_threshold.index.get_indexer_for(
            (df_threshold.loc[(df_threshold.index.month == date.month) & (df_threshold.index.day == date.day)].index))
        
        # Condition i). More than 5% of the domain (US) has daily averaged SAT exceeding the threshold value
        cond_1 = t_US[pos_summer[i], :, :] > threshold[threshold_pos, :, :]
        grid_cont = np.count_nonzero(cond_1)

        # Condition ii): Centre of these warm points does not move faster than 5◦ latitude or longitude per day
        # Center defined as the point with the max temperature in the domain
        if pos_summer[i] + vel[1] == t_US.shape[0]: break
        max_i = np.max(t_US[pos_summer[i], :, :])
        pos_max = np.where(t_US[pos_summer[i], :, :] == max_i)
        max_i2 = np.max(t_US[pos_summer[i] + vel[1], :, :])
        pos_max2 = np.where(t_US[pos_summer[i] + vel[1], :, :] == max_i2)

        dif_lats = np.abs(lats_US[pos_max2[0]] - lats_US[pos_max[0]])
        dif_lons = np.abs(lons_US[pos_max2[1]] - lons_US[pos_max[1]])

        if grid_cont > min_grid_5 and dif_lats.all() <= vel[0] and dif_lons.all() <= vel[0]:
            pos_heat_wavesi.append(pos_summer[i])  # dates_d[pos_summer[i]] is the date that meets both conditions



elif (name != 'NCEP') and (seasons == False): 
    t = t[delete_days:,:,:]
    df_t = pd.DataFrame(data=np.reshape(t_US, [t_US.shape[0], t_US.shape[1] * t_US.shape[2]]))
    threshold = np.percentile(t_US, threshold_value, axis=0)
    
    pos_heat_wavesi = []
    for pos in range(t_US.shape[0]):
        # Condition i). More than 5% of the domain (US) has daily averaged SAT exceeding the threshold value
        cond_1 = t_US[pos, :, :] > threshold
        grid_cont = np.count_nonzero(cond_1)

        # Condition ii): Centre of these warm points does not move faster than 5◦ latitude or longitude per day
        # Center defined as the point with the max temperature in the domain     
        if pos + vel[1] == t_US.shape[0]: break
        max_i = np.max(t_US[pos, :, :])
        pos_max = np.where(t_US[pos, :, :] == max_i)
        max_i2 = np.max(t_US[pos + vel[1], :, :])
        pos_max2 = np.where(t_US[pos + vel[1], :, :] == max_i2)

        dif_lats = np.abs(lats_US[pos_max2[0]] - lats_US[pos_max[0]])
        dif_lons = np.abs(lons_US[pos_max2[1]] - lons_US[pos_max[1]])

        if grid_cont > min_grid_5 and dif_lats.all() <= vel[0] and dif_lons.all() <= vel[0]:
            pos_heat_wavesi.append(pos)  # dates_d[pos_summer[i]] is the date that meets both conditions



duration_hw, pos_day1_hw = duration_heat_waves(pos_heat_wavesi, min_duration)
print(f'Number of heat waves events: {len(duration_hw)}')

# Probability distribution function for duration
bins = np.arange(np.unique(duration_hw)[0], np.unique(duration_hw)[-1] + 1, 1)
Hist, bins1 = np.histogram(duration_hw, len(bins))
PDF_temp = Hist / len(duration_hw)


fig = plt.figure(figsize=[4, 4])
df = pd.DataFrame({'x': bins, 'PDF': PDF_temp})
df.plot.bar(x='x', y='PDF', rot=0, color='dimgray', width=.4)
plt.ylabel('PDF', fontsize=13, color = 'dimgray')
plt.xlabel('Duration (d)', fontsize=13, color = 'dimgray')
plt.xticks(fontsize=12, color = 'dimgray')
plt.yticks(fontsize=12, color = 'dimgray')
plt.legend(fontsize=12, labelcolor='linecolor')
#plt.show()
plt.savefig(path_figures + f'SAT_PDF_{methodology}_{vel_str}.png', dpi=200)
plt.close()


pos_heat_waves = []
for day1, dur in zip(pos_day1_hw, duration_hw): 
    pos_heat_waves.append([day1 + i for i in range(dur)])
pos_heat_waves = [item for sublist in pos_heat_waves for item in sublist]
print(f'Number heat waves days: {len(pos_heat_waves)}')


# Saving the resume of statistics in a .csv
resume = pd.read_csv(f'{path_outputs}/../resume_heatWaves_statistics.csv', index_col = 0)
resume.loc[name, f'HW days'] = len(pos_heat_waves)
resume.loc[name, f'HW events'] = len(duration_hw)
resume.to_csv(f'{path_outputs}/../resume_heatWaves_statistics.csv')

# Saving the position of HW days in a .csv
pos_heat_waves_serie = pd.Series(pos_heat_waves)
pos_heat_waves_serie.to_csv(f'{path_outputs}/resume_positions_HWdays_{model}_{methodology}.csv')

t_US_hw = t_US[pos_heat_waves, :, :]


# Calculating anomalies
if seasons == True: anom_t_US = anomalies_seasons(df_t)
elif seasons == False: anom_t_US = anomalies_noseasons(df_t)
else: print('Variable "seasons" must be defined') 

anom_t_US = anom_t_US.values.reshape(anom_t_US.shape[0], lats_US.shape[0], lons_US.shape[0])
anom_t_US_hw = anom_t_US[pos_heat_waves, :, :]


# Ploting intensity
RdYlBu_v2_list = ['rgb(224,243,248)','rgb(171,217,233)','rgb(116,173,209)','rgb(69,117,180)','rgb(49,54,149)']
my_cmap = matplotlib.colors.ListedColormap(RdYlBu_v2_list, name='RdYlBu')

li, ls, intervalos, limite, color = 0.7, 3.3, 15, 1.6, 'RdYlBu_r'
bounds = np.round(np.linspace(li, ls, intervalos), 3)
colormap = center_white_anom(color, intervalos, bounds, limite)
maps2(lons_US, np.ndarray.round(lats_US,2), 0.7, 3.3 + np.abs(bounds[0] - bounds[1]), np.nanmean(anom_t_US_hw, axis=0), r'SAT anomalies [°C]',
      colormap, path_figures + f'SAT_mean_anomt_hw_{methodology}_{vel_str}.png', topography)
