import datetime as dt
from dateutil.relativedelta import relativedelta
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import pandas as pd
import numpy as np
from cartopy import crs
import cartopy
import matplotlib.ticker as mticker
import math
import os
import matplotlib.colors as mcolors
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.cm as cm

name = 'exp7_Held_Suarez'
topography = False
delete_days = 1000
#lat_minHW = 30; lat_maxHW = 50; lon_minHW = 245; lon_maxHW = 290; midlat=45   #EXP7, EXP8 - V1
#lat_minHW = 30; lat_maxHW = 50; lon_minHW = 255; lon_maxHW = 280; midlat=45  #EXP9

#lat_minHW = 30; lat_maxHW = 55; lon_minHW = 240; lon_maxHW = 270; midlat=50   #EXP10
lat_minHW = 30; lat_maxHW = 55; lon_minHW = 245; lon_maxHW = 290; midlat=50  #EXP7

def center_white_anom(cmap, num, bounds, limite):
    import matplotlib as mpl
    barcmap = mpl.cm.get_cmap(cmap, num)
    barcmap.set_bad(color='white', alpha=0.5)
    bar_vals = barcmap(np.arange(num))  # extract those values as an array
    pos = np.arange(num)
    centro = pos[(bounds >= -limite) & (bounds <= limite)]
    for i in centro:
        bar_vals[i] = [1, 1, 1, 1]  # change the middle value
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("new" + cmap, bar_vals)
    return newcmap

def maps1(x, y, minn, maxx, matriz,  cmap, path, norm, topography='',  units=''):
    fig = plt.figure(figsize=[8, 6])
    
    if topography == True:
        ax = fig.add_subplot(1, 1, 1, projection=crs.PlateCarree())
        #ax.outline_patch.set_edgecolor('None')
        ax.add_feature(cartopy.feature.COASTLINE, lw=0.5, zorder=11)
    else: ax= fig.add_subplot(1, 1, 1)

    im = ax.contourf(x, y[:], matriz[:,:], levels=20, cmap = cmap, vmin=minn, vmax=maxx, norm=norm)
    r = ax.contour(x, y, matriz, levels=20, colors='k', linewidths=0.5)
    #ax.invert_yaxis()

    #cbaxes = fig.add_axes([0.2, -0.1, 0.6, 0.030])
    #cb = plt.colorbar(im, orientation="horizontal", pad=0.1, cax=cbaxes, format='%.1f')
    cb = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), orientation="horizontal", pad=0.1, format='%.1f', shrink=0.8)
    cb.set_label(units, fontsize=9, color='dimgrey')
    cb.outline.set_edgecolor(None)
    cb.ax.tick_params(labelcolor='dimgrey', color='dimgrey', labelsize=9)
    plt.savefig(path, dpi=200)
    #plt.show()
    plt.close()


def calculate_composites(pos_HW, matriz):
    dic_pos_composites = {}

    time_lags = np.arange(-20, 21, 1)
    for pos in pos_HW.index:
        if pos == 0: 
            dic_pos_composites[0] = []
            dic_pos_composites[0].append(pos_HW.iloc[0][0])
            continue
        elif pos == pos_HW.index[-1]:
            break
        elif pos_HW.iloc[pos][0]-1 != pos_HW.iloc[pos-1][0]:
            for num in time_lags:  
                if num in dic_pos_composites.keys(): dic_pos_composites[num].append(pos_HW.iloc[pos][0]+num)  
                else: 
                    dic_pos_composites[num] = []
                    dic_pos_composites[num].append(pos_HW.iloc[pos][0]+num)


    day_0 = np.mean(matriz[dic_pos_composites[0]], axis=0)
    day_5 = np.mean(matriz[dic_pos_composites[-5]], axis=0)
    day_10 = np.mean(matriz[dic_pos_composites[-10]], axis=0)
    day_15 = np.mean(matriz[dic_pos_composites[-15]], axis=0)
    day_20 = np.mean(matriz[dic_pos_composites[-20]], axis=0)
    day_5_post = np.mean(matriz[dic_pos_composites[5]], axis=0)

    composites_matriz = np.zeros((6, day_0.shape[0], day_0.shape[1]))
    for i, matriz_i in enumerate([day_20, day_15, day_10, day_5, day_0, day_5_post]):
        composites_matriz[i,:,:] = matriz_i
    
    return composites_matriz


def calculate_composites2(pos_HW, matriz):
    dic_composites = {}

    time_lags = np.arange(-20, 21, 1)
    for pos in pos_HW.index:
        if pos == 0: 
            dic_composites[0] = []
            dic_composites[0].append(pos_HW.iloc[0][0])
            continue
        elif pos == pos_HW.index[-1]:
            break
        elif pos_HW.iloc[pos][0]-1 != pos_HW.iloc[pos-1][0]:
            for num in time_lags:  
                if num in dic_composites.keys(): dic_composites[num].append(pos_HW.iloc[pos][0]+num)  
                else: 
                    dic_composites[num] = []
                    dic_composites[num].append(pos_HW.iloc[pos][0]+num)

    composites_matrix_complete =  np.zeros((len(time_lags), matriz.shape[1], matriz.shape[2]))
    for ii, lag in enumerate(dic_composites.keys()):
        composites_matrix_complete[ii] = np.mean(matriz[dic_composites[lag]], axis = 0)

    day_0 = np.mean(matriz[dic_composites[0]], axis=0)
    day_5 = np.mean(matriz[dic_composites[-5]], axis=0)
    day_10 = np.mean(matriz[dic_composites[-10]], axis=0)
    day_15 = np.mean(matriz[dic_composites[-15]], axis=0)
    day_20 = np.mean(matriz[dic_composites[-20]], axis=0)
    day_5_post = np.mean(matriz[dic_composites[5]], axis=0)

    composites_matrix = np.zeros((6, day_0.shape[0], day_0.shape[1]))
    for i, matriz_i in enumerate([day_20, day_15, day_10, day_5, day_0, day_5_post]):
        composites_matrix[i,:,:] = matriz_i
    
    return composites_matrix, composites_matrix_complete


def composites(lats, lons, Matriz_t, min_t, max_t, Matriz_sf, levels_sf, Matriz_env, levels_env, lat_minHW, lat_maxHW, lon_minHW, lon_maxHW, cmap, path, var = ''):
    fig = plt.figure(figsize=[15, 21])
    for i, tit in enumerate(['Day -20', 'Day -15', 'Day -10', 'Day -5', 'Day 0', 'Day 5']):
        ax = fig.add_subplot(6, 1, i + 1, projection=crs.PlateCarree(central_longitude=180))
        ax.outline_patch.set_edgecolor('None')
        #ax.add_feature(cartopy.feature.COASTLINE, lw=0.5, zorder=11)

        im = ax.contourf(lons, lats, Matriz_t[i, :, :], cmap=cmap, extend='both', levels=np.arange(min_t, max_t, 0.1),transform=crs.PlateCarree())
        im2 = ax.contour(lons, lats, Matriz_env[i, :, :], extend='both', levels=levels_env, colors='k', linewidths=1.5,transform=crs.PlateCarree())
        #im2 = ax.contour(lons, lats, Matriz_env[i, :, :], extend='both', levels=[25, 30, 35], colors='k', linewidths=1.5,transform=crs.PlateCarree())
        ax.clabel(im2, inline=True, fontsize=10, fmt='%1.1f')
        im3 = ax.contour(lons, lats, Matriz_sf[i, :, :], extend='both', levels=levels_sf, colors='limegreen', negative_linestyles = 'dashed', linewidths=1.7,transform=crs.PlateCarree())
      
        ax.plot([lon_minHW, lon_maxHW], [lat_minHW, lat_minHW], transform=crs.PlateCarree(), color='b', lw=1.5)
        ax.plot([lon_minHW, lon_maxHW], [lat_maxHW, lat_maxHW], transform=crs.PlateCarree(), color='b', lw=1.5)
        ax.plot([lon_minHW, lon_minHW], [lat_minHW, lat_maxHW], transform=crs.PlateCarree(), color='b', lw=1.5)
        ax.plot([lon_maxHW, lon_maxHW], [lat_minHW, lat_maxHW], transform=crs.PlateCarree(), color='b', lw=1.5)
        plt.ylim(25,70)
        #plt.ylim(15,60)

        ax.set_title(f'{tit}', fontsize=15, color='k')

    cbaxes = fig.add_axes([0.3, 0.06, 0.4, 0.015])
    cb = plt.colorbar(im, orientation="horizontal", pad=0.2, cax=cbaxes, format='%.2f')
    cb.set_label(var, fontsize=15, color='dimgrey')
    cb.outline.set_edgecolor('dimgrey')
    cb.ax.tick_params(labelcolor='dimgrey', color='dimgrey', labelsize=14)
    plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.9,
                    top=0.92,
                    wspace=0.4,
                    hspace=0.4)
    #plt.show()
    plt.savefig(path, dpi=200)
    plt.close()


def composites_coastlines(lats, lons, Matriz_t, min_t, max_t, Matriz_sf, levels_sf, Matriz_env, levels_env, lat_minHW, lat_maxHW, lon_minHW, lon_maxHW, cmap, path, var = ''):
    fig = plt.figure(figsize=[15, 20])
    for i, tit in enumerate(['Day -20', 'Day -15', 'Day -10', 'Day -5', 'Day 0', 'Day 5']):
        ax = fig.add_subplot(6, 1, i + 1, projection=crs.PlateCarree(central_longitude=180))
        ax.outline_patch.set_edgecolor('None')
        ax.add_feature(cartopy.feature.COASTLINE, lw=0.5, zorder=11)

        im = ax.contourf(lons, lats, Matriz_t[i, :, :], cmap=cmap, extend='both', levels=np.arange(min_t, max_t, 0.1),transform=crs.PlateCarree())
        im2 = ax.contour(lons, lats, Matriz_env[i, :, :], extend='both', levels=levels_env, colors='k', linewidths=1.5,transform=crs.PlateCarree())
        ax.clabel(im2, inline=True, fontsize=10, fmt='%1.1f')
        #im2 = ax.contour(lons, lats, Matriz_env[i, :, :], extend='both', levels=[25, 30, 35], colors='k', linewidths=1.5,transform=crs.PlateCarree())
        im3 = ax.contour(lons, lats, Matriz_sf[i, :, :], extend='both', levels=levels_sf, colors='mediumorchid', negative_linestyles = 'dashed', linewidths=1.4,transform=crs.PlateCarree())
        ax.clabel(im2, inline=True, fontsize=10, fmt='%1.1f')

        ax.plot([lon_minHW, lon_maxHW], [lat_minHW, lat_minHW], transform=crs.PlateCarree(), color='b', lw=1.5)
        ax.plot([lon_minHW, lon_maxHW], [lat_maxHW, lat_maxHW], transform=crs.PlateCarree(), color='b', lw=1.5)
        ax.plot([lon_minHW, lon_minHW], [lat_minHW, lat_maxHW], transform=crs.PlateCarree(), color='b', lw=1.5)
        ax.plot([lon_maxHW, lon_maxHW], [lat_minHW, lat_maxHW], transform=crs.PlateCarree(), color='b', lw=1.5)
        plt.ylim(25,70)
        #plt.ylim(15,60)

        gl = ax.gridlines(crs=crs.PlateCarree(), draw_labels=True,linewidth=0.8, color='gray', alpha=0.5, linestyle='--')
        gl.top_labels = False
        # gl.ylocator = mticker.FixedLocator([45, 60])
        # gl.left_labels = True
        # gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 14, 'color': 'dimgrey'}
        gl.ylabel_style = {'size': 14, 'color': 'dimgrey'}
        ax.set_title(f'{tit}', fontsize=15, color='k')

    cbaxes = fig.add_axes([0.3, 0.06, 0.4, 0.015])
    cb = plt.colorbar(im, orientation="horizontal", pad=0.2, cax=cbaxes, format='%.2f')
    cb.set_label(var, fontsize=15, color='dimgrey')
    cb.outline.set_edgecolor('dimgrey')
    cb.ax.tick_params(labelcolor='dimgrey', color='dimgrey', labelsize=14)
    plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.9,
                    top=0.92,
                    wspace=0.3,
                    hspace=0.3)
    #plt.show()
    plt.savefig(path, dpi=200)
    plt.close()



from scipy import linalg, fft as sp_fft
def hilbert_filtered(x, k_min, k_max, N=None, axis=-1, path=''):
    """
    Compute the analytic signal, using the Hilbert transform.

    The transformation is done along the last axis by default.

    Parameters
    ----------
    x : array_like
        Signal data.  Must be real.
    N : int, optional
        Number of Fourier components.  Default: ``x.shape[axis]``
    axis : int, optional
        Axis along which to do the transformation.  Default: -1.

    Returns
    -------
    xa : ndarray
        Analytic signal of `x`, of each 1-D array along `axis`

    Notes
    -----
    The analytic signal ``x_a(t)`` of signal ``x(t)`` is:

    .. math:: x_a = F^{-1}(F(x) 2U) = x + i y

    where `F` is the Fourier transform, `U` the unit step function,
    and `y` the Hilbert transform of `x`. [1]_

    In other words, the negative half of the frequency spectrum is zeroed
    out, turning the real-valued signal into a complex signal.  The Hilbert
    transformed signal can be obtained from ``np.imag(hilbert(x))``, and the
    original signal from ``np.real(hilbert(x))``.

    Examples
    --------
    In this example we use the Hilbert transform to determine the amplitude
    envelope and instantaneous frequency of an amplitude-modulated signal.

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy.signal import hilbert, chirp

    >>> duration = 1.0
    >>> fs = 400.0
    >>> samples = int(fs*duration)
    >>> t = np.arange(samples) / fs

    We create a chirp of which the frequency increases from 20 Hz to 100 Hz and
    apply an amplitude modulation.

    >>> signal = chirp(t, 20.0, t[-1], 100.0)
    >>> signal *= (1.0 + 0.5 * np.sin(2.0*np.pi*3.0*t) )

    The amplitude envelope is given by magnitude of the analytic signal. The
    instantaneous frequency can be obtained by differentiating the
    instantaneous phase in respect to time. The instantaneous phase corresponds
    to the phase angle of the analytic signal.

    >>> analytic_signal = hilbert(signal)
    >>> amplitude_envelope = np.abs(analytic_signal)
    >>> instantaneous_phase = np.unwrap(np.angle(analytic_signal))
    >>> instantaneous_frequency = (np.diff(instantaneous_phase) /
    ...                            (2.0*np.pi) * fs)

    >>> fig, (ax0, ax1) = plt.subplots(nrows=2)
    >>> ax0.plot(t, signal, label='signal')
    >>> ax0.plot(t, amplitude_envelope, label='envelope')
    >>> ax0.set_xlabel("time in seconds")
    >>> ax0.legend()
    >>> ax1.plot(t[1:], instantaneous_frequency)
    >>> ax1.set_xlabel("time in seconds")
    >>> ax1.set_ylim(0.0, 120.0)
    >>> fig.tight_layout()

    References
    ----------
    .. [1] Wikipedia, "Analytic signal".
           https://en.wikipedia.org/wiki/Analytic_signal
    .. [2] Leon Cohen, "Time-Frequency Analysis", 1995. Chapter 2.
    .. [3] Alan V. Oppenheim, Ronald W. Schafer. Discrete-Time Signal
           Processing, Third Edition, 2009. Chapter 12.
           ISBN 13: 978-1292-02572-8

    """
    x = np.asarray(x)
    if np.iscomplexobj(x):
        raise ValueError("x must be real.")
    if N is None:
        N = x.shape[axis]
    if N <= 0:
        raise ValueError("N must be positive.")

    Xf = sp_fft.fft(x, N, axis=axis)
    pow = np.abs(Xf/len(x))**2
    freq = np.fft.fftfreq(len(x), 1)

    # fig = plt.figure(figsize = [10,4])
    # plt.plot(freq*360, pow*100/np.sum(pow),'o-', color='dimgray', markersize=1.8)
    # plt.xscale('log')
    # plt.legend(loc="upper right", frameon=True)
    # plt.ylabel('[%] Potencia')
    # plt.xlabel('Frecuencia')
    # plt.title('Espectro  de Fourier')
    # #plt.show()
    # #plt.savefig(path, dpi=200) 
    # plt.close()

    h = np.zeros(N)
    if N % 2 == 0:
        h[0] = h[N // 2] = 1
        h[1:N // 2] = 2
    else:
        h[0] = 1
        h[1:(N + 1) // 2] = 2

    filter = np.where((np.abs(freq*360) < k_min) ^ (np.abs(freq*360) > k_max))[0]
    Xf_filtered = np.copy(Xf)
    Xf_filtered[filter] = 0

    if x.ndim > 1:
        ind = [np.newaxis] * x.ndim
        ind[axis] = slice(None)
        h = h[tuple(ind)]
    x = sp_fft.ifft(Xf_filtered * h, axis=axis)
    return x

# def anomalias_GCMdrycore(df_VAR):
#     ANOMA = df_VAR * np.nan
#     STDAR = df_VAR * np.nan
#     media = df_VAR.mean()
#     std = df_VAR.std()
#     for i in range(df_VAR.shape[0]):
#         dia = df_VAR.loc[i]
#         anoma = dia - media
#         ANOMA.iloc[i] = anoma
#         STDAR.iloc[i] = anoma / std
#     return ANOMA, STDAR

def anomalias_GCMdrycore(df_VAR):
    ANOMA = df_VAR * np.nan
    media = df_VAR.mean()
    for i in range(df_VAR.shape[0]):
        dia = df_VAR.loc[i]
        anoma = dia - media
        ANOMA.iloc[i] = anoma
    return ANOMA


import matplotlib.colors
RdYlBu_list = ['rgb(165,0,38)','rgb(215,48,39)','rgb(244,109,67)','rgb(253,174,97)','rgb(254,224,144)','rgb(255,255,191)','rgb(224,243,248)','rgb(171,217,233)','rgb(116,173,209)','rgb(69,117,180)','rgb(49,54,149)']
my_cmap = matplotlib.colors.ListedColormap(RdYlBu_list, name='RdYlBu')
def hovmoller(time_lags, lons, Matriz_t, cmap_t, norm_t, min_t, max_t, Matriz_sf, cmap_sf, norm_sf, min_sf, max_sf, Matriz_env, levels_env, path, var_1 = '', var_2 = ''):
    fig = plt.figure(figsize=[8,4])
    ax = fig.add_subplot(1, 2, 1)
    im = ax.contourf(lons, time_lags, Matriz_t[:,:], cmap=cmap_t, extend='both', levels=np.linspace(min_t, max_t, 15))
    plt.ylim(-15,15)
    plt.ylabel(r'Time lag', fontsize=8)
    plt.xlabel(r'Relative longitude', fontsize=8)
    ax.set_yticks([-15,-10,-5,0,5,10,15])
    ax.set_xticks([-180,-90,0,90,180])
    ax.tick_params(labelsize=7)
    ax.axhline(0, ls = '--', color='dimgray', lw = 0.3)
    ax.axvline(0, ls = '--', color='dimgray', lw = 0.3)
    cbaxes = fig.add_axes([0.11, 0.12, 0.35, 0.03])
    #cb = plt.colorbar(cm.ScalarMappable(norm=norm_t, cmap=cmap_t), orientation="horizontal", pad=0.2, cax=cbaxes, format='%.2f')
    cb = plt.colorbar(im, orientation="horizontal", pad=0.2, cax=cbaxes, format='%.1f')
    cb.set_label(var_1, fontsize=7, color='dimgrey')
    cb.outline.set_edgecolor('dimgrey')
    cb.ax.tick_params(labelcolor='dimgrey', color='dimgrey', labelsize=8)

    ax2 = fig.add_subplot(1, 2, 2)
    im2 = ax2.contour(lons, time_lags, Matriz_env[:, :], extend='both', levels=levels_env, colors='k', linewidths=0.9)
    ax2.clabel(im2, inline=True, fontsize=6, fmt='%1.1f')
    im3 = ax2.contourf(lons, time_lags, Matriz_sf[:, :],  cmap=cmap_sf, extend='both', levels=np.linspace(min_sf, max_sf, 15))
    plt.ylim(-15,15)
    plt.xlabel(r'Relative longitude', fontsize=8)
    ax2.set_yticks([-15,-10,-5,0,5,10,15])
    ax2.set_xticks([-180,-90,0,90,180])
    ax2.tick_params(labelsize=7)
    ax2.axhline(0, ls = '--', color='dimgray', lw = 0.3)
    ax2.axvline(0, ls = '--', color='dimgray', lw = 0.3)
    cbaxes = fig.add_axes([0.55, 0.12, 0.35, 0.03])
    #cb = plt.colorbar(cm.ScalarMappable(norm=norm_sf, cmap=cmap_sf), orientation="horizontal", pad=0.2, cax=cbaxes, format='%.2f')
    cb = plt.colorbar(im3, orientation="horizontal", pad=0.2, cax=cbaxes, format='%1.2f')
    cb.set_label(var_2, fontsize=7, color='dimgrey')
    cb.outline.set_edgecolor('dimgrey')
    cb.ax.tick_params(labelcolor='dimgrey', color='dimgrey', labelsize=7)

    plt.subplots_adjust(left=0.08,
                    bottom=0.3,
                    right=0.92,
                    top=0.92,
                    wspace=0.15,
                    hspace=0.2)
    plt.savefig(path, dpi=200)
    plt.close()


path_home = '/home/castanev/'
path_data = f'/scratch/brown/castanev/DryCore_Wu/output/{name}/post_processed/output/'
path_figures = f'{path_home}Data Analysis/DryCore_Wu/{name}/Figures/'
path_outputs = f'{path_home}Data Analysis/DryCore_Wu/{name}/'

# Lats and lons 
nc_name = f't.atmos_daily.nc'
ncfile = Dataset(f'{path_data}{nc_name}')
lats = np.array(ncfile['lat'])  # 64. 
lons = np.array(ncfile['lon'])  # 128. 

pos_lats_NH = np.where((lats >= 0))[0]
lats_NH = lats[pos_lats_NH]

pos_lats_hw = np.where((lats_NH >= lat_minHW) & (lats_NH <= lat_maxHW))[0]   
pos_lons_hw = np.where((lons >= lon_minHW) & (lons <= lon_maxHW))[0]   


# File with the positions of heat waves 
name_file_posHW = f'resume_positions_HWdays_GCM_Teng.csv'
pos_HWdays = pd.read_csv(f'{path_outputs}{name_file_posHW}', index_col=0)


# =================================================== Temperature ==========================================================
t_k   = np.array(ncfile['temp'])  #(t, lat, lon) °K

t = pd.DataFrame(data=np.reshape(t_k, [t_k.shape[0], t_k.shape[1] * t_k.shape[2]]))
t = t - 273.15 # °C
t = t.values.reshape(t_k.shape[0], t_k.shape[1], t_k.shape[2]) #(t, lat, lon) °C

# Deleting the first days and obtaining NH
t = t[delete_days:,:,:]
t_NH = t[:,pos_lats_NH,:]

VAR 		= pd.DataFrame(data = np.reshape(t_NH,[t_NH.shape[0],t_NH.shape[1]*t_NH.shape[2]]))
t_NH_anom = anomalias_GCMdrycore(VAR)
t_NH_anom = t_NH_anom.values.reshape(t_NH_anom.shape[0], t_NH.shape[1], t_NH.shape[2])

composites_matriz_t, composites_matrix_complete_t = calculate_composites2(pos_HWdays, t_NH_anom)

t_mean_2d = np.mean(t_k[delete_days:,:,:], axis = 0)
# if topography == False:
#     path = path_figures + f'Mean_t_surface.png'
#     norm = mcolors.DivergingNorm(vmin=184, vcenter=273.15, vmax = 300)
#     maps1(lons, lats, 184, 300, t_mean_2d, 'coolwarm', path, norm, topography, 'T [K]')
# elif topography == True:
#     path = path_figures + f'Mean_t_surface.png'
#     norm = mcolors.DivergingNorm(vmin=225, vcenter=273.15, vmax = 314)
#     maps1(lons, lats, 225, 314, t_mean_2d, 'coolwarm', path, norm, topography, 'T [K]')

# ================================================== Streamfunction ==========================================================
# name_file = f'anom_streamfunction_d_north.nc'
# file = Dataset(f'{path_data}{name_file}')
# sf_anom = np.array(file['anom_sf'])

name_file = f'streamfunction_daily.nc'
file = Dataset(f'{path_data}{name_file}')
sf   = np.array(file['SF'])  #(t, lat, lon) 

# Deleting the first days and obtaining NH
sf = sf[delete_days:,:,:]
sf_NH = sf[:,pos_lats_NH,:]

VAR 		= pd.DataFrame(data = np.reshape(sf_NH,[sf_NH.shape[0],sf_NH.shape[1]*sf_NH.shape[2]]))
sf_anom = anomalias_GCMdrycore(VAR)
sf_anom = sf_anom.values.reshape(sf_anom.shape[0], sf_NH.shape[1], sf_NH.shape[2])

composites_matriz_sf, composites_matrix_complete_sf = calculate_composites2(pos_HWdays, sf_anom)
composites_matriz_sf = composites_matriz_sf/10000000
composites_matrix_complete_sf =composites_matrix_complete_sf/10000000


# =================================================== RWP Envelope ==========================================================
name_file = f'anom_v300_d_north.nc'
file = Dataset(f'{path_data}{name_file}')
v300_anom = np.array(file['anom_v300'])


kmin = 1
kmax = 30
RWP_envelope = np.zeros_like(v300_anom)*np.nan
for t in range(v300_anom.shape[0]):
    for num, lat in enumerate(lats_NH):
        serie = v300_anom[t,num,:]
        
        analytic_signal_filt = hilbert_filtered(serie, kmin, kmax)
        analytic_signal_hilbert_filt = np.abs(analytic_signal_filt)

        RWP_envelope[t,num,:] = analytic_signal_hilbert_filt

composites_matriz_envelope, composites_matrix_complete_envelope = calculate_composites2(pos_HWdays, RWP_envelope)


# ==================================================== Geopotential height ==============================================================================
nc_name = f'z500.atmos_daily.nc'  
ncfile = Dataset(f'{path_data}/{nc_name}') 
z500_NH = np.array(ncfile['z500'])[:,pos_lats_NH,:]
z500_NH = z500_NH[delete_days:]

composites_matriz_z500, composites_matrix_complete_z500 = calculate_composites2(pos_HWdays, z500_NH)


nc_name = f'z1000.atmos_daily.nc'  
ncfile = Dataset(f'{path_data}/{nc_name}') 
z1000_NH = np.array(ncfile['z1000'])[:,pos_lats_NH,:]
z1000_NH = z1000_NH[delete_days:]

VAR = pd.DataFrame(data=np.reshape(z1000_NH, [z1000_NH.shape[0], z1000_NH.shape[1] * z1000_NH.shape[2]]))
z1000_NH_anom= anomalias_GCMdrycore(VAR)
z1000_NH_anom = z1000_NH_anom.values.reshape(z1000_NH_anom.shape[0], lats_NH.shape[0], lons.shape[0])

composites_matriz_z1000, composites_matrix_complete_z1000 = calculate_composites2(pos_HWdays, z1000_NH_anom)



# ================================================ COMPOSITES AND HOVMOLLER DIAGRAM (FIGURES) ==========================================================
# ================================================================================================================
pos_midlat = np.where(abs(lats_NH - midlat) == np.min(abs(lats_NH - midlat)))[0][0]
pos_middle_HW = np.where(abs(lons - (lon_minHW+lon_maxHW)/2) == np.min(abs(lons - (lon_minHW+lon_maxHW)/2)))[0][0]
lons_hovmoller = np.roll(lons, - (round(abs(pos_middle_HW-len(lons)/2))))
time_lags = np.arange(-20, 21, 1)
pos_50 = np.where(abs(lats_NH - 50) == np.min(abs(lats_NH - 50)))[0][0]



composites_matrix_complete_t_midlat = composites_matrix_complete_t[:,pos_midlat,:]
hovmoller_t = np.roll(composites_matrix_complete_t_midlat,  -(round(abs(pos_middle_HW-len(lons)/2))), axis = 1)

composites_matrix_complete_sf_midlat = composites_matrix_complete_sf[:,pos_midlat,:]
hovmoller_sf = np.roll(composites_matrix_complete_sf_midlat,  -(round(abs(pos_middle_HW-len(lons)/2))), axis = 1)

composites_matrix_complete_env_midlat = composites_matrix_complete_envelope[:,pos_midlat,:]
hovmoller_envelope = np.roll(composites_matrix_complete_env_midlat,  -(round(abs(pos_middle_HW-len(lons)/2))), axis = 1)




composites_matrix_complete_t_midlat_mean = np.mean(composites_matrix_complete_t[:,pos_lats_hw,:], 1)
hovmoller_t_mean = np.roll(composites_matrix_complete_t_midlat_mean,  -(round(abs(pos_middle_HW-len(lons)/2))), axis = 1)

composites_matrix_complete_sf_midlat_mean = np.mean(composites_matrix_complete_sf[:,pos_lats_hw,:], 1)
hovmoller_sf_mean = np.roll(composites_matrix_complete_sf_midlat_mean,  -(round(abs(pos_middle_HW-len(lons)/2))), axis = 1)

composites_matrix_complete_env_midlat_mean = np.mean(composites_matrix_complete_envelope[:,pos_lats_hw,:], 1)
hovmoller_envelope_mean = np.roll(composites_matrix_complete_env_midlat_mean,  -(round(abs(pos_middle_HW-len(lons)/2))), axis = 1)



if name == 'exp7_Held_Suarez':
    #EXP7
    li, ls, intervalos, limite, color = -5, 5, 15, 1, 'coolwarm'
    bounds = np.round(np.linspace(li, ls, intervalos), 3)
    colormap = center_white_anom(color, intervalos, bounds, limite)
    composites(lats_NH, lons, composites_matriz_t, -5, 5 + np.abs(bounds[0] - bounds[1]), composites_matriz_z1000, [-2,-1.1,1.1, 2], composites_matriz_z500, [5200, 5300, 5400], lat_minHW, lat_maxHW, lon_minHW, lon_maxHW, colormap, f'{path_figures}composites_evolution_t_z500_z1000.png', var = 'T anomaly [°C]') 
    composites(lats_NH, lons, composites_matriz_t, -5, 5 + np.abs(bounds[0] - bounds[1]), composites_matriz_sf, [-0.17,0.07,0.17], composites_matriz_envelope, [28, 30, 36, 39], lat_minHW, lat_maxHW, lon_minHW, lon_maxHW, colormap, f'{path_figures}composites_evolution_t_sf_E.png', var = 'T anomaly [°C]')

    # EXP7
    li, ls, intervalos, limite, color = -5, 5, 15, 1.2, 'RdYlBu_r'
    bounds_t = np.round(np.linspace(li, ls, intervalos), 3)
    colormap_t = center_white_anom(color, intervalos, bounds_t, limite)
    norm_t = mcolors.DivergingNorm(vmin=-5, vcenter=0, vmax = 5)

    li, ls, intervalos, limite, color = -1.2, 1.2, 15, 0.3, 'RdGy_r'
    bounds_sf = np.round(np.linspace(li, ls, intervalos), 3)
    colormap_sf = center_white_anom(color, intervalos, bounds_sf, limite)
    norm_sf = mcolors.DivergingNorm(vmin=-1.2, vcenter=0, vmax = 1.2)

    hovmoller(time_lags, np.linspace(-180,180,len(lons)), hovmoller_t, colormap_t, norm_t, -5, 5, hovmoller_sf, colormap_sf, norm_sf, -1.2, 1.2 + np.abs(bounds_sf[0] - bounds_sf[1]), hovmoller_envelope, [27,30,33], path_figures + f'Hovmoller.png', var_1 = 'T anomalies [K]', var_2 = r'Streamfunction anomalies [10 x $10^{7}$ $m^{2}$/s]')


    #EXP7
    li, ls, intervalos, limite, color = -4, 4, 15, 0.7, 'RdYlBu_r'
    bounds_t = np.round(np.linspace(li, ls, intervalos), 3)
    colormap_t = center_white_anom(color, intervalos, bounds_t, limite)
    norm_t = mcolors.DivergingNorm(vmin=-5, vcenter=0, vmax = 5)

    li, ls, intervalos, limite, color = -0.6, 0.6, 15, 0.09, 'RdGy_r'
    bounds_sf = np.round(np.linspace(li, ls, intervalos), 3)
    colormap_sf = center_white_anom(color, intervalos, bounds_sf, limite)
    norm_sf = mcolors.DivergingNorm(vmin=-0.2, vcenter=0, vmax = 0.18)

    hovmoller(time_lags, np.linspace(-180,180,len(lons)), hovmoller_t_mean, colormap_t, norm_t, -4, 4, hovmoller_sf_mean, colormap_sf, norm_sf, -0.6, 0.6 + np.abs(bounds_sf[0] - bounds_sf[1]), hovmoller_envelope_mean, [25, 27, 29], path_figures + f'Hovmoller_mean_midlat.png', var_1 = 'T anomalies [K]', var_2 = r'Streamfunction anomalies [10 x $10^{7}$ $m^{2}$/s]')



elif name == 'exp8_NCEPsymm_noSeason_noTop':
    #EXP8
    li, ls, intervalos, limite, color = -2, 2, 15, 0.4, 'coolwarm'
    bounds = np.round(np.linspace(li, ls, intervalos), 3)
    colormap = center_white_anom(color, intervalos, bounds, limite)
    composites(lats_NH, lons, composites_matriz_t, -2, 2 + np.abs(bounds[0] - bounds[1]), composites_matriz_z1000, [-1.5,-0.8,0.8, 1.5], composites_matriz_z500, [5300, 5400, 5450], lat_minHW, lat_maxHW, lon_minHW, lon_maxHW, colormap, f'{path_figures}composites_evolution_t_z500_z1000.png', var = 'T anomaly [°C]') 
    composites(lats_NH, lons, composites_matriz_t, -2, 2 + np.abs(bounds[0] - bounds[1]), composites_matriz_sf, [-0.02,-0.01,0.01,0.02], composites_matriz_envelope, [8, 10], lat_minHW, lat_maxHW, lon_minHW, lon_maxHW, colormap, f'{path_figures}composites_evolution_t_sf_E.png', var = 'T anomaly [°C]')

    #EXP8
    li, ls, intervalos, limite, color = -2, 2, 15, 0.4, 'RdYlBu_r'
    bounds_t = np.round(np.linspace(li, ls, intervalos), 3)
    colormap_t = center_white_anom(color, intervalos, bounds_t, limite)
    norm_t = mcolors.DivergingNorm(vmin=-2, vcenter=0, vmax = 2)

    li, ls, intervalos, limite, color = -0.3, 0.3, 15, 0.09, 'RdGy_r'
    bounds_sf = np.round(np.linspace(li, ls, intervalos), 3)
    colormap_sf = center_white_anom(color, intervalos, bounds_sf, limite)
    norm_sf = mcolors.DivergingNorm(vmin=-0.3, vcenter=0, vmax = 0.3)

    hovmoller(time_lags, np.linspace(-180,180,len(lons)), hovmoller_t, colormap_t, norm_t, -2, 2, hovmoller_sf, colormap_sf, norm_sf, -0.3, 0.3 + np.abs(bounds_sf[0] - bounds_sf[1]), hovmoller_envelope, [8.8,9.5], path_figures + f'Hovmoller.png', var_1 = 'T anomalies [K]', var_2 = r'Streamfunction anomalies [10 x $10^{7}$ $m^{2}$/s]')

    #EXP8
    li, ls, intervalos, limite, color = -1.4, 1.4, 15, 0.25, 'RdYlBu_r'
    bounds_t = np.round(np.linspace(li, ls, intervalos), 3)
    colormap_t = center_white_anom(color, intervalos, bounds_t, limite)
    norm_t = mcolors.DivergingNorm(vmin=-1.4, vcenter=0, vmax = 1.4)

    li, ls, intervalos, limite, color = -0.18, 0.18, 15, 0.06, 'RdGy_r'
    bounds_sf = np.round(np.linspace(li, ls, intervalos), 3)
    colormap_sf = center_white_anom(color, intervalos, bounds_sf, limite)
    norm_sf = mcolors.DivergingNorm(vmin=-0.18, vcenter=0, vmax = 0.18)

    hovmoller(time_lags, np.linspace(-180,180,len(lons)), hovmoller_t, colormap_t, norm_t, -1.4, 1.4, hovmoller_sf, colormap_sf, norm_sf, -0.18, 0.18 + np.abs(bounds_sf[0] - bounds_sf[1]), hovmoller_envelope, [5.3, 5.5], path_figures + f'Hovmoller_mean_midlat.png', var_1 = 'T anomalies [K]', var_2 = r'Streamfunction anomalies [10 x $10^{7}$ $m^{2}$/s]')


elif name == 'exp10_NCEPasymm_noSeason_noTop':
    #EXP10
    li, ls, intervalos, limite, color = -2, 2, 15, 0.4, 'coolwarm'
    bounds = np.round(np.linspace(li, ls, intervalos), 3)
    colormap = center_white_anom(color, intervalos, bounds, limite)
    composites(lats_NH, lons, composites_matriz_t, -2, 2 + np.abs(bounds[0] - bounds[1]), composites_matriz_z1000, [-1,-0.2,0.2, 1], composites_matriz_z500, [5300, 5400, 5500], lat_minHW, lat_maxHW, lon_minHW, lon_maxHW, colormap, f'{path_figures}composites_evolution_t_z500_z1000.png', var = 'T anomaly [°C]') 
    composites(lats_NH, lons, composites_matriz_t, -2, 2 + np.abs(bounds[0] - bounds[1]), composites_matriz_sf, [-0.07,-0.025,0.025,0.07], composites_matriz_envelope, [7,8.5], lat_minHW, lat_maxHW, lon_minHW, lon_maxHW, colormap, f'{path_figures}composites_evolution_t_sf_E.png', var = 'T anomaly [°C]')

    #EXP10
    li, ls, intervalos, limite, color = -2, 2, 15, 0.4, 'RdYlBu_r'
    bounds_t = np.round(np.linspace(li, ls, intervalos), 3)
    colormap_t = center_white_anom(color, intervalos, bounds_t, limite)
    norm_t = mcolors.DivergingNorm(vmin=-2, vcenter=0, vmax = 2)

    li, ls, intervalos, limite, color = -0.3, 0.3, 15, 0.09, 'RdGy_r'
    bounds_sf = np.round(np.linspace(li, ls, intervalos), 3)
    colormap_sf = center_white_anom(color, intervalos, bounds_sf, limite)
    norm_sf = mcolors.DivergingNorm(vmin=-0.3, vcenter=0, vmax = 0.3)

    hovmoller(time_lags, np.linspace(-180,180,len(lons)), hovmoller_t, colormap_t, norm_t, -2, 2, hovmoller_sf, colormap_sf, norm_sf, -0.3, 0.3 + np.abs(bounds_sf[0] - bounds_sf[1]), hovmoller_envelope, [8.8,9.5], path_figures + f'Hovmoller.png', var_1 = 'T anomalies [K]', var_2 = r'Streamfunction anomalies [10 x $10^{7}$ $m^{2}$/s]')

    #EXP10
    li, ls, intervalos, limite, color = -1.4, 1.4, 15, 0.25, 'RdYlBu_r'
    bounds_t = np.round(np.linspace(li, ls, intervalos), 3)
    colormap_t = center_white_anom(color, intervalos, bounds_t, limite)
    norm_t = mcolors.DivergingNorm(vmin=-1.4, vcenter=0, vmax = 1.4)

    li, ls, intervalos, limite, color = -0.18, 0.18, 15, 0.06, 'RdGy_r'
    bounds_sf = np.round(np.linspace(li, ls, intervalos), 3)
    colormap_sf = center_white_anom(color, intervalos, bounds_sf, limite)
    norm_sf = mcolors.DivergingNorm(vmin=-0.18, vcenter=0, vmax = 0.18)

    hovmoller(time_lags, np.linspace(-180,180,len(lons)), hovmoller_t, colormap_t, norm_t, -1.4, 1.4, hovmoller_sf, colormap_sf, norm_sf, -0.18, 0.18 + np.abs(bounds_sf[0] - bounds_sf[1]), hovmoller_envelope, [5.3, 5.5], path_figures + f'Hovmoller_mean_midlat.png', var_1 = 'T anomalies [K]', var_2 = r'Streamfunction anomalies [10 x $10^{7}$ $m^{2}$/s]')
