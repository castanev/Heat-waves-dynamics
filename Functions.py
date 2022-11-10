

# ================================= PLOTS ===============================================
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



def maps1(x, y, minn, maxx, matriz,  cmap, path, norm, units=''):
    fig = plt.figure(figsize=[8, 6])
    ax = fig.add_subplot(1, 1, 1)
    im = ax.contourf(x, y[:], matriz[:,:], levels=20, cmap = cmap, vmin=minn, vmax=maxx, norm=norm)
    r = ax.contour(x, y, matriz, levels=20, colors='k', linewidths=0.5)
    #ax.invert_yaxis()
    cb = plt.colorbar(im, orientation="horizontal", pad=0.1, format='%.1f', shrink=0.8)
    cb.set_label(units, fontsize=9, color='dimgrey')
    cb.outline.set_edgecolor(None)
    cb.ax.tick_params(labelcolor='dimgrey', color='dimgrey', labelsize=9)
    plt.savefig(path, dpi=200)
    plt.close()



def maps2(Lons, Lats, minn, maxx, matriz, var, cmap, path, topography):
    fig = plt.figure(figsize=[6.5, 4.5])
    ax = fig.add_subplot(1, 1, 1, projection=crs.PlateCarree(central_longitude=180))
    
    if topography == True:
        ax.add_feature(cartopy.feature.BORDERS, lw=0.5)
        ax.outline_patch.set_edgecolor('None')
        ax.add_feature(cartopy.feature.COASTLINE, lw=0.5, zorder=11)

    im = ax.contourf(Lons, Lats, matriz, cmap=cmap, extend='both', \
                     levels=np.arange(minn, maxx, 0.2), transform=crs.PlateCarree())

    gl = ax.gridlines(crs=crs.PlateCarree(), draw_labels=True,
                      linewidth=0.7, color='gray', alpha=0.2, linestyle='--')

    gl.top_labels = False
    gl.right_labels = False
    gl.ylocator = mticker.FixedLocator(Lats[np.arange(0, len(Lats), 3)])
    gl.xlabel_style = {'size': 10, 'color': 'dimgrey'}
    gl.ylabel_style = {'size': 10, 'color': 'dimgrey'}

    cbaxes = fig.add_axes([0.2, 0.12, 0.6, 0.030])
    cb = plt.colorbar(im, orientation="horizontal", pad=0.1, cax=cbaxes, format='%.1f')
    cb.set_label(var, fontsize=10, color='dimgrey')
    cb.outline.set_edgecolor(None)
    cb.ax.tick_params(labelcolor='dimgrey', color='dimgrey', labelsize=9)
    plt.savefig(path, dpi=200)
    plt.close()


def composites_coastlines(lats, lons, Matriz_t, min_t, max_t, Matriz_sf, levels_sf, Matriz_env, levels_env, lat_minHW, lat_maxHW, lon_minHW, lon_maxHW, cmap, path, var = ''):
    fig = plt.figure(figsize=[15, 20])
    for i, tit in enumerate(['Day -20', 'Day -15', 'Day -10', 'Day -5', 'Day 0', 'Day 5']):
        ax = fig.add_subplot(6, 1, i + 1, projection=crs.PlateCarree(central_longitude=180))
        ax.outline_patch.set_edgecolor('None')
        ax.add_feature(cartopy.feature.COASTLINE, lw=0.5, zorder=11)

        im = ax.contourf(lons, lats, Matriz_t[i, :, :], cmap=cmap, extend='both', levels=np.arange(min_t, max_t, 0.2),transform=crs.PlateCarree())
        #im2 = ax.contour(lons, lats, Matriz_env[i, :, :], extend='both', levels=[1,2,2.5], colors='k', linewidths=1.5,transform=crs.PlateCarree())
        im2 = ax.contour(lons, lats, Matriz_env[i, :, :], extend='both', levels=levels_env, colors='k', linewidths=1.5,transform=crs.PlateCarree())
        ax.clabel(im2, inline=True, fontsize=10, fmt='%1.1f')
        im3 = ax.contour(lons, lats, Matriz_sf[i, :, :], extend='both', levels=levels_sf, colors='limegreen', negative_linestyles = 'dashed', linewidths=1.7,transform=crs.PlateCarree())
        #ax.clabel(im3, inline=True, fontsize=10, fmt='%1.1f')

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


def hovmoller(time_lags, lons, Matriz_t, cmap_t, norm_t, min_t, max_t, Matriz_sf, cmap_sf, norm_sf, min_sf, max_sf, Matriz_env, levels_env, path, var_1 = '', var_2 = ''):
    import matplotlib.colors
    RdYlBu_list = ['rgb(165,0,38)','rgb(215,48,39)','rgb(244,109,67)','rgb(253,174,97)','rgb(254,224,144)','rgb(255,255,191)','rgb(224,243,248)','rgb(171,217,233)','rgb(116,173,209)','rgb(69,117,180)','rgb(49,54,149)']
    my_cmap = matplotlib.colors.ListedColormap(RdYlBu_list, name='RdYlBu')

    fig = plt.figure(figsize=[7,4])
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
    cb.ax.tick_params(labelcolor='dimgrey', color='dimgrey', labelsize=7)

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


# ================================= DETECTION OF HEAT WAVES ===============================================

# To calculate the duration and the position of the first day of each event
def duration_heat_waves(pos_hw, min_duration):
    count = 1
    duration_hw = []
    pos_day1_hw = []
    for i in range(len(pos_hw[:-1])):
        if pos_hw[i] + 1 == pos_hw[i + 1]:
            count += 1
        else:
            if count >= min_duration and len(pos_day1_hw) == 0:
                duration_hw.append(count)
                pos_day1_hw.append(pos_hw[i - count + 1])
            elif count >= min_duration and abs((pos_day1_hw[-1] + duration_hw[-1]) - pos_hw[i - count + 1]) > 20:
                duration_hw.append(count)
                pos_day1_hw.append(pos_hw[i - count + 1])
            count = 1
    return duration_hw, pos_day1_hw


def anomalies_seasons(df_VAR):
    ANOMA = df_VAR * np.nan
    for i in np.arange(1, 13):
        mes = df_VAR.loc[df_VAR.index.month == i]
        temp = mes * np.nan
        Nodays = mes.index.day.max()
        if np.isnan(Nodays) == True: continue
        for j in np.arange(1, Nodays + 1):
            dia = mes.loc[mes.index.day == j]
            media = dia.mean()
            anoma = dia - media
            temp[temp.index.day == j] = anoma
        ANOMA.loc[temp.index] = temp
    return ANOMA


def anomalies_noseasons(df_VAR):
    ANOMA = df_VAR * np.nan
    media = df_VAR.mean()
    for i in range(df_VAR.shape[0]):
        dia = df_VAR.loc[i]
        anoma = dia - media
        ANOMA.iloc[i] = anoma
    return ANOMA


# ================================= EVOLUTION OF HEAT WAVES ===============================================

def calculate_composites2(pos_HW, matriz):
    dic_composites = {}

    time_lags = np.arange(-20, 21, 1)
    for pos in pos_HW.index:
        if pos == 0: 
            dic_composites[0] = []
            dic_composites[0].append(pos_HW.iloc[0][0])
            continue
        elif pos_HW.iloc[pos][0] >= 27270: break  # len(other variables)
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


def subseasonal_anomalies(df_VAR):
    years = np.unique(np.array([ii.year for ii in df_VAR.index]))
    ANOMA = df_VAR * np.nan
    for i in years:
        pos_y = np.where(df_VAR.index.year == i)[0]
        year = df_VAR.iloc[pos_y]
        temp = year * np.nan
        for j in ([12,1,2], [3,4,5], [6,7,8], [9,10,11]):
            pos_season = np.where((pd.to_datetime(year.index).month == j[0]) ^ (pd.to_datetime(year.index).month == j[1]) ^ (pd.to_datetime(year.index).month == j[2]))[0]
            season = year.iloc[pos_season]
            media = season.mean()
            anoma = season - media
            temp.iloc[pos_season] = anoma
        ANOMA.iloc[pos_y] = temp
    return ANOMA



def hilbert_filtered(x, k_min, k_max, N=None, axis=-1, path=''):
    from scipy import linalg, fft as sp_fft
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


def calculate_composites2(pos_HW, matriz):
    dic_composites = {}

    time_lags = np.arange(-20, 21, 1)
    for pos in pos_HW.index:
        if pos == 0: 
            dic_composites[0] = []
            dic_composites[0].append(pos_HW.iloc[0][0])
            continue
        elif pos_HW.iloc[pos][0] >= 27270: break  # len(other variables)
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
