##########################################################################################################
# functions to estimate beach slope
##########################################################################################################

import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colorbar
from matplotlib import lines
from scipy import integrate as sintegrate
from scipy import signal as ssignal
from astropy.stats import LombScargle
import geopandas as gpd
import pdb

# plotting params
plt.style.use('default')
plt.rcParams['font.size'] = 14
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['legend.fontsize'] = 1

def remove_duplicates(output):
    """
    Function to remove from the output dictionnary entries containing shorelines for the same date
    and satellite mission. This happens when there is an overlap between adjacent
    satellite images.

    Arguments:
    -----------
        output: dict
            contains output dict with shoreline and metadata

    Returns:
    -----------
        output_no_duplicates: dict
            contains the updated dict where duplicates have been removed

    """

    # nested function
    def duplicates_dict(lst):
        "return duplicates and indices"
        def duplicates(lst, item):
                return [i for i, x in enumerate(lst) if x == item]
        return dict((x, duplicates(lst, x)) for x in set(lst) if lst.count(x) > 1)

    dates = output['dates']
    # make a list with year/month/day
    dates_str = [_.strftime('%Y%m%d') for _ in dates]
    # create a dictionnary with the duplicates
    dupl = duplicates_dict(dates_str)
    # if there are duplicates, only keep the first element
    if dupl:
        output_no_duplicates = dict([])
        idx_remove = []
        for k,v in dupl.items():
            idx_remove.append(v[0])
        idx_remove = sorted(idx_remove)
        idx_all = np.linspace(0, len(dates_str)-1, len(dates_str))
        idx_keep = list(np.where(~np.isin(idx_all,idx_remove))[0])
        for key in output.keys():
            output_no_duplicates[key] = [output[key][i] for i in idx_keep]
        print('%d duplicates' % len(idx_remove))
        return output_no_duplicates
    else:
        print('0 duplicates')
        return output

def remove_inaccurate_georef(output, accuracy):
    """
    Function to remove from the output dictionnary entries containing shorelines that were mapped
    on images with inaccurate georeferencing (RMSE > accuracy or flagged with -1)

    Arguments:
    -----------
        output: dict
            contains the extracted shorelines and corresponding metadata
        accuracy: int
            minimum horizontal georeferencing accuracy (metres) for a shoreline to be accepted

    Returns:
    -----------
        output_filtered: dict
            contains the updated dictionnary

    """

    # find indices of shorelines to be removed
    idx = np.where(~np.logical_or(np.array(output['geoaccuracy']) == -1,
                                  np.array(output['geoaccuracy']) >= accuracy))[0]
    output_filtered = dict([])
    for key in output.keys():
        output_filtered[key] = [output[key][i] for i in idx]
    print('%d bad georef' % (len(output['geoaccuracy']) - len(idx)))
    return output_filtered

def transects_from_geojson(filename):
    """
    Reads transect coordinates from a .geojson file.

    Arguments:
    -----------
        filename: str
            contains the path and filename of the geojson file to be loaded

    Returns:
    -----------
        transects: dict
            contains the X and Y coordinates of each transect.

    """

    gdf = gpd.read_file(filename)
    transects = dict([])
    for i in gdf.index:
        transects[gdf.loc[i,'name']] = np.array(gdf.loc[i,'geometry'].coords)

    print('%d transects have been loaded' % len(transects.keys()))

    return transects

def compute_intersection(output, transects, settings):
    """
    Computes the intersection between the 2D mapped shorelines and the transects, to generate
    time-series of cross-shore distance along each transect.

    Arguments:
    -----------
        output: dict
            contains the extracted shorelines and corresponding dates.
        transects: dict
            contains the X and Y coordinates of the transects (first and last point needed for each
            transect).
        settings: dict
                along_dist: float
                alongshore distance to caluclate the intersection (median of points
                within this distance).
                max_std: float
                if the standard deviation of the points is above this threshold a nan is returned
                max_range: float
                if the range of the points is above this threshold a nan is returned                
                min_val: float
                largest negative value along transect (landwards of transect origin)
                nan/max: str
                'nan', 'max' or 'auto', how to deal with multiple intersections, 
                either put a nan or take the maximum (most seawards intersection),
                or automatically decide based on the occurence of multiple intersections
                (for example if there is a lagoon behind the beach, there are always 2 intersections)
                prc_std: percentage of occurrence to use in 'auto' mode to switch from 'nan' to 'max'

    Returns:
    -----------
        cross_dist: dict
            time-series of cross-shore distance along each of the transects. These are not tidally
            corrected.

    """

    # initialise dictionary with intersections for each transect
    cross_dist = dict([])

    shorelines = output['shorelines']
    along_dist = settings['along_dist']

    # loop through each transect
    for key in transects.keys():

        # initialise variables
        std_intersect = np.zeros(len(shorelines))
        med_intersect = np.zeros(len(shorelines))
        max_intersect = np.zeros(len(shorelines))
        min_intersect = np.zeros(len(shorelines))
        n_intersect = np.zeros(len(shorelines))

        # loop through each shoreline
        for i in range(len(shorelines)):

            sl = shorelines[i]

            # compute rotation matrix
            X0 = transects[key][0,0]
            Y0 = transects[key][0,1]
            temp = np.array(transects[key][-1,:]) - np.array(transects[key][0,:])
            phi = np.arctan2(temp[1], temp[0])
            Mrot = np.array([[np.cos(phi), np.sin(phi)],[-np.sin(phi), np.cos(phi)]])

            # calculate point to line distance between shoreline points and the transect
            p1 = np.array([X0,Y0])
            p2 = transects[key][-1,:]
            d_line = np.abs(np.cross(p2-p1,sl-p1)/np.linalg.norm(p2-p1))
            # calculate the distance between shoreline points and the origin of the transect
            d_origin = np.array([np.linalg.norm(sl[k,:] - p1) for k in range(len(sl))])
            # find the shoreline points that are close to the transects and to the origin
            # the distance to the origin is hard-coded here to 1 km
            idx_dist = np.logical_and(d_line <= along_dist, d_origin <= 1000)
            idx_close = np.where(idx_dist)[0]

            # in case there are no shoreline points close to the transect
            if len(idx_close) == 0:
                std_intersect[i] = np.nan
                med_intersect[i] = np.nan
                max_intersect[i] = np.nan
                min_intersect[i] = np.nan
                n_intersect[i] = np.nan
            else:
                # change of base to shore-normal coordinate system
                xy_close = np.array([sl[idx_close,0],sl[idx_close,1]]) - np.tile(np.array([[X0],
                                   [Y0]]), (1,len(sl[idx_close])))
                xy_rot = np.matmul(Mrot, xy_close)
                # remove points that are too far landwards relative to the transect origin (i.e., negative chainage)
                xy_rot[0, xy_rot[0,:] < settings['min_val']] = np.nan

                # compute std, median, max, min of the intersections
                std_intersect[i] = np.nanstd(xy_rot[0,:])
                med_intersect[i] = np.nanmedian(xy_rot[0,:])
                max_intersect[i] = np.nanmax(xy_rot[0,:])
                min_intersect[i] = np.nanmin(xy_rot[0,:])
                n_intersect[i] = len(xy_rot[0,:])

        # quality control the intersections using dispersion metrics (std and range)
        condition1 = std_intersect <= settings['max_std']
        condition2 = (max_intersect - min_intersect) <= settings['max_range']
        condition3 = n_intersect > 2
        idx_good = np.logical_and(np.logical_and(condition1, condition2), condition3)

        # decide what to do with the intersections with high dispersion
        if settings['nan/max'] == 'auto':
            # compute the percentage of data points where the std is larger than the user-defined max
            prc_over = np.sum(std_intersect > settings['max_std'])/len(std_intersect)
            # if more than a certain percentage is above, use the maximum intersection
            if prc_over > settings['prc_std']:
                med_intersect[~idx_good] = max_intersect[~idx_good]
                med_intersect[~condition3] = np.nan
            # otherwise put a nan
            else:
                med_intersect[~idx_good] = np.nan

        elif settings['nan/max'] == 'max':
            med_intersect[~idx_good] = max_intersect[~idx_good]
            med_intersect[~condition3] = np.nan

        elif settings['nan/max'] == 'nan':
            med_intersect[~idx_good] = np.nan

        else:
            raise Exception('the nan/max parameter can only be: nan, max or auto')

        # store in dict
        cross_dist[key] = med_intersect

    return cross_dist

def reject_outliers(cross_distance, output, settings):
    """

    Arguments:
    -----------
        cross_distance: dict
            time-series of shoreline change
        output: dict
            mapped shorelines with metadata
        settings: dict

    Returns:
    -----------
        chain_dict: dict
            contains the updated time-series of cross-shore distance with the corresponding dates

    """

    chain_dict = dict([])

    for i,key in enumerate(list(cross_distance.keys())):

        chainage = cross_distance[key].copy()
        if sum(np.isnan(chainage)) == len(chainage):
            continue

        # 1. remove nans and negative chainages
        idx_nonan = np.where(~np.isnan(chainage))[0]
        chainage1 = [chainage[k] for k in idx_nonan]
        dates1 = [output['dates'][k] for k in idx_nonan]

        # 3. remove outliers based on despiking [iterative method]
        chainage3, dates3 = identify_outliers(chainage1, dates1, settings['max_cross_change'])

        # fill with nans the indices to be removed from cross_distance
        idx_kept = []
        for date in output['dates']: idx_kept.append(date in dates3)
        chainage[~np.array(idx_kept)] = np.nan
        # store in chain_dict
        chain_dict[key] = chainage

        print('%s  - outliers removed %d'%(key, len(dates1) - len(dates3)))

    return chain_dict

def identify_outliers(chainage, dates, cross_change):
    """
    Remove outliers based on despiking [iterative method]

    Arguments:
    -----------
    chainage: list
        time-series of shoreline change
    dates: list of datetimes
        correspondings dates
    cross_change: float
        threshold distance to identify a point as an outlier

    Returns:
    -----------
    chainage_temp: list
        time-series of shoreline change without outliers
    dates_temp: list of datetimes
        dates without outliers

    """

    # make a copy of the inputs
    chainage_temp = chainage.copy()
    dates_temp = dates.copy()

    # loop through the time-series always starting from the start
    # when an outlier is found, remove it and restart
    # repeat until no more outliers are found in the time-series
    k = 0
    while k < len(chainage_temp):

        for k in range(len(chainage_temp)):

            # check if the first point is an outlier
            if k == 0:
                # difference between 1st and 2nd point in the time-series
                diff = chainage_temp[k] - chainage_temp[k+1]
                if np.abs(diff) > cross_change:
                    chainage_temp.pop(k)
                    dates_temp.pop(k)
                    break

            # check if the last point is an outlier
            elif k == len(chainage_temp)-1:
                # difference between last and before last point in the time-series
                diff = chainage_temp[k] - chainage_temp[k-1]
                if np.abs(diff) > cross_change:
                    chainage_temp.pop(k)
                    dates_temp.pop(k)
                    break

            # check if a point is an isolated outlier or in a group of 2 consecutive outliers
            else:
                # calculate the difference with the data point before and after
                diff_m1 = chainage_temp[k] - chainage_temp[k-1]
                diff_p1 = chainage_temp[k] - chainage_temp[k+1]
                # remove point if isolated outlier, distant from both neighbours
                condition1 = np.abs(diff_m1) > cross_change
                condition2 = np.abs(diff_p1) > cross_change
                # check that distance from neighbours has the same sign
                condition3 = np.sign(diff_p1) == np.sign(diff_m1)
                if np.logical_and(np.logical_and(condition1,condition2),condition3):
                    chainage_temp.pop(k)
                    dates_temp.pop(k)
                    break

                # check for 2 consecutive outliers in the time-series
                if k >= 2 and k < len(chainage_temp)-2:

                    # calculate difference with the data around the neighbours of the point
                    diff_m2 = chainage_temp[k-1] - chainage_temp[k-2]
                    diff_p2 = chainage_temp[k+1] - chainage_temp[k+2]
                    # remove if there are 2 consecutive outliers (see conditions below)
                    condition4 = np.abs(diff_m2) > cross_change
                    condition5 = np.abs(diff_p2) > cross_change
                    condition6 = np.sign(diff_m1) == np.sign(diff_p2)
                    condition7 = np.sign(diff_p1) == np.sign(diff_m2)
                    # check for both combinations (1,5,6 and ,2,4,7)
                    if np.logical_and(np.logical_and(condition1,condition5),condition6):
                        chainage_temp.pop(k)
                        dates_temp.pop(k)
                        break
                    elif np.logical_and(np.logical_and(condition2,condition4),condition7):
                        chainage_temp.pop(k)
                        dates_temp.pop(k)
                        break

                    # also look for clusters of 3 outliers
                    else:
                        # increase the distance to make sure these are really outliers
                        condition4b = np.abs(diff_m2) > 1.5*cross_change
                        condition5b = np.abs(diff_p2) > 1.5*cross_change
                        condition8 = np.sign(diff_m2) == np.sign(diff_p2)
                        # if point is close to immediate neighbours but
                        # the neighbours are far from their neighbours, point is an outlier
                        if np.logical_and(np.logical_and(np.logical_and(condition4b,condition5b),
                                                         np.logical_and(~condition1,~condition2)),
                                                         condition8):
                            print('*', end='')
                            chainage_temp.pop(k)
                            dates_temp.pop(k)
                            break

        # if one full loop is completed (went through all the time-series without removing outlier)
        # then increment k to get out of the loop
        k = k + 1


    # return the time-series where the outliers have been removed
    return chainage_temp, dates_temp

def plot_cross_distance(dates, cross_distance):
    'plot the time-series of shoreline change from CoastSat'

    for i,key in enumerate(cross_distance.keys()):
        idx_nan = np.isnan(cross_distance[key])
        chain = cross_distance[key][~idx_nan]
        dates_temp = [dates[k] for k in np.where(~idx_nan)[0]]
        if len(chain)==0 or sum(idx_nan) > 0.5*len(idx_nan): continue
        fig,ax=plt.subplots(1,1,figsize=[12,3])
        fig.set_tight_layout(True)
        ax.grid(linestyle=':', color='0.5')
        ax.plot(dates_temp, chain - np.mean(chain), '-o', ms=3, mfc='w', mec='C0')
        ax.set(title='Transect %s - %d points'%(key,len(chain)), ylabel='distance [m]', ylim=get_min_max_dict(cross_distance))

def get_min_max_dict(cross_distance):
    'get min and max of a dictionary of time-series'
    xmin = 1e10
    xmax = -1e10
    for key in cross_distance.keys():
        ts = cross_distance[key] - np.nanmedian(cross_distance[key])
        if np.nanmin(ts) < xmin:
            xmin = np.nanmin(ts)
        if np.nanmax(ts) > xmax:
            xmax = np.nanmax(ts)
    xmax = np.max([np.abs(xmin),np.abs(xmax)])
    xmin = -np.max([np.abs(xmin),np.abs(xmax)])
    return [xmin, xmax]

def get_min_max(y):
    'get min and max of a time-series'
    ymin = np.nanmin(y)
    ymax = np.nanmax(y)
    ymax = np.max([np.abs(ymin),np.abs(ymax)])
    ymin = -np.max([np.abs(ymin),np.abs(ymax)])
    return [ymin,ymax]

###################################################################################################
# Tide functions
###################################################################################################

def compute_tide(coords,date_range,time_step,ocean_tide,load_tide):
    'compute time-series of water level for a location and dates using a time_step'
    # list of datetimes (every timestep)
    dates = []
    date = date_range[0]
    while date <= date_range[1]:
        dates.append(date)
        date = date + timedelta(seconds=time_step)
    # convert list of datetimes to numpy dates
    dates_np = np.empty((len(dates),), dtype='datetime64[us]')
    for i,date in enumerate(dates):
        dates_np[i] = datetime(date.year,date.month,date.day,date.hour,date.minute,date.second)
    lons = coords[0]*np.ones(len(dates))
    lats = coords[1]*np.ones(len(dates))
    # compute heights for ocean tide and loadings
    ocean_short, ocean_long, min_points = ocean_tide.calculate(lons, lats, dates_np)
    load_short, load_long, min_points = load_tide.calculate(lons, lats, dates_np)
    # sum up all components and convert from cm to m
    tide_level = (ocean_short + ocean_long + load_short + load_long)/100

    return dates, tide_level

def compute_tide_dates(coords,dates,ocean_tide,load_tide):
    'compute time-series of water level for a location and dates (using a dates vector)'
    dates_np = np.empty((len(dates),), dtype='datetime64[us]')
    for i,date in enumerate(dates):
        dates_np[i] = datetime(date.year,date.month,date.day,date.hour,date.minute,date.second)
    lons = coords[0]*np.ones(len(dates))
    lats = coords[1]*np.ones(len(dates))
    # compute heights for ocean tide and loadings
    ocean_short, ocean_long, min_points = ocean_tide.calculate(lons, lats, dates_np)
    load_short, load_long, min_points = load_tide.calculate(lons, lats, dates_np)
    # sum up all components and convert from cm to m
    tide_level = (ocean_short + ocean_long + load_short + load_long)/100

    return tide_level

def find_tide_peak(dates,tide_level,settings):
    'find the high frequency peak in the tidal time-series'
    # create frequency grid
    t = np.array([_.timestamp() for _ in dates]).astype('float64')
    days_in_year = 365.2425
    seconds_in_day = 24*3600
    time_step = settings['n_days']*seconds_in_day
    freqs = frequency_grid(t,time_step,settings['n0'])
    # compute power spectrum
    ps_tide,_,_ = power_spectrum(t,tide_level,freqs,[])
    # find peaks in spectrum
    idx_peaks,_ = ssignal.find_peaks(ps_tide, height=0)
    y_peaks = _['peak_heights']
    idx_peaks = idx_peaks[np.flipud(np.argsort(y_peaks))]
    # find the strongest peak at the high frequency (defined by freqs_cutoff[1])
    idx_max = idx_peaks[freqs[idx_peaks] > settings['freqs_cutoff']][0]
    # compute the frequencies around the max peak with some buffer (defined by buffer_coeff)
    freqs_max = [freqs[idx_max] - settings['delta_f'], freqs[idx_max] + settings['delta_f']]
    # make a plot of the spectrum
    fig = plt.figure()
    fig.set_size_inches([12,4])
    fig.set_tight_layout(True)
    ax = fig.add_subplot(111)
    ax.grid(linestyle=':', color='0.5')
    ax.plot(freqs,ps_tide)
    ax.set_title('$\Delta t$ = %d days'%settings['n_days'], x=0, ha='left')
    ax.set(xticks=[(days_in_year*seconds_in_day)**-1, (30*seconds_in_day)**-1, (16*seconds_in_day)**-1, (8*seconds_in_day)**-1],
                   xticklabels=['1y','1m','16d','8d']);
    # show top 3 peaks
    for k in range(2):
        ax.plot(freqs[idx_peaks[k]], ps_tide[idx_peaks[k]], 'ro', ms=4)
        ax.text(freqs[idx_peaks[k]], ps_tide[idx_peaks[k]]+1, '%.1f d'%((freqs[idx_peaks[k]]**-1)/(3600*24)),
                ha='center', va='bottom', fontsize=8, bbox=dict(boxstyle='square', ec='k',fc='w', alpha=0.5))
    ax.axvline(x=freqs_max[1], ls='--', c='0.5')
    ax.axvline(x=freqs_max[0], ls='--', c='0.5')
    ax.axvline(x=(2*settings['n_days']*seconds_in_day)**-1, ls='--', c='k')
    return freqs_max

def frequency_grid(time,time_step,n0):
    'define frequency grid for Lomb-Scargle transform'
    T = np.max(time) - np.min(time)
    fmin = 1/T
    fmax = 1/(2*time_step) # Niquist criterium
    df = 1/(n0*T)
    N = np.ceil((fmax - fmin)/df).astype(int)
    freqs = fmin + df * np.arange(N)
    return freqs

def power_spectrum(t,y,freqs,idx_cut):
    'compute power spectrum and integrate'
    model = LombScargle(t, y, dy=None, fit_mean=True, center_data=True, nterms=1, normalization='psd')
    ps = model.power(freqs)
    # integrate the entire power spectrum
    E = sintegrate.simps(ps, x=freqs, even='avg')
    if len(idx_cut) == 0:
        idx_cut = np.ones(freqs.size).astype(bool)
    # integrate only frequencies above cut-off
    Ec = sintegrate.simps(ps[idx_cut], x=freqs[idx_cut], even='avg')
    return ps, E, Ec

###################################################################################################
# Slope functions
###################################################################################################

def range_slopes(min_slope, max_slope, delta_slope):
    'create list of beach slopes to test'
    beach_slopes = []
    slope = min_slope
    while slope < max_slope:
        beach_slopes.append(slope)
        slope = slope + delta_slope
    beach_slopes.append(slope)
    beach_slopes = np.round(beach_slopes,len(str(delta_slope).split('.')[1]))
    return beach_slopes

def tide_correct(chain,tide_level,beach_slopes):
    'apply tidal correction with a range of slopes'
    tsall = []
    for i,slope in enumerate(beach_slopes):
        # apply tidal correction
        tide_correction = (tide_level)/slope
        ts = chain + tide_correction
        tsall.append(ts)
    return tsall

def plot_spectrum_all(dates_rand,composite,tsall,settings, title):
    'plot the spectrum of the tidally-corrected time-series of shoreline change'
    t = np.array([_.timestamp() for _ in dates_rand]).astype('float64')
    seconds_in_day = 24*3600
    days_in_year = 365.2425
    time_step = settings['n_days']*seconds_in_day
    freqs = frequency_grid(t,time_step,settings['n0'])
    beach_slopes = range_slopes(settings['slope_min'], settings['slope_max'], settings['delta_slope'])

    # make figure 1
    fig = plt.figure()
    fig.set_size_inches([12,5])
    fig.set_tight_layout(True)
    fig.suptitle(title, x=0.1, ha='left', fontweight='bold',
     bbox=dict(boxstyle='square', ec='k',fc='w', alpha=0.5))
    cmap = cm.get_cmap('RdYlGn')
    color_list = cmap(np.linspace(0,1,len(beach_slopes)))
    indices = np.arange(0,len(beach_slopes))
    # axis labels
    freq_1month = 1/(days_in_year*seconds_in_day/12)
    xt = freq_1month/np.array([12,3,1, 20/(days_in_year/12), 16/(days_in_year/12)])
    xl = ['1y', '3m', '1m', '20d', '16d']
    # loop for plots
    ax = fig.add_subplot(111)
    ax.grid(which='major', linestyle=':', color='0.5')
    ax.set(xticks=xt, xticklabels=xl, title='Power Spectrum of tidally-corrected time-series', ylabel='amplitude')
    for i,idx in enumerate(indices):
        # compute spectrum
        ps,_,_ = power_spectrum(t,tsall[idx],freqs,[])
        ax.plot(freqs, ps, '-', color=color_list[idx,:], lw=1)
    # draw some references
    ax.axvline(x=settings['freqs_max'][0], ls='--', c='0.5')
    ax.axvline(x=settings['freqs_max'][1], ls='--', c='0.5')
    ax.axvspan(xmin=settings['freqs_max'][0], xmax=settings['freqs_max'][1], color='0.85')
    ax.axvline(x=(16*seconds_in_day)**-1, ls='--', c='k')

    # make figure 2
    fig = plt.figure()
    fig.set_size_inches([12,5])
    fig.set_tight_layout(True)
    # axis labels
    xt = 1./(np.flipud(np.arange(settings['n_days']*2,21,1))*24*3600)
    xl = ['%d d'%(_) for _ in np.flipud(np.arange(settings['n_days']*2,21,1))]
    # loop for plots
    ax = fig.add_subplot(111)
    ax.axvline(x=settings['freqs_max'][0], ls='--', c='0.5')
    ax.axvline(x=settings['freqs_max'][1], ls='--', c='0.5')
    ax.axvspan(xmin=settings['freqs_max'][0], xmax=settings['freqs_max'][1], color='0.85')
    ax.grid(which='major', linestyle=':', color='0.5')
    ax.set(xticks=xt, xticklabels=xl, ylabel='amplitude', title='Inset into the tidal peak frequency bands')
    idx_interval = np.logical_and(freqs >= settings['freqs_max'][0], freqs <= settings['freqs_max'][1])
    for i,idx in enumerate(indices):
        # compute spectrum
        ps, _,_ = power_spectrum(t,tsall[idx],freqs,[])
        ax.plot(freqs[idx_interval], ps[idx_interval], '-', color=color_list[idx,:], lw=1)
    # non-corrected time-series
    ps,_,_ = power_spectrum(t,composite,freqs,[])
    ax.plot(freqs[idx_interval], ps[idx_interval], '--', color='b', lw=1.5)
    # add legend
    nc_line = lines.Line2D([],[],ls='--', c='b', lw=1.5, label='non-corrected time-series')
    ax.legend(handles=[nc_line], loc=2)

def integrate_power_spectrum(dates_rand,tsall,settings):
    'integrate power spectrum at the frequency band of peak tidal signal'
    t = np.array([_.timestamp() for _ in dates_rand]).astype('float64')
    seconds_in_day = 24*3600
    time_step = settings['n_days']*seconds_in_day
    freqs = frequency_grid(t,time_step,settings['n0'])
    beach_slopes = range_slopes(settings['slope_min'], settings['slope_max'], settings['delta_slope'])
    # integrate power spectrum
    idx_interval = np.logical_and(freqs >= settings['freqs_max'][0], freqs <= settings['freqs_max'][1])
    E = np.zeros(beach_slopes.size)
    for i in range(len(tsall)):
        ps, _, _ = power_spectrum(t,tsall[i],freqs,[])
        E[i] = sintegrate.simps(ps[idx_interval], x=freqs[idx_interval], even='avg')
    # plot energy vs slope curve
    fig = plt.figure()
    fig.set_size_inches([12,4])
    fig.set_tight_layout(True)
    ax = fig.add_subplot(111)
    ax.grid(linestyle=':', color='0.5')
    ax.set(title='Energy in tidal frequency band', xlabel='slope values',ylabel='energy')
    ax.plot(beach_slopes,E,'-k',lw=1.5)
    cmap = cm.get_cmap('RdYlGn')
    color_list = cmap(np.linspace(0,1,len(beach_slopes)))
    for i in range(len(beach_slopes)): ax.plot(beach_slopes[i], E[i],'o',ms=8,mec='k',mfc=color_list[i,:])
    ax.plot(beach_slopes[np.argmin(E)],np.min(E),'bo',ms=14,mfc='None',mew=2, label='%.3f'%beach_slopes[np.argmin(E)])

    return beach_slopes[np.argmin(E)]
