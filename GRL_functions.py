# function to read in all the model files, **note it assumes a time period from 1850-2014 with 1980 timesteps
def read_models(institution_dir, variable_dir, start_date, end_date):
    """ Read in all the CMIP6 histtorical run model files of specified variable from specified start and end date.
        **note it assumes a time period from 1850-2014 with 1980 timesteps 
        
        Args:
        institution_directory (str): directory of CMIP6 institutions
        variable_directory (str): climate variable (e.g. tas, pr, sst, ts etc.)
        start_date (date_str): start date for data
        end_date (date_str): end date for data
    """
    
    # import relevant libraries
    import os
    import xarray as xr
    import pandas as pd
    # store all institutions found in the institution directory
    institution_list = os.listdir(f'{institution_dir}')
    
    # creates a dictionary containing the model and model path
    models = {}

    # find the models from each instituion and store them in a list
    for institution in institution_list:
        model_list = os.listdir(f'{institution_dir}{institution}')
        # find the 
        for model in model_list:
            # check if the historical model with the right variable exists and if so save the version number for the file
            if os.path.exists(f'{institution_dir}{institution}/{model}{variable_dir}'):
                version = os.listdir(f'{institution_dir}{institution}/{model}{variable_dir}')
                # for each version, save the path and then store with the model in a dictionary 'models'
                for v in version:
                    path = f'{institution_dir}{institution}/{model}{variable_dir}{v}/*.nc'
                    models[model] = path
                     
    # Prints the number of models loaded into the dictionary
    print(f'{len(models)} model paths found and loaded into the dictionary "models"') 
    
    # Now open each dataset, and store the dataset and model name in two arrays, 'names' and 'ds'.  
    # The try-except method allows all other files to be read even if one file does not exist or has some issues. 
    # By deleting the time and bnds coordinates we are removing any differences in time.  
    # (We add the time back later). 
    names = []
    ds = []

    for name, path in models.items():
        try:
            d = xr.open_mfdataset(path, combine='by_coords', chunks={'time': 12})
            # checks if there is data for each month in the time period, and if so stores the model file data
            if len(d['time'])==1980:
                # remove any differences in time
                del d['time_bnds']
                del d['time']
                # select times from the specified period
                time_month = pd.date_range(start = f'{start_date}',end = f'{end_date}', freq ='M')
                # add the time coordinate back in
                d.coords['time'] = time_month
                ds.append(d)
                names.append(name)
                #print(name, path)
            else:
                print(f'Model {name} has different time so is now removed')
        except OSError:
            # No files read, move on to the next
            print(f'Path for {name} does not exist')
            continue

    # Combine the individual data sets into a single xarray dataset, with the coordinate
    # 'model' representing the source model
    multi_model = xr.concat(ds, dim='model', coords = 'minimal')
    multi_model.coords['model'] = names
    
    # del multi_model.model['height']
                
    # print the number of models that have been successfully loaded into the xarray
    print(f'{len(multi_model.model)} models have been successfully loaded into an xarray')
    
    return multi_model


def SEA_combine(K,S,A,E,P,names):
    """ Reset the time axis for data used in the superposed epoch analysis  
        
        Args:
        K (xarray): data set of first eruption 
        S (xarray): data set of second eruption 
        A (xarray): data set of third eruption 
        E (xarray): data set of fourth eruption 
        P (xarray): data set of fifth eruption 
        names (dict): dictionary of names for each of the five eruptions
    """
    import numpy as np
    import xarray as xr
    ds=[]
    ds.append(K)
    ds.append(S)
    ds.append(A)
    ds.append(E)
    ds.append(P)

    # store all eruptions in an array
    composite = xr.concat(ds, dim='volcano', coords = 'minimal')
    composite.coords['volcano'] = names
    
    return composite


# define function to calculate monthly anomalies for a multidimensional array of models
def monthly_anomaly(dataset, start_date, end_date):
    
    """ Calculate monthly anomalies for a multidimensional array of models.  
        
        Args:
        dataset (xarray): data set of climate variable (e.g tas)
        start_date (date_str): start date of climatology to calculate monthly anomaly
        end_date (date_str): end date of climatology to calculate monthly anomaly
    """
    
    # group the data into months
    variable_monthly = dataset.groupby('time.month')

    # calculate the mean climatology along each month for the time period 1850-1900 
    clim_monthly = dataset.sel(time = slice(f'{start_date}', f'{end_date}')).groupby('time.month').mean(dim = 'time')

    # caclulate the anomalies for each month and return it as an array
    multi_monthly_anom = (variable_monthly - clim_monthly)

    return multi_monthly_anom


# define function to calculate the seasonal mean used in seasonal anomaly calculation:
def seasonal_mean(data):
    """ Calculate the seasonal mean used in seasonal anomaly calculation.  
        
        Args:
        data (xarray): data set of climate variable (e.g tas)
    """
    return data.groupby('time.season').mean()


# function to calculate a seasonal anomaly for a multidimensional xarray over a time period entered by user
def seasonal_anomaly(dataset, start_date, end_date):
    """ Calculate a seasonal anomaly for a multidimensional xarray over a time period entered by user.  
        
        Args:
        dataset (xarray): data set of climate variable (e.g tas)
        start_date (date_str): start date to calculate seasonal anomaly
        end_date (date_str): end date to calculate seasonal anomaly
    """
    # first I need to define a new coordinate (seasonyear) so that december gets counted with the adjoining jan and feb
    seasonyear = (dataset.time.dt.year + (dataset.time.dt.month//12)) 
    dataset.coords['seasonyear'] = seasonyear
    
        
    # group data into seasons and calculate the seasonal mean for each year in the dataset 
    yearly_seasonal = dataset.groupby('seasonyear').apply(seasonal_mean)

    # calculate the mean climatology along each season for the time period 
    clim_seasonal = yearly_seasonal.sel(seasonyear = slice(f'{start_date}',f'{end_date}')).mean(dim = 'seasonyear')

    # calculate the anomaly and returns it as an xarray
    multi_seasonal_anom = (yearly_seasonal - clim_seasonal)
        
    return multi_seasonal_anom


# defines an array of titles for seasonal spatial graphs - shorter format to seasonal_title(2)
def seasonal_title_short(K_dates, season_name, season):
    """Create titles for graphs by combining strings for each year, season post-eruption.  
    
    Args:
        K_dates (list): list of years to be plotted
        title_label (list): list of summer relative to eruption
        season_name(list): season name (e.g. 'summer')
        season (list): season (e.g. 'DJF')
    """
    title_label = [f'pre-eruption', f'1st {season_name}', 
                   f'2nd {season_name}', f'3rd {season_name}', 
                   f'4th {season_name}']
    
    titles=[]
    for i,vals in enumerate(K_dates):
        t = title_label[i] + ' (' + season + f' {K_dates[i]})'
        titles.append(t)
    
    return titles


# find where the anomalies are outside a threshold of +/- 2 standard deviations 
def stat_sig(dataset):
    """Find where the anomalies are outside a threshold of +/- 2 standard deviations.  Standard deviation calculated based on an 1850-1880 climatology.  
    
    Args:
        dataset (xarray): xarray of climate variable(s)
    """    
    import xarray as xr
    
    # calculate the standard deviation 
    if hasattr(dataset, 'time'):
        std = dataset.sel(time = slice('1850-01', '1879-12')).std(dim = ['time'])
    elif hasattr(dataset, 'seasonyear'):
        std = dataset.sel(seasonyear = slice('1850', '1879')).std(dim = ['seasonyear'])
    
    # mark points oustide the 2 standard deviation threshold with a 100 (and non significant points with a zero)
    sig = xr.where((dataset < - 2*std) | (dataset > 2*std), 100, 0)
    
    return sig 


# calculate the nino 3.4 index 
def nino34(sst_dataset, start_date, end_date, std):
    """ Calculate the NINO34 index from SST values and normalise by dividing by the standard deviation calculate over user specified time period.   
        
        Args:
        sst_dataset (xarray): data set of sea surface temperature values
        start_date (date_str): start date of std climatology
        end_date (date_str): end date of std climatology
        std (int): if std==1, calculate the std and divide NINO34 index by std
    """
    # select out the region for nino34 definition
    region = sst_dataset.sel(lat=slice(-5,5), lon=slice(190,240))
    
    # calculate the mean climatology along each month
    clim = region.sel(time = slice(f'{start_date}', f'{end_date}')).groupby('time.month').mean(dim = ['time','lat','lon'])
    
    # calculate the anomaly using input dates for climatology and take the lat lon mean 
    #anom = monthly_anom_xr(region, f'{start_date}', f'{start_date}').mean(dim=['lat','lon'])
    anom = (region.groupby('time.month') - clim).mean(dim=['lat','lon'])
    
    # chunk the data into groups of 5 timepoints so I can then use rolling mean 
    anom = anom.chunk({'time': 5})
    
    if std == 1:
        # calculate the standard deviation so we can normalise the model data 
        std = region.sel(time = slice(f'{start_date}', f'{end_date}')).mean(dim=['lat', 'lon']).std(dim = ['time'])
        
        # calculate the nino3.4 index using a rolling 5 month mean and normalised by the std
        nino34_index = anom.rolling(time=5).mean() / std
    elif std == 0:
            nino34_index = anom.rolling(time=5).mean()
    
    return nino34_index



# -------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------- PLOTTING FUNCTIONS ---------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------

# define a function for spatial plots, plotting a dataset at different time intervals
def spatial_plot_cv(rows, cols, dataset, cmax, times, titles, colours, units, std=None, label_loc=None, cbar=None):
    """Create a figure of spatial graphs with subplots for each time snapshot as specified in the dataset and times array. 
    Can use stippling to show areas where anomaly exceeds 2 standard deviations. 
    
    Args:
        rows (int): number of rows of subplots
        cols (int): number of columns of subplots
        dataset (xarray): data array of climate variable to be plotted
        cmax (float): 2 element arry with [minimum, maximum] value for colourbar
        times (date_str): dictionary of dates (date_str) for each time to be plotted
        titles (date_str): dictionary of titles (str) for each subplot
        colours (dict): colour palette for cmap
        units (str): units for axes label
        std (int): if std==1: use stippling
    """
    import matplotlib.pyplot as plt, cartopy.crs as ccrs, numpy as np
    
    fig = plt.figure()
    axs = []
    
    # calculate the standard deviation (if desired, and input to 'std' is 1)
    if std == 1:
        sig_dataset = stat_sig(dataset)
      
    if label_loc is None:
        label_loc = [4.8, 0.15]
        
    # set cbar size and padding
    if cbar is None:
        cbar = [0.35, 0.02]
    
    # set discrete colourbar with 15 intervals
    cmap = plt.get_cmap(f'{colours}')#, 15)
    
    for i, d in enumerate(times):    
        # Add a subplot with a projection    
        ax = fig.add_subplot(rows, cols, i+1, projection=ccrs.PlateCarree(180))        
        # Select the date and corresponding data and plot it    
        # We'll add a separate colour bar, but make sure all plots share the same min and max colour value
        if hasattr(dataset, 'time'):
            data = dataset.sel(time = times[i])   
        elif hasattr(dataset, 'seasonyear'):
            data = dataset.sel(seasonyear = times[i]) 
        
        C = data.plot(ax=ax, add_colorbar=False, transform=ccrs.PlateCarree(), cmap = cmap, vmin=cmax[0], vmax=cmax[1])
        # hatching where anomalies exceed a threshold of 2 standard deviations
        if std == 1:
            if hasattr(dataset, 'time'):
                data2 = sig_dataset.sel(time = times[i]).mean(dim='time')
                data2.plot.contourf(levels=[99, 1e10], hatches=[None,'..'], colors='none', add_colorbar=False, 
                                    transform=ccrs.PlateCarree())
            elif hasattr(dataset, 'seasonyear'): # try this instead if the dataset is seasonal
                data2 = sig_dataset.sel(seasonyear = times[i])
                data2.plot.contourf(levels=[99, 1e10], hatches=[None,'..'], colors='none', add_colorbar=False, 
                                    transform=ccrs.PlateCarree())
        
        # axes
        ax.coastlines()
        # set the axis limits to be slightly larger (2.5 degrees wither way) than the upper and lower bounds of the dataset 
        if (len(data.lon) < int(175/1.5)) & (len(data.lat) < int(175/1.5)):
            ax.set_extent([data.lon[0] - 2.5, data.lon[-1] + 2.5, data.lat[0] - 2.5, data.lat[-1] + 2.5], crs=ccrs.PlateCarree())

        # add titles for each subplot
        ax.set_xlabel(None)
        ax.set_title(titles[i])
        # Gather all axes for when we make the colour bar    
        axs.append(ax)    

    # Put the colour bar to the left of all axes
    cbar = plt.colorbar(C, orientation='vertical', ax=axs, shrink=cbar[0], pad=cbar[1])
    cbar.ax.text(label_loc[0], label_loc[1],f'{units}', fontsize=12, rotation=90)
    
    return fig


# define function to plot figures for composite graphs 
def SEA_plots(mmm_dataset, comp_dataset, p10 = None, p90 = None, color_cycle = None, ax = None, line_style = None, subplot_legend = None, **plt_kwargs):
    """Create subplots for a superposed epoch analysis (SEA) graph.  SEA graph is composed of time series of each eruption contained in the mmm_dataset and the composite (of all eruptions in the mmm_dataset).  Shading is used to show the 10th and 90th percentiles of the composite.   
    Return the axis.  
    
    Args:
        mmm_dataset (xarray): xarray of eruptions and the values (multi-model mean of climate variable) for each to be plotted
        comp_dataset (xarray): xarray of composite values (multi-eruption multi-model mean of climate variable) to be plotted
        color_cycle (dict): dictionary of colours (as strings) 
        p10 (array): array of values of 10th percentile
        p90 (array): array of values of 90th percentile
        ax (axis): axis
        **kwargs
    """
    import xarray as xr, matplotlib.pyplot as plt, numpy as np, seaborn as sns
    
    # checking if an axis has been defined and if not creates one with function "get current axes"
    if ax is None:
        ax = plt.gca()  
        
    if line_style is None:
        line_style = ['-', ':', '-.', '--', '-']
    
    # SUBPLOT 1
    i=0
    # loop over all eruptions and plot the seasonal anomalies on one graph
    for v in mmm_dataset.volcano:
        if v.data == 'Krakatoa':
            mmm_dataset.sel(volcano=v).plot(ax=ax, label = v.data, color = color_cycle[i], ls=line_style[i], lw='1.5') # plot anomalies 

        else:
            mmm_dataset.sel(volcano=v).plot(ax=ax, label = v.data, color = color_cycle[i], ls=line_style[i]) # plot anomalies 
        i = i+1
    
    if p10 is not None:
        ax.fill_between(p10.time.data, p10.data, p90.data, color='lightgrey')

    comp_dataset.plot(color = 'k', ax=ax, label = 'Composite') 
    
    # add red line for 0 line of eruptions 
    ax.axvline(x=0, color = 'r', linestyle = '--', alpha = 0.9, linewidth='1.5')
    
    # put legend on each suplot
    if subplot_legend is not None:
        ax.legend(loc="upper left")
    
    ax.set_facecolor('white')
    ax.grid(which='major', linestyle='-', linewidth='0.5', color='grey') # customise major grid
    ax.minorticks_on() # need this line in order to get the minor grid lines 
    ax.grid(which='minor', linestyle=':', linewidth='0.5', color='grey')
    
    # make horizontal zero line thicker 
    ax.axhline(y=0, color = 'k', alpha = 0.9, linewidth='0.8')
    
    # set's axis ticks to every 12 months 
    ax.set_xticks(np.arange(min(mmm_dataset.time), max(mmm_dataset.time)+1, 12))
     
    # set lables
    ax.set_xlabel(None) 
    ax.set_ylabel(None)
    
    return ax



# define a function for subplots of the nino3.4 index over time 
def nino34_plot(ds, e_date, thold, ax = None, **kwargs):
    """Create subplot of timeseries of SST anomalies for NINO34 index.  
    Values that exceed the specified threshold are shaded. 
    Shows the timing of the Krakatoa 1883 eruption (dotted vertical line).  
    Return the axis.  
    
    Args:
        ds (xarray): dataset of SST anomalies (use output from function "nino34")
        e_date (dict): dict of eruption dates
        thold (float): threshold value of NINO34 index
        ax (axis): axis
        **kwargs
    """
    import matplotlib.pyplot as plt, numpy as np, pandas as pd
    
    # checking if an axis has been defined and if not creates one with function "get current axes"
    if ax is None:
        ax = plt.gca()
        
    # SUBPLOT
    #plot data and fill if it's over the thresholds
    ds.plot(color='k', lw=1, ax=ax)
    ax.fill_between(ds.time.values, ds.values, thold, where=ds.values>thold, interpolate =True, color='crimson', alpha=0.6)
    ax.fill_between(ds.time.values, ds.values, -thold, where=ds.values<-thold, interpolate =True, color='royalblue', alpha=0.8)
    ax.axhline(0, color='k', lw=0.8)
    ax.axhline(thold, color='k', lw=0.8, linestyle = ':')
    ax.axhline(-thold, color='k', lw=0.8, linestyle = ':')

    # plot gridlines
    ax.grid(which='major', ls=':', lw='0.5', color='grey') # customise major grid
    ax.set_axisbelow(True) # sets the gridlines behind the data

    #set the frequency of the xticks 
    years = pd.date_range(ds.time.data[0], ds.time.data[-1], freq='YS')
    ax.set_xticks(years.values)
    ax.set_xticklabels(years.year) # .year shows only the year (not month)

    # remove xlabels and set ylabel
    ax.set_xlabel(None)
    ax.set_title(None)
    ax.set_ylabel(f'Sea surface temperature anomaly [$\degree$C]', size=12)

    # add dotted lineshowing the year of the krakatoa eruption
    ax.axvline(x=e_date[0], color = 'red', linestyle = '--', alpha = 0.9, linewidth='1.5')
    
    return ax


    
    
    
    
    
    
    
    
    
    
    
    

