import xarray as xr
import sys
import random
from scipy import stats
import glob
from resampling import _resample_iterations_idx

random.seed(0)

def g_kde(y, x):
    """Firstly, kernel density estimation of the probability density function of randomized anomalies.
    Secondly, evaluates the estimated pdf on a set of points.
    
    Args:
        y (np.array): datapoints to estimate from (randomized anomalies)
        x (np.array): datapoints to be evaluated (composite values)
    Returns:
        np.array: the estimated pdf on composite values
    """
    kde = stats.gaussian_kde(y)
    return kde(x)

var = sys.argv[1] # ta, ua ; input variable
time_scale = sys.argv[2] # 20 or 30 ; input timescale
its = int(sys.argv[3]) # 10000 ; number of samples
what = sys.argv[4] # anomalies ; what kind of anomalies
DJF_bool = sys.argv[5] # DJF only (bool)
rechunk = True # allows rechunking in xr.apply_ufunc

comp_name_ls = ['himalayas', 'eastasia', 'westamer'] # list of hitspots (composites)
# numbers of events in composites
if DJF_bool == 'DJF':                                                          
    size_dict = {'20': [37,37,25], '30': []}     
    w_clim = f'_{DJF_bool}only'
else:
    size_dict = {'20': [45,74,36], '30': [38,66,35]}
    w_clim = ''

line_width = 5
lev_sys_fo = ''

# special multiplication factor
if var in ['sink','lwatend','dzmuadt']:
    factor = 24*3600
else:
    factor = 1


# open time-series for anomaly calculation
root_path = '/mnt/nas4.meop2/meop40.data.model/CMAM/0A.daily/'  
if var in ['sink', 'TEM-res3-new', 'TEM-res2-new','lwatend']:
    infiles = f'{root_path}{var}_197901-201012.zarr'
    ds = xr.open_zarr(infiles)
else:
    infiles = glob.glob(f'{root_path}{var}/{var}_6hr*_CMAM_CMAM30-SD_r1i1p1_????0101??-*.nc')
    print(len(infiles))
    ds = xr.open_mfdataset(infiles, concat_dim='time', parallel=True)#, engine = 'h5netcdf')

# open climatology (pre-calculated) for anomaly calculation
clim = xr.open_dataset(f'{root_path}{var}/{lev_sys_fo}/{var}_climatology_woSSW.nc')

# anomaly calculation
if what == 'percentages':                                                                                
    ts_sel_anom = (ds[var].groupby('time.month') - clim[var]).groupby('time.month')/clim[var]*100.
elif what == 'absolute':
    ts_sel_anom = ds[var]
elif what == 'anomalies':
    ts_sel_anom = (ds[var].groupby('time.month') - clim[var])*factor
print('input data opened')

# for-loop via composites
for comp_name,size in zip(comp_name_ls, size_dict[time_scale]):
    print(comp_name, size)
    # samples generation (loaded from external function)
    rnd_arr = _resample_iterations_idx(ts_sel_anom, 
                                    its, 'time', replace=True, chunk=False, dim_max = size) 
    print("".ljust(line_width)+'{} samples generated'.format(its))
    # load of composite dataarray
    comp_file = f'{root_path}composites_woSSW{w_clim}/{var}_{what}_comp_{comp_name}_{time_scale}days.nc'
    ds_comp = xr.open_dataarray(comp_file)*factor
    print("".ljust(line_width)+'{} opened'.format(comp_file))
    
    # statistical significance calculation (vectorized g_kde)
    da_kde = xr.apply_ufunc(g_kde, rnd_arr, ds_comp,\
                       input_core_dims=[['iteration'], []],\
                       vectorize=True, dask='parallelized',\
                       exclude_dims=set(("iteration",)), \
                       dask_gufunc_kwargs=dict(allow_rechunk=True), \
                       output_core_dims=[[]], \
                       output_dtypes=[ds_comp.dtype])
    print("".ljust(line_width)+'p-values calculated')
    # output the calculation                        
    outfile_name = f'{root_path}composites_woSSW{w_clim}/{var}_pvalues_from{its}_comp_{comp_name}_{time_scale}days.nc'
    da_kde.name = var
    da_kde.to_netcdf(outfile_name)
    print("".ljust(line_width)+'{} saved'.format(outfile_name))
    #del da_kde, rnd_arr, ds_comp
    print()


print('done')
