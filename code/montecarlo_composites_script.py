# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"></ul></div>

# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"></ul></div>

import xarray as xr
import sys
import random
from scipy import stats
import glob

random.seed(0)

def g_kde(y, x):
    kde = stats.gaussian_kde(y)
    return kde(x)

var = sys.argv[1] # ta, ua
time_scale = sys.argv[2] # 20 or 30
its = int(sys.argv[3]) # 10000
what = sys.argv[4]
DJF_bool = sys.argv[5]

comp_name_ls = ['himalayas', 'eastasia', 'westamer']
if DJF_bool == 'DJF':                                                          
    size_dict = {'20': [37,37,25], '30': []}     
    w_clim = f'_{DJF_bool}only'
else:
    size_dict = {'20': [45,74,36], '30': [38,66,35]}
    w_clim = ''

line_width = 5
lev_sys_fo = ''

if var in ['sink','lwatend']:
    factor = 24*3600
else:
    factor = 1


root_path = '/mnt/4data/CMAM/0A.daily/'
if var in ['sink', 'TEM-res3-new', 'TEM-res2-new','lwatend']:
    infiles = f'{root_path}{var}_197901-201012.zarr'
    ds = xr.open_zarr(infiles)
else:
    infiles = glob.glob(f'{root_path}{var}/{var}_6hr*_CMAM_CMAM30-SD_r1i1p1_????0101??-*.nc')
    print(len(infiles))
    ds = xr.open_mfdataset(infiles, concat_dim='time', parallel=True)#, engine = 'h5netcdf')
clim = xr.open_dataset(f'{root_path}{var}/{lev_sys_fo}/{var}_climatology_woSSW.nc')

if what == 'percentages':                                                                                
    ts_sel_anom = (ds[var].groupby('time.month') - clim[var]).groupby('time.month')/clim[var]*100.
elif what == 'absolute':
    ts_sel_anom = ds[var]
elif what == 'anomalies':
    ts_sel_anom = (ds[var].groupby('time.month') - clim[var])*factor
print('input data opened')

for comp_name,size in zip(comp_name_ls, size_dict[time_scale]):
    print(comp_name, size)
    rnd_means = xr.concat([ts_sel_anom.isel(time = random.sample(range(ds.time.shape[0]), size)).mean('time') \
                       for n in range(its)], dim = 'its')
    print("".ljust(line_width)+'{} samples generated'.format(its))
    comp_file = f'{root_path}composites_woSSW{w_clim}/{var}_{what}_comp_{comp_name}_{time_scale}days.nc'
    ds_comp = xr.open_dataarray(comp_file)*factor
    print("".ljust(line_width)+'{} opened'.format(comp_file))
    da_kde = xr.apply_ufunc(g_kde, rnd_means, ds_comp,\
                       input_core_dims=[['its'], []],\
                       vectorize=True, dask='allowed')
    print("".ljust(line_width)+'p-values calculated')
                        
    outfile_name = f'{root_path}composites_woSSW{w_clim}/{var}_pvalues_from{its}_comp_{comp_name}_{time_scale}days.nc'
    da_kde.to_netcdf(outfile_name)
    print("".ljust(line_width)+'{} saved'.format(outfile_name))
    del da_kde, rnd_means, ds_comp
    print()


print('done')
