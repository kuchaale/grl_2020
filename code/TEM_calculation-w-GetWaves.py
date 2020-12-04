# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.3'
#       jupytext_version: 1.0.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# + {"toc": true, "cell_type": "markdown"}
# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"></ul></div>
# -

import xarray as xr
#from aostools.climate import *
from tem_calculation import ComputeEPfluxDiv, ComputeEPflux #* ComputeVstar, ComputeWstar, , ComputeEPfluxDiv, 
import numpy as np
from pathlib import Path
import sys
import datetime
#from tqdm import tqdm
from dask.diagnostics import ProgressBar
import dask 


what = sys.argv[1]
scheduler='single-threaded'


# +
def process_out(out_da, out_var, units = None, standard_name = None, long_name = None):
    #out_da = out_da.to_dataset(name = out_var)
    
    out_da['plev'].attrs['units'] = 'Pa'
    out_da['plev'].attrs['long_name'] = "pressure" 
    out_da['plev'].attrs['standard_name'] = "air_pressure"  
    
    if units is not None:
        out_da[out_var].attrs['units'] = units
    if standard_name is not None:
        out_da[out_var].attrs['standard_name'] = standard_name
    if long_name is not None:
        out_da[out_var].attrs['long_name'] = long_name
        
    now = datetime.datetime.utcnow()
    out_da.attrs['description'] = 'Created ' + now.strftime("%Y-%m-%dT%H:%M:%SZ")  + ' using aostools'        

    out_file = root_path / out_var / (out_var+infile.name.split(var)[-1])
    print(out_file)
    delayed_obj = out_da.to_netcdf(out_file, compute = False)
    
def fluxes_calc(ds_all, w):
    ep1, ep2 = ComputeEPflux(ds_all, 'ua', 'va', 'ta', 'wap', wave = w, do_ubar = True)   
    
    var_name = f'ep_phi-wn{w}'
    ep1 = process_out(ep1.to_dataset(name = var_name), var_name, units = 'm2 s-2', \
        standard_name = f"meridional EP-flux component", \
        long_name = f'Meridional EP-flux component')#, zarr_bool = True)

    var_name = f'ep_p-wn{w}'
    ep2 = process_out(ep2.to_dataset(name = var_name), var_name, units = 'hPa m s-2', \
        standard_name = f"vertical EP-flux component", \
        long_name = f'Vertical EP-flux component')
    
    return ep1, ep2
    
def epfd_calc(ds_all, w):
    div = ComputeEPfluxDiv(ds_all[f'ep_phi-wn{w}'], ds_all[f'ep_p-wn{w}'])

    var_name = f'acceldiv-calc_wn{w}'
    div = process_out(div.to_dataset(name = var_name), var_name, units = 'm/s/day', \
                        standard_name = f"ep flux divergence", \
                        long_name = f'EP Flux Divergence')
    #print(div)
    return div



# -

root_path = Path('/mnt/4data/CMAM/0A.daily/')
var_ls = ['ua', 'va', 'wap', 'ta']    
wave_range = [-1]#range(1,4) #21
for ifile in range(32):
    if what != 'epfd20': 
        ds_ls = []
        for var in var_ls:
            print(var,)
            infiles = root_path / var 
            infile = sorted(list(infiles.glob(f'{var}_6hrPlev_CMAM_CMAM30-SD_r1i1p1_*010100-*123118.nc')))[ifile]
            print(infile)
            ds = xr.open_dataset(infile)
            ds_ls.append(ds)

        ds_all = xr.merge(ds_ls)
        ds_all['wap'] = ds_all['wap']/100.
        levs = ds_all.plev/100. #transfer to hPa
        lats = ds_all.lat
    #print(ds_all)

    if what == 'vstar':
        p0=1e3
        v_bar, t_bar = ComputeVertEddy(ds_all['va'].values, ds_all['ta'].values, levs.values, p0, wave=-1)
        dp  = np.gradient(levs)[np.newaxis,:,np.newaxis]
        vstar = v_bar - np.gradient(t_bar, edge_order=2)[1]/dp
        out_da = xr.DataArray(vstar, coords=[('time', ds_all.time), \
                                                         ('plev', ds_all.plev), \
                                                         ('lat', ds_all.lat)])
        out_var = 'vstar-calc'
        process_out(out_da, out_var, units = 'm s-1', \
                    standard_name = "northward_transformed_eulerian_mean_air_velocity", \
                    long_name = 'Northward Transformed Eulerian Mean Air Velocity')

        
    elif what == 'epfd':
        _, _, div1, div2 = ComputeEPfluxDiv(lats.values, levs.values, ds_all['ua'].values, \
                        ds_all['va'].values, ds_all['ta'].values, ds_all['wap'].values/100., \
                        do_ubar = True, wave = -1)
        out_da = xr.DataArray(div1+div2, coords=[('time', ds_all.time), \
                                                         ('plev', ds_all.plev), \
                                                         ('lat', ds_all.lat)])
        out_var = 'acceldiv-calc'
        process_out(out_da, out_var, units = 'm s-2', \
                    standard_name = "ep_flux_divergence", \
                    long_name = 'EP Flux Divergence')
    
    elif what == 'fluxes':
        delayed_results = []
        for w in wave_range:
            delayed_obj   = delayed_obj = dask.delayed(fluxes_calc)(ds_all, w) # ep2 
            delayed_results.append(delayed_obj)                
        
        with ProgressBar():
            results = dask.compute(*delayed_results, scheduler = scheduler)

    
    elif what == 'epfd20':

        delayed_results = []
        for w in wave_range:
            ds_ls = []
            var_ls = [f'ep_phi-wn{w}', f'ep_p-wn{w}']
            for var in var_ls:
                print(var,)
                infiles = root_path / var 
                infile = sorted(list(infiles.glob(f'{var}_6hrPlev_CMAM_CMAM30-SD_r1i1p1_*010100-*123118.nc')))[ifile]
                print(infile)
                ds = xr.open_dataset(infile)
                ds_ls.append(ds)

            ds_all = xr.merge(ds_ls)
            #print(ds_all)
            delayed_obj = dask.delayed(epfd_calc)(ds_all, w)
            delayed_results.append(delayed_obj)                
        
        with ProgressBar():
            results = dask.compute(*delayed_results)#, scheduler = scheduler)
        #sys.exit()

    else:
        print(f'{what} nonapplicable')
        sys.exit()

