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
import glob
import numpy as np
import sys 

invar = sys.argv[1]
root_path = '/mnt/nas4.meop2/meop40.data.model/CMAM/0A.daily/' #'/mnt/4data/CMAM/0A.daily/'
infiles = sorted(glob.glob(f'{root_path}{invar}/{invar}_6hrPlev_CMAM_CMAM30-SD_r1i1p1_*-*18.nc'))

# +
var = f'dzm{invar}dt'
cesta_out = f'{root_path}{var}/'
for i, infile in enumerate(infiles):
      
    suffix = infile.split(invar)[-1]#infile_u.split('/lwa_')[1]
    outfile = f'{cesta_out}{var}{suffix}'
    
    da = xr.open_dataset(infile)[invar].mean('lon')
    
  
    da_out = da.differentiate('time', datetime_unit='s')
    da_out.name = var
    print(outfile)    
    da_out.to_netcdf(outfile)


