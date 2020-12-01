
import xarray as xr
import glob
import numpy as np

infiles = sorted(glob.glob('/mnt/4data/CMAM/0A.daily/lwa/logH/lwa_logH_CMAM_CMAM30-SD_r1i1p1_*-*.nc'))

cesta_out = '/mnt/4data/CMAM/0A.daily/lwatend/logH/lwatend_'
for infile in infiles:
    print(infile)
    suffix = infile.split('/lwa_')[1]
    ds=xr.open_dataset(infile, autoclose = True)
    dt = ds.time.diff('time').values[:,np.newaxis,np.newaxis,np.newaxis]/1e9 # in seconds
    tendency = ds.lwa.diff('time')/dt.astype(float)
    ds_out = tendency.to_dataset()
    outfile = '{0}{1}'.format(cesta_out, suffix)
    print(outfile)
    ds_out.to_netcdf(outfile)
