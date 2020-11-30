import xarray as xr
import glob
import numpy as np
from hn2016_falwa.wrapper import qgpv_eqlat_lwa
from math import pi

def compute_LWA_ufunc(qgpv, area):
    clat = np.abs(np.cos(ylat*pi/180.))
    _, LWA = qgpv_eqlat_lwa(ylat, qgpv, area, Earth_radius*clat*dphi)
    return LWA

dims = ['z', 'lat', 'lon']
cesta = '/mnt/4data/CMAM/0A.daily/'
qgpv_list = sorted(glob.glob(cesta+'qgpv/logH/qgpv_logH*18.nc'))
Earth_radius = 6.378e+6

for i,qgpv_file in enumerate(qgpv_list[:]):
    print(qgpv_file)
    suffix = qgpv_file[len(cesta)+14:]
    da_qgpv = xr.open_dataarray(qgpv_file)
    
    if i==0:
        ylat = da_qgpv.lat.values
        nlon = da_qgpv.lon.shape[0]
        nlat = da_qgpv.lat.shape[0]
        dphi = (ylat[2]-ylat[1])*pi/180. # Equal spacing between latitude grid points, in radian
        area = 2.*pi*Earth_radius**2 *(np.cos(ylat[:,np.newaxis]*pi/180.)*dphi)/float(nlon) * np.ones((nlat,nlon))
        area = np.abs(area)
        area_xa = xr.DataArray(area, coords=[('lat', da_qgpv.lat), \
                               ('lon', da_qgpv.lon)])
        
    LWA_xa = xr.apply_ufunc(compute_LWA_ufunc, da_qgpv, area_xa, vectorize=True,\
                        input_core_dims = [dims[1:], dims[1:]], \
                        output_core_dims = [dims[1:]], \
                        dask = 'parallelized')
    
    LWA_xa.name = 'lwa'
    outfile = '{0}lwa/logH/lwa{1}'.format(cesta, suffix)
    print(outfile)
    LWA_xa.to_netcdf(outfile)

