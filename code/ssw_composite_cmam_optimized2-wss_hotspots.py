import xarray as xr
import sys
import pandas as pd
import glob
from scipy import stats
import xarray.ufuncs as xrf
from itertools import product
from cftime import DatetimeNoLeap
#from dask.distributed import Client

#client = Client(set_as_default = True)

def open_date_file(file_path):
    df = pd.read_csv(file_path, index_col=0, parse_dates=True)
    df['BeginDate'] = df.BeginDate.apply(lambda t: pd.to_datetime(t, format='%Y-%m-%d'))
    return df

def ttest_1samp(a, popmean, dim):
    """
    This is a two-sided test for the null hypothesis that the expected value
    (mean) of a sample of independent observations `a` is equal to the given
    population mean, `popmean`
    
    Inspired here: https://github.com/scipy/scipy/blob/v0.19.0/scipy/stats/stats.py#L3769-L3846
    
    Parameters
    ----------
    a : xarray
        sample observation
    popmean : float or array_like
        expected value in null hypothesis, if array_like than it must have the
        same shape as `a` excluding the axis dimension
    dim : string
        dimension along which to compute test
    
    Returns
    -------
    mean : xarray
        averaged sample along which dimension t-test was computed
    pvalue : xarray
        two-tailed p-value
    """
    n = a[dim].shape[0]
    df = n - 1
    a_mean = a.mean(dim)
    d = a_mean - popmean
    v = a.var(dim, ddof=1)
    denom = xrf.sqrt(v / float(n))

    t = d /denom
    prob = stats.distributions.t.sf(xrf.fabs(t), df) * 2
    prob_xa = xr.DataArray(prob, coords=a_mean.coords, name = a_mean.name)
    return a_mean, prob_xa
    

var = sys.argv[1]
try:
    w_clim = sys.argv[2]
    w_clim = '_'+w_clim
    print('climatology wo SSW')
except  IndexError:
    w_clim = ''
    print('climatology w SSW')


what = sys.argv[3]
what_ls = ['anomalies', 'absolute', 'percentages']
if what not in what_ls:
    raise ValueError('could not find {0} within [{1},{2},{3}]'.format(what, *what_ls))

if var[:3].lower() in ['lwa']:
    lev_sys_fo = 'log_coord/'
    lev_sys_fi = '_logH'
elif var.lower() == 'acceldivmrho':
    lev_sys_fo = ''
    lev_sys_fi = '_6hrPlev'
else:
    lev_sys_fo = ''
    lev_sys_fi = ''



if var in ['lwatend', 'TEM-res2', 'TEM-res', 'TEM-res3', 'TEM-res3-new','TEM-res2-new', 'sink']:
    zarr = True
    #hourly_index =  pd.date_range('1979-01-01-06', '2010-12-31-18', freq='6H')
else:
    zarr = False
    #hourly_index = pd.date_range('1979-01-01', '2010-12-31-18', freq='6H')


signif = False
DJF_bool = True
cesta = '/mnt/4data/CMAM/0A.daily/'
max_lag = 10
line_width = 5 
ch_lev = 70
timescale = int(sys.argv[4])
type_ls = ['himalayas', 'westamer', 'eastasia']

if zarr:
    infiles = '{}{}{}_197901-201012.zarr'.format(cesta,var,lev_sys_fi)
    print(infiles)
else:
    files_path = cesta+var+'/'+lev_sys_fo+var+'_*_CMAM_CMAM30-SD_r1i1p1_??????????-??????????.nc'
    print('opening '+files_path)
    infiles = sorted(glob.glob(files_path))
    print(len(infiles))

def open_infile(infiles):
    if isinstance(infiles, str) & ('zarr' in infiles):
        ds = xr.open_zarr(infiles)
    else:
        ds = xr.open_mfdataset(infiles, concat_dim = 'time', parallel = True, combine='nested')
    return ds

try:
    DJF_bool = sys.argv[5]
    print('{} only'.format(DJF_bool))
except  IndexError:
    DJF_bool = ''

with open_infile(infiles) as ds:
    #datetimeindex = ds.indexes['time'].to_datetimeindex()
    #ds['time'] = datetimeindex
    with xr.open_dataset( cesta+var+'/'+lev_sys_fo+var+'_climatology'+w_clim+'.nc') as clim:
        print('composites construction')
        if len(DJF_bool) > 0:
            w_clim += '_{}only'.format(DJF_bool)


        for di, ssw_type in enumerate(type_ls):
            print("".ljust(line_width) + ssw_type)
            df_dates = open_date_file('accelogw_{}_hotspot@{}hPa_{}dayts_indexes.csv'.format(ssw_type, ch_lev, timescale))
        #sys.exit()
            xa_ls = []
            pv_ls = []
            for il, lag in enumerate(range(-max_lag,max_lag+1)):
            #print("".ljust(line_width*3)+str(lag)+' lag')
                dates = df_dates.set_index('BeginDate')
                dates = dates.index +pd.Timedelta(str(lag)+' days')
            #filter lags withi 29th February
            #dates = dates[dates.apply(lambda x: not (x.day in [29] and x.month in [2]))]
                dates = dates[~((dates.month == 2) & (dates.day == 29))]
            #filter dates shited to year out of range
            #dates = dates[dates.apply(lambda x: not (x.year in [1978,2011]))]
                dates =  dates[(dates.year != 2011) & (dates.year != 1978)]
                if DJF_bool == 'DJF':
                    dates = dates[(dates.month == 12) | (dates.month == 1) | (dates.month == 2)]
                elif DJF_bool == 'JF':
                    dates = dates[(dates.month == 1) | (dates.month == 2)]
                elif DJF_bool == 'J':
                    dates = dates[(dates.month == 1)]
                elif DJF_bool == 'F':
                    dates = dates[(dates.month == 2)]
                elif DJF_bool == 'D':
                    dates = dates[(dates.month == 12)]

            #choose all values within particular day
                #hourly_index_temp = hourly_index[hourly_index.floor('D').isin(dates)]
                hourly_index_temp = [DatetimeNoLeap(*x, hour) for x, hour in product(zip(dates.year, dates.month, dates.day), range(0,24,6))]
                if var in ['lwatend']:
                    ds['time'] = xr.cftime_range('1979-01-01T06', '2010-12-31T18', freq='6H', calendar='noleap')
                #print(lag, dates.shape,)
                #print(hourly_index_temp)
                ds_sel = ds.sel(time = hourly_index_temp[:]) 
                print(ds_sel[var].shape)
                if what == 'percentages':
                    comp = (ds_sel[var].groupby('time.month') - clim[var]).groupby('time.month')/clim[var]*100.
                elif what == 'absolute':
                    comp = ds_sel[var]
                elif what == 'anomalies':
                    comp = ds_sel[var].groupby('time.month') - clim[var]

            
                if signif:
                    comp_m, pvalues = ttest_1samp(comp, 0.0, dim = 'time') # calculates stat. signaficance using t-test
                    pv_ls.append(pvalues)
                else:
                    comp_m = comp.mean('time')
                xa_ls.append(comp_m)

        
            print("".ljust(line_width*2) + 'concatenation')
            xa_comp = xr.concat(xa_ls, dim = 'lag')
            xa_ls = []
            xa_comp['lag'] = range(-max_lag, max_lag+1)
            outfile = "{}composites{}/{}{}_{}_comp_{}_{}days.nc".format(cesta, w_clim, var, lev_sys_fi, what, type_ls[di], timescale)
            print("".ljust(line_width*3) + outfile)
            xa_comp.to_netcdf(outfile)

            if signif:
                xa_pval = xr.concat(pv_ls, dim = 'lag')
                pv_ls = []
                xa_pval['lag'] = range(-max_lag, max_lag+1)
                print("".ljust(line_width*2) + 'saving')
                xa_pval.to_netcdf(cesta+'composites'+w_clim+'/'+var+lev_sys_fi+'_pvalues_comp_'+type_ls[di]+'.nc')

print('done')
