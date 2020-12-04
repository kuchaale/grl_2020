import xarray as xr
import sys
import pandas as pd
import glob
import numpy as np

def open_date_file(file_path):
    df = pd.read_csv(file_path, index_col=0, parse_dates=True)
    df['BeginDate'] = df.BeginDate.apply(lambda t: pd.to_datetime(t, format='%Y-%m-%d'))
    df['EndDate'] = df.BeginDate.apply(lambda t: pd.to_datetime(t, format='%Y-%m-%d'))
    return df
    
def extract_ssw_years(df):
    b_years = df.BeginDate.apply(lambda x: x.year)
    e_years = df.EndDate.apply(lambda x: x.year)
    return merge_years(b_years, e_years)

def merge_years(split_years, displ_years):
    all_years = np.unique(list(split_years)+list(displ_years))
    return all_years


def open_infile(infiles):
    if isinstance(infiles, str) & ('zarr' in infiles):
        ds = xr.open_zarr(infiles)
    else:
        ds = xr.open_mfdataset(infiles, concat_dim = 'time', combine='by_coords')
    return ds


var = sys.argv[1]
how = sys.argv[2]
lev_sys = sys.argv[3]
period=sys.argv[4]
if period == 'month':
    period_out = ''
else:
    period_out = '_'+period

cesta = '/mnt/10T.data4/data/CMAM/0A.daily/'


if var in ['lwatend', 'TEM-res2', 'TEM-res', 'TEM-res3', 'TEM-res3-new','TEM-res2-new', 'sink']:
    zarr = True
else:
    zarr = False



if zarr:
    files_path = '{}{}{}_197901-201012.zarr'.format(cesta,var,lev_sys)
    print(files_path)
else:
    if lev_sys == '':
        files_path = cesta+var+'/'+var+'_6hr*_CMAM_CMAM30-SD_r1i1p1_????????*-????????*.nc'#????????00-*.nc'
    elif lev_sys == 'logH':
        files_path = cesta+var+'/log_coord/'+var+'_*_CMAM_CMAM30-SD_r1i1p1_????????00-*.nc'
    files_path = sorted(glob.glob(files_path))
    print('opening ')#+files_path)
    print(len(files_path))
#sys.exit()

file_path = 'ssw_dates_split_enso.csv'
df_dates_split = open_date_file(file_path)

file_path = 'ssw_dates_displ.csv'
df_dates_displ = open_date_file(file_path)

my = merge_years(extract_ssw_years(df_dates_split), extract_ssw_years(df_dates_displ))

with open_infile(files_path) as ds:
    fy = list(filter(lambda x: x.year not in my, ds.time.to_index()))
    #print(fy)
    ds = ds.sel(time = fy)
    print(ds.chunks)
    if lev_sys == 'logH':
        lev_sys = 'log_coord'
    outfile = cesta+var+'/'+lev_sys+'/'+var
    if how == 'std':
        print('std calculation')
        climatology = ds.groupby('time.{}'.format(period)).std('time')
        print('saving')
        outfile += '_variability'+period_out+'_woSSW.nc'
        print(outfile)
        climatology.to_netcdf(outfile)
    elif how == 'mean':
        print('mean calculation')
        climatology = ds.groupby('time.{}'.format(period)).mean('time')
        print('saving')
        outfile += '_climatology'+period_out+'_woSSW.nc'
        print(outfile)
        climatology.to_netcdf(outfile)

print('done')
