import xarray as xr
import numpy as np
import pandas as pd

"""
Save the cleaned multi-year monthly mean ET data as the output from gdalward break time index
"""
def save_gleam_et():
    et = xr.open_dataset('../data/E_2008-2017_GLEAM_v3.5b_ymonmean_360x720.nc')
    temp = np.dstack([et.Band1.values,
                      et.Band2.values,
                      et.Band3.values,
                      et.Band4.values,
                      et.Band5.values,
                      et.Band6.values,
                      et.Band7.values,
                      et.Band8.values,
                      et.Band9.values,
                      et.Band10.values,
                      et.Band11.values,
                      et.Band12.values,
                     ])
    
    myet=xr.DataArray(np.swapaxes(temp,2,0)[:,::-1,:],
            coords=[range(1,13), et.lon[::-1], et.lat],
                         dims=['month','lat','lon'],
                         name='E')
    
    myet.to_netcdf('../data/E_2008-2017_GLEAM_v3.5b_ymonmean_360x720_clean.nc')
    print('The cleaned GLEAM ET data saved ')

# Save monthly cleaned ET for 2000 to 2020 at 0.5 deg
def save_gleam_et_mon():
    d = xr.open_dataset('../../data/GLEAM/E_2000-2020_GLEAM_v3.5a_MO_360x720.nc')
    d_temp = np.zeros([252,360,720])
    for i in range(252):
        d_temp[i:,:,]=np.flipud(np.moveaxis(d['Band'+str(i+1)].values,0,1))
    
    myet=xr.DataArray(d_temp,
            coords=[pd.date_range('2000','2021',freq='M'),d.lon[::-1], d.lat],
                         dims=['time','lat','lon'],
                         name='E')
    myet.to_netcdf('../data/E_2000-2020_GLEAM_v3.5a_MO_360x720_clean.nc')
    print('The cleaned GLEAM ET data saved ')

# Save data to TP extent
def save_ChinaForcing_prec():
    et = xr.open_dataset('../../data/ChinaForcing/Data_forcing_01mo_010deg/prec_CMFD_V0106_B-01_01mo_050deg_2008-2017_ymonmean.nc')
    temp = np.dstack([et.Band1.values,
                      et.Band2.values,
                      et.Band3.values,
                      et.Band4.values,
                      et.Band5.values,
                      et.Band6.values,
                      et.Band7.values,
                      et.Band8.values,
                      et.Band9.values,
                      et.Band10.values,
                      et.Band11.values,
                      et.Band12.values,
                     ])

    myet=xr.DataArray(np.moveaxis(temp,2,0)[:,::-1,:],
            coords=[range(1,13), et.lat[::-1], et.lon],
                         dims=['month','lat','lon'],
                         name='prec')

    myet.to_netcdf('../data/prec_CMFD_V0106_B-01_01mo_050deg_2008-2017_ymonmean_clean.nc')
    print('The cleaned ChinaForcing prec data saved ')

# return a subset of global data for TP, for lon,lat starting at 0.25 deg
def subset_tb(df):
    ma_boundary=(df.lat>=0) & (df.lat<=70.25) & (df.lon>=50) & (df.lon<=180.25)
    return df.where(ma_boundary,drop=True)

# return the et fraction based on MODIS data for different land cover groups
# return 1 for all type (default)
def et_lc_fraction(lc_type='all'):
    if lc_type=='all':
        det_lc=xr.DataArray([1])
    else:
        if lc_type=='forest':
            lc_list=[1,2,3,4,5]# forest
        if lc_type=='shrub':
            lc_list=[6,7,8,9]# shrub
        if lc_type=='grass':
            lc_list=[10]# grass
        if lc_type=='baresnow':
            lc_list=[15,16]# bare and snow
        if lc_type=='other':
            lc_list=[11,12,13,14]# other

        det = xr.open_dataset('../data/processed/et_fraction_clean.nc')
        det_tb = subset_tb(det['ET_fraction_for_land_cover'])
        det_lc = det_tb.sel(landcover=lc_list).sum(dim='landcover')
    return det_lc.values

# load et data values for different land cover groups at different time scales
def load_et_region(source_region='TP',scale='year',lc_type='all'):
    if scale=='year':
        dfe = xr.open_dataset('../data/E_2008-2017_GLEAM_v3.5b_ymonmean_360x720_clean.nc')
    if scale=='month':
        dfe = xr.open_dataset('../data/E_2000-2020_GLEAM_v3.5a_MO_360x720_clean.nc')
    if source_region=='TP':
      #  ma_boundary=(dfe.lat>=0) & (dfe.lat<=70.25) & (dfe.lon>=50) & (dfe.lon<=180.25)    
        etf = et_lc_fraction(lc_type=lc_type) # et fraction for land cover group
        dfe=subset_tb(dfe) * etf # multiple total et by et fraction
    return dfe.E

# Trend estimation using polyfit for different months
# df is a dataset; requirs new xarray version, mync env
def get_linear_trend(df):
    df['time'] = np.arange(df.time.shape[0])
    return df.polyfit(dim='time',deg=1,skipna=True)

def save_et_trend_mon():
    det = load_et_region(source_region='TP',scale='month')
    d_temp = np.zeros([12,det.shape[1],det.shape[2]])
    for i in range(12):
        d_temp[i,:,:]=get_linear_trend(det.isel(time=det.time.dt.month==(i+1)))['polyfit_coefficients'][0].values
    
    mytrend=xr.DataArray(d_temp,
            coords=[range(1,13),det.lat, det.lon],
                         dims=['month','lat','lon'],
                         name='E_trend')

    mytrend.to_netcdf('../data/processed/Etrend_2000-2020_GLEAM_v3.5a_TP_mon.nc')
    print('monthly ET trend saved')


if __name__=="__main__":
    save_gleam_et()
#    save_gleam_et_mon()
    save_ChinaForcing_prec()
    save_et_trend_mon()
