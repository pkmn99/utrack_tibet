import xarray as xr
import numpy as np
import pandas as pd
from process_et_data import load_et_region

# convert int8 to valid data values
def conversion_value(x):
    return np.exp(x*-0.1)

# make a source mask and return valid grid lat/lon within the mask
def make_region_mask(region='TP',track_type='source',save_index=False):
    # Tibet mask
    if region=='TP':
        tb=xr.open_dataset('../data/DBATP_360x720-touch.nc') # Tibet as an example
    if region=='lindibaohu':
        tb=xr.open_dataset('../data/%s_360x720.nc'%region)
    if region=='caodibaohu':
        tb=xr.open_dataset('../data/%s_360x720.nc'%region)
    if region=='shahuazhili':
        tb=xr.open_dataset('../data/%s_360x720.nc'%region)
    if region=='shuituliushi':
        tb=xr.open_dataset('../data/%s_360x720.nc'%region)
        
    tb_mask=xr.DataArray(np.flipud(tb.Band1.values==1),
                        # coords=[d.coords[track_type+'lat'],d.coords[track_type+'lon']],
                         coords=[np.arange(90,-90,-0.5),np.arange(0,360,0.5)],
                         dims=[track_type+'lat',track_type+'lon'])

    # get lat and lon for every grid within the mask
    # tip from https://stackoverflow.com/questions/40592630/get-coordinates-of-non-nan-values-of-xarray-dataset#
    # ma3=d.moisture_flow[:,:,0,0].where(tb_mask==1) # use the 0,0 grid as target region to produce source mask
    ma_stacked = tb_mask.stack(grid=['sourcelat','sourcelon'])
    myindex=ma_stacked[ma_stacked] # retain only valid grids
    if save_index:
        myindex.grid.to_pandas().reset_index().iloc[:,0:2].to_csv('../data/processed/grid_index_source_region_%s.csv'%region,index=False)
        print('index saved to file')
    return myindex, tb_mask

"""
The global utrack data is too large. save the subset data for a defined source region extent
For example, a big extent to cover Tibet (source) and downstream areas 
Longitude range is 0 to 360
"""
def save_subset_utrack_data(mon,source_region='TP',lat_south=0,lat_north=70, lon_west=50, lon_east=180):
    d = xr.open_dataset('../data/utrack_climatology_0.5_%s.nc'%mon)

    myindex=make_region_mask(region=source_region,track_type='source')[0]
#    myindex.grid.to_pandas().reset_index().iloc[:,0:2].to_csv('../data/processed/grid_index_source_region_%s.csv'%region)

    # set export lat lon boundary
    ma_boundary=(d.targetlat>=0) & (d.targetlat<=70) & (d.targetlon>=50) & (d.targetlon<=180)
    # detemine export boundary size 
    [n_lat,n_lon]=ma_boundary.where(ma_boundary,drop=True).shape
    
    # create 3d array to save subset data within export boundary [grid, targetlat, targelon]
    frc_export=xr.DataArray(np.zeros([myindex.shape[0],n_lat,n_lon]),
                            coords=[range(myindex.shape[0]), 
                                    ma_boundary.where(ma_boundary,drop=True).targetlat,
                                    ma_boundary.where(ma_boundary,drop=True).targetlon],
                            dims=['grid','targetlat','targetlon'],
                            name='moisture_flow'
                            )
    for i in range(myindex.shape[0]):
#    for i in range(10):
#        temp=conversion_value(d.moisture_flow.sel(sourcelat=myindex.sourcelat.values[i],
#                        sourcelon=myindex.sourcelon.values[i]).where(ma_boundary,drop=True))

        temp=conversion_value(d.moisture_flow.sel(sourcelat=myindex.sourcelat.values[i],
                                                  sourcelon=myindex.sourcelon.values[i]))
        temp2=temp/temp.sum() # renormalize to sum as one
        frc_export[i]=temp2.where(ma_boundary,drop=True)
        print('grid %d has percent sum of %f'%(i,frc_export[i].sum()))

    frc_export.to_netcdf('../data/processed/utrack_climatology_0.5_%s_%s.nc'%(mon,source_region))
    print('file saved')

# Create region grid index
def region_grid_index(region):
# region='lindibaohu'
    ind=pd.read_csv('../data/processed/grid_index_source_region_%s.csv'%region)
    ind.loc[:,region]=1
#     return ind.set_index(['sourcelat','sourcelon'])
    return ind

# Create subregion index for all grids in TP, with each subregion denoted by 1 in their column
def make_subregion_index():
    ind = region_grid_index('TP').merge(region_grid_index('lindibaohu'),how='left') \
    .merge(region_grid_index('caodibaohu'),how='left').merge(region_grid_index('shahuazhili'),how='left') \
    .merge(region_grid_index('shuituliushi'),how='left')
    return ind

"""
Calculate ET contribution to total precipitation over the region
"""
def save_prec_contribution(source_region='TP'):
    subregion_index=make_subregion_index()
    dfe = load_et_region()
    n_lat, n_lon = dfe.shape[1::]
    pre = np.zeros([12, n_lat, n_lon])
    for i in range(12):
        dfm = xr.open_dataset('../data/processed/utrack_climatology_0.5_%02d_TP.nc'%(i+1))
        # extract ET values for all grids at month i, dim = (13xx,)
        temp_et = dfe.sel(month=i+1,
                          lat=xr.DataArray(subregion_index.loc[subregion_index[source_region].dropna().index,'sourcelat'].values+0.25),
                          lon=xr.DataArray(subregion_index.loc[subregion_index[source_region].dropna().index,'sourcelon'].values+0.25))
        # expand 1d ET to size of 3d dfm, dim= (13xx, 14x, 26x)
        # based on https://stackoverflow.com/questions/60044087/expand-and-copy-1d-numpy-array-to-3d
        # prec contri = moisture pct for selected grid (3d: grid, lat, lon) * grid ET (expanded to 3d, grid, lat, lon)
        temp_p=dfm.moisture_flow.sel(grid=subregion_index[source_region].dropna().index.values) * np.tile(temp_et.values[:, np.newaxis, np.newaxis], (1, n_lat, n_lon))
        # calculate by subregion
        pre[i,:,:] = temp_p.sum(axis=0)
    
    # create 3d array to save subset data within export boundary [grid, targetlat, targelon]
    pre_export=xr.DataArray(pre, coords=[range(1,13), dfe.lat,dfe.lon],
                            dims=['month','lat','lon'],
                            name='prec')
    pre_export.to_netcdf('../data/processed/utrack_climatology_prec_0.5_mon_%s.nc'%(source_region))
    print('prec contribution file for region %s saved'%source_region)

def save_prec_contribution_by_etrend(source_region='TP'):
    subregion_index=make_subregion_index()
    dfe = xr.open_dataset('../data/processed/Etrend_2000-2020_GLEAM_v3.5a_TP_mon.nc')['E_trend']*21 #2000 to 2020, 21 years
    n_lat, n_lon = dfe.shape[1::]
    pre = np.zeros([12, n_lat, n_lon])
    for i in range(12):
        dfm = xr.open_dataset('../data/processed/utrack_climatology_0.5_%02d_TP.nc'%(i+1))
        # extract ET values for all grids at month i, dim = (13xx,)
        temp_et = dfe.sel(month=i+1,
                          lat=xr.DataArray(subregion_index.loc[subregion_index[source_region].dropna().index,'sourcelat'].values+0.25),
                          lon=xr.DataArray(subregion_index.loc[subregion_index[source_region].dropna().index,'sourcelon'].values+0.25))
        # expand 1d ET to size of 3d dfm, dim= (13xx, 14x, 26x)
        # based on https://stackoverflow.com/questions/60044087/expand-and-copy-1d-numpy-array-to-3d
        # prec contri = moisture pct for selected grid (3d: grid, lat, lon) * grid ET (expanded to 3d, grid, lat, lon)
        temp_p=dfm.moisture_flow.sel(grid=subregion_index[source_region].dropna().index.values) * np.tile(temp_et.values[:, np.newaxis, np.newaxis], (1, n_lat, n_lon))
        # calculate by subregion
        pre[i,:,:] = temp_p.sum(axis=0)
    
    # create 3d array to save subset data within export boundary [grid, targetlat, targelon]
    pre_export=xr.DataArray(pre, coords=[range(1,13), dfe.lat,dfe.lon],
                            dims=['month','lat','lon'],
                            name='prec')
    pre_export.to_netcdf('../data/processed/prec_change_by_et_change_2000-2020_%s.nc'%(source_region))
    print('prec contribution by et change file for region %s saved'%source_region)

if __name__=="__main__":
#    for i in range(12):
#        save_subset_utrack_data('%02d'%(i+1),source_region='TP')

#    make_region_mask(region='TP',track_type='source',save_index=True)
#    make_region_mask(region='lindibaohu',track_type='source',save_index=True)
#    make_region_mask(region='caodibaohu',track_type='source',save_index=True)
#    make_region_mask(region='shahuazhili',track_type='source',save_index=True)
#    make_region_mask(region='shuituliushi',track_type='source',save_index=True)

    save_prec_contribution(source_region='TP')
#    save_prec_contribution(source_region='lindibaohu')
#    save_prec_contribution(source_region='caodibaohu')
#    save_prec_contribution(source_region='shahuazhili')
#    save_prec_contribution(source_region='shuituliushi')
    save_prec_contribution_by_etrend()
