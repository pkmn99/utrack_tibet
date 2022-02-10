import xarray as xr
import pandas as pd
import numpy as np
from rasterstats import zonal_stats
from affine import Affine

from plot_prec_con_map import load_zonal_prec
from plot_et_prec_change import get_china_list
from process_et_data import get_tb_mask

"""
Make affine. Parameters for TP 
"""
def make_affine(res=0.5,lat_north=70.5,lon_west=50):
    # a = width of a pixel
    # b = row rotation (typically zero)
    # c = x-coordinate of the upper-left corner of the upper-left pixel
    # d = column rotation (typically zero)
    # e = height of a pixel (typically negative)
    # f = y-coordinate of the of the upper-left corner of the upper-left pixel

    # Affine(0.041666666667, 0.0, -125.0208333333335,
    #        0.0, -0.041666666667, 49.9375000000025)
    a = res        # change in x with x
    b = 0           # change in y with x
    c = lon_west # 50      # x offset
    d = 0           # change in y with x
    e = -res        # change in y with y
    f = lat_north # 70.5 # y offset

    af = Affine(a, b, c, d, e, f)
    
    return af

# Use Qgis to export shpfile csv and then load in pyhton with modifications
# only need code and name (province for China and country for rest of the world)
def get_region_list(rerun=False):
    if rerun:
        regionlist=pd.read_csv('../../data/China_gis/china_and_around/China_provinces_with_around_countries.csv')
        regionlist['code']=regionlist['省代码']
    #     regionlist['name']=regionlist['省']
        regionlist.loc[35::,'code'] = regionlist.loc[35::,'ADM0_CODE']
    #     regionlist.loc[35::,'name'] = regionlist.loc[35::,'ADM0_NAME']
        regionlist.loc[:,'name'] = regionlist.loc[:,'ADM0_NAME']
        regionlist[['code','name']].to_csv('../data/processed/region_list.csv',index=False)
        print('csv file saved')
    else:
        regionlist = pd.read_csv('../data/processed/region_list.csv')
    return regionlist[['code','name']]

"""
Function to calculate zonal statistis, use region list as index
"""
def calculate_zonal(shape_fn, data, affine, name):
    re=get_region_list()
    zs = zonal_stats(shape_fn, data, stats=['mean'], affine=affine)
    return re.join(pd.DataFrame(zs)).set_index('name')['mean'].rename(name)

# save zonal resutls to csv file
# precipitation contribution by ET
def save_zonal_prec_con(source_region='TP',lc_type='all'):
    shape_fn = '../data/shp/China_provinces_with_around_countries.shp'
    af = make_affine()
    re=get_region_list()
    if lc_type=='all':
        df = xr.open_dataset('../data/processed/utrack_climatology_prec_0.5_mon_%s.nc'%source_region)
    else:
        df = xr.open_dataset('../data/processed/utrack_climatology_prec_0.5_mon_%s_%s.nc'%(source_region,lc_type))
    prec_temp = [calculate_zonal(shape_fn, df.prec[i].values, af, 'prec'+str(i+1)) for i in range(12)]
    
    df_result= pd.DataFrame(prec_temp).transpose()
    
    # Calculate year sum prec
    df_result.loc[:,'precYear']=df_result.iloc[:,0:12].sum(axis=1)

    if lc_type=='all':
        df_result.to_csv('../data/processed/prec_con_mon_%s_zonal.csv'%source_region)
    else:
        df_result.to_csv('../data/processed/prec_con_mon_%s_%s_zonal.csv'%(source_region,lc_type))
    print('zonal results saved for source_region %s and lc_type %s'%(source_region,lc_type))

# save zonal resutls for prec change induced by et change
def save_zonal_prec_et(source_region='TP'):
    shape_fn = '../data/shp/China_provinces_with_around_countries.shp'
    
    df = xr.open_dataset('../data/processed/prec_change_by_et_change_2000-2020_%s.nc'%source_region)
    af = make_affine()
    re=get_region_list()
    prec_temp = [calculate_zonal(shape_fn, df.prec[i].values, af, 'prec'+str(i+1)) for i in range(12)]
    
    df_result= pd.DataFrame(prec_temp).transpose()
    
    # Calculate year sum prec
    df_result.loc[:,'precYear']=df_result.iloc[:,0:12].sum(axis=1)
    df_result.to_csv('../data/processed/prec_change_by_et_mon_%s_zonal.csv'%source_region)
    print('zonal results saved')

# zonal annual precipitation, used to calculate relative contribution 
def save_zonal_prec(source_region='TP'):
    shape_fn = '../data/shp/China_provinces_with_around_countries.shp'
    df = xr.open_dataset('../data/prec_CMFD_V0106_B-01_01mo_050deg_2008-2017_ymonmean_clean.nc')
    af = make_affine()
    re=get_region_list()
    prec_temp = [calculate_zonal(shape_fn, df.prec[i].values, af, 'prec'+str(i+1)) for i in range(12)]
    
    df_result= pd.DataFrame(prec_temp).transpose()
    
    # Calculate year sum prec
    df_result.loc[:,'precYear']=df_result.iloc[:,0:12].sum(axis=1)
    df_result.to_csv('../data/processed/prec_mon_%s_zonal.csv'%source_region)
    print('zonal results saved')

# create summerized table for TP, subregion, and different 
def save_table():
    china_list=get_china_list('china')
    
    ds_tp_abs = load_zonal_prec(type='absolute',time_scale='season',rank=1000)
    ds_tp_rel = load_zonal_prec(type='relative',time_scale='season',rank=1000)
    
    ds_lin = load_zonal_prec(rank=1000,source_region='lindibaohu')
    ds_cao = load_zonal_prec(rank=1000,source_region='caodibaohu')
    ds_shui = load_zonal_prec(rank=1000,source_region='shuituliushi')
    ds_sha = load_zonal_prec(rank=1000,source_region='shahuazhili')
    
    ds_forest = load_zonal_prec(rank=1000,lc_type='forest')
    ds_shrub = load_zonal_prec(rank=1000,lc_type='shrub')
    ds_grass = load_zonal_prec(rank=1000,lc_type='grass')
    ds_baresnow = load_zonal_prec(rank=1000,lc_type='baresnow')
    ds_other = load_zonal_prec(rank=1000,lc_type='other')
    
    table1=ds_tp_abs.reindex(china_list).rename(columns={'precYear':'Annual'}).join(
           ds_tp_rel.reindex(china_list).rename(columns={'precYear':'Annual_rel'}),rsuffix='_rel').sort_values(by='Annual',ascending=False)
    
    table2=ds_cao.reindex(china_list).rename(columns={'precYear':'caodizhili'}).join(
           ds_lin.reindex(china_list).rename(columns={'precYear':'lindizhili'})).join(
           ds_shui.reindex(china_list).rename(columns={'precYear':'shuituliushi'})).join(
           ds_sha.reindex(china_list).rename(columns={'precYear':'shahuazhili'}))
    
    table3=ds_forest.reindex(china_list).rename(columns={'precYear':'forest'}).join(
           ds_shrub.reindex(china_list).rename(columns={'precYear':'shrub'})).join(
           ds_grass.reindex(china_list).rename(columns={'precYear':'grass'})).join(
           ds_baresnow.reindex(china_list).rename(columns={'precYear':'baresnow'})).join(
           ds_other.reindex(china_list).rename(columns={'precYear':'other'}))


    # within-TP contribution
    dpc = xr.open_dataset('../data/processed/utrack_climatology_prec_0.5_mon_TP.nc')
    dp = xr.open_dataset('../data/prec_CMFD_V0106_B-01_01mo_050deg_2008-2017_ymonmean_clean.nc')
    tb=get_tb_mask(scale='TP')
    # TP prec seasonal prec contribution 
    temp_pc=[dpc.prec.where(tb).mean(dim=['lat','lon']).sel(month=[3,4,5]).sum().values,
            dpc.prec.where(tb).mean(dim=['lat','lon']).sel(month=[6,7,8]).sum().values,
            dpc.prec.where(tb).mean(dim=['lat','lon']).sel(month=[9,10,11]).sum().values,
            dpc.prec.where(tb).mean(dim=['lat','lon']).sel(month=[12,1,2]).sum().values,
            dpc.prec.where(tb).mean(dim=['lat','lon']).sum().values]
    
   # TP prec 
    temp_p=[dp.prec.where(tb).mean(dim=['lat','lon']).sel(month=[3,4,5]).sum().values,
            dp.prec.where(tb).mean(dim=['lat','lon']).sel(month=[6,7,8]).sum().values,
            dp.prec.where(tb).mean(dim=['lat','lon']).sel(month=[9,10,11]).sum().values,
            dp.prec.where(tb).mean(dim=['lat','lon']).sel(month=[12,1,2]).sum().values,
            dp.prec.where(tb).mean(dim=['lat','lon']).sum().values]
    
   # construct table 4 
    d = {"TP": np.append(temp_pc,np.array(temp_pc)/np.array(temp_p))}
    table4=pd.DataFrame(d,index=['MAM','JJA','SON','DJF','Annual','MAM_rel','JJA_rel','SON_rel','DJF_rel','Annual_rel']).transpose()
    
    table1.join(table2).join(table3).append(table4,sort=False).to_csv('../data/processed/summary_table.csv')
    print('summary table saved')

if __name__=="__main__":
#    save_zonal_prec_con()
#    save_zonal_prec()
#    save_zonal_prec_et()
# save zonal for different TP subregions
#    save_zonal_prec_con(source_region='lindibaohu')
#    save_zonal_prec_con(source_region='caodibaohu')
#    save_zonal_prec_con(source_region='shahuazhili')
#    save_zonal_prec_con(source_region='shuituliushi')

# save zonal for different TP land cover groups 
#    save_zonal_prec_con(lc_type='forest')
#    save_zonal_prec_con(lc_type='grass')
#    save_zonal_prec_con(lc_type='shrub')
#    save_zonal_prec_con(lc_type='baresnow')
#    save_zonal_prec_con(lc_type='other')
    save_table()
