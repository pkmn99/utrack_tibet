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
def save_zonal_prec_con(source_region='TP',lc_type='all',et_data='GLEAM_v3.5a'):
    shape_fn = '../data/shp/China_provinces_with_around_countries.shp'
    af = make_affine()
    re=get_region_list()
    if lc_type=='all':
        df = xr.open_dataset('../data/processed/utrack_climatology_prec_0.5_mon_%s_%s.nc'%(source_region,et_data))
    else:
        df = xr.open_dataset('../data/processed/utrack_climatology_prec_0.5_mon_%s_%s_%s.nc'%(source_region,lc_type,et_data))
    prec_temp = [calculate_zonal(shape_fn, df.prec[i].values, af, 'prec'+str(i+1)) for i in range(12)]
    
    df_result= pd.DataFrame(prec_temp).transpose()
    
    # Calculate year sum prec
    df_result.loc[:,'precYear']=df_result.iloc[:,0:12].sum(axis=1)

    if lc_type=='all':
        df_result.to_csv('../data/processed/prec_con_mon_%s_%s_zonal.csv'%(source_region,et_data))
    else:
        df_result.to_csv('../data/processed/prec_con_mon_%s_%s_%s_zonal.csv'%(source_region,lc_type,et_data))
    print('zonal results saved for source_region %s and lc_type %s, et_data %s'%(source_region,lc_type,et_data))

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

# calculate in-TP prec contribution
def intp_con(type='both',source_region='TP',scale='seasonal',lc_type='all',et_data='GLEAM_v3.5a'):
    # within-TP contribution
    if lc_type=='all':
        dpc = xr.open_dataset('../data/processed/utrack_climatology_prec_0.5_mon_%s_%s.nc'%(source_region,et_data))
    else:
        dpc = xr.open_dataset('../data/processed/utrack_climatology_prec_0.5_mon_%s_%s_%s.nc'%(source_region,lc_type,et_data))
    dp = xr.open_dataset('../data/processed/prec_CMFD_V0106_B-01_01mo_050deg_2008-2017_ymonmean_clean.nc')
    tb=get_tb_mask(scale='TP')

    # TP prec seasonal prec contribution 
    if (scale=='seasonal')&(type=='both'):
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

    if (scale=='year')&(type=='absolute'):
        temp_pc=dpc.prec.where(tb).mean(dim=['lat','lon']).sum().values
        temp_p=temp_pc # dummy variable
    return temp_pc,temp_p

# create summerized table for TP, subregion, and different 
# the et_data option only works for TP all land cover  
def save_table(et_data='GLEAM_v3.5a'):
    china_list=get_china_list('china')
    
    ds_tp_abs = load_zonal_prec(type='absolute',time_scale='season',rank=1000,et_data=et_data)
    ds_tp_rel = load_zonal_prec(type='relative',time_scale='season',rank=1000,et_data=et_data)
    
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
    
    table2=ds_cao.reindex(china_list).rename(columns={'precYear':'caodibaohu'}).join(
           ds_lin.reindex(china_list).rename(columns={'precYear':'lindibaohu'})).join(
           ds_shui.reindex(china_list).rename(columns={'precYear':'shuituliushi'})).join(
           ds_sha.reindex(china_list).rename(columns={'precYear':'shahuazhili'}))
    
    table3=ds_forest.reindex(china_list).rename(columns={'precYear':'forest'}).join(
           ds_shrub.reindex(china_list).rename(columns={'precYear':'shrub'})).join(
           ds_grass.reindex(china_list).rename(columns={'precYear':'grass'})).join(
           ds_baresnow.reindex(china_list).rename(columns={'precYear':'baresnow'})).join(
           ds_other.reindex(china_list).rename(columns={'precYear':'other'}))


    # within-TP contribution
    [temp_pc,temp_p]=intp_con(et_data=et_data)
    
   # construct table 4 
    d = {"TP": np.append(temp_pc,np.array(temp_pc)/np.array(temp_p))}
    table4=pd.DataFrame(d,index=['MAM','JJA','SON','DJF','Annual','MAM_rel','JJA_rel','SON_rel','DJF_rel','Annual_rel']).transpose()
    
    # table 5: inTP contribution for subregion
    temp_lin=intp_con(type='absolute',source_region='lindibaohu',scale='year')[0]
    temp_cao=intp_con(type='absolute',source_region='caodibaohu',scale='year')[0]
    temp_shui=intp_con(type='absolute',source_region='shuituliushi',scale='year')[0]
    temp_sha=intp_con(type='absolute',source_region='shahuazhili',scale='year')[0]
    d = {"TP": np.array([temp_lin,temp_cao,temp_shui,temp_sha])}
    table5=pd.DataFrame(d,index=['lindibaohu','caodibaohu','shuituliushi','shahuazhili']).transpose()

    # table 6: inTP contribution for different lc groups 
    temp_forest=intp_con(type='absolute',lc_type='forest',scale='year')[0]
    temp_shrub=intp_con(type='absolute',lc_type='shrub',scale='year')[0]
    temp_grass=intp_con(type='absolute',lc_type='grass',scale='year')[0]
    temp_baresnow=intp_con(type='absolute',lc_type='baresnow',scale='year')[0]
    temp_other=intp_con(type='absolute',lc_type='other',scale='year')[0]
    d = {"TP": np.array([temp_forest,temp_shrub,temp_grass,temp_baresnow,temp_other])}
    table6=pd.DataFrame(d,index=['forest','shrub','grass','baresnow','other']).transpose()

   # table1.join(table2).join(table3).append(table4,sort=False).append(table5,sort=False) \
   #         .append(table6,sort=False).to_csv('../data/processed/summary_table.csv')
    table1.join(table2).join(table3).append(table4.join(table5).join(table6),sort=False) \
             .to_csv('../data/processed/summary_table_%s.csv'%et_data)
    print('summary table saved; et_data=%s'%et_data)


if __name__=="__main__":
#    save_zonal_prec_con()
#    save_zonal_prec()
#    save_zonal_prec_et()

# save zonal for different TP subregions
#    et_data='GLEAM_v3.5a'
#    save_zonal_prec_con(source_region='lindibaohu',et_data=et_data)
#    save_zonal_prec_con(source_region='caodibaohu',et_data=et_data)
#    save_zonal_prec_con(source_region='shahuazhili',et_data=et_data)
#    save_zonal_prec_con(source_region='shuituliushi',et_data=et_data)

# save for different ET data
#    save_zonal_prec_con(et_data='MODIS')
#    save_zonal_prec_con(et_data='ERA5')

# save zonal for different TP land cover groups 
#    save_zonal_prec_con(lc_type='forest',et_data=et_data)
#    save_zonal_prec_con(lc_type='grass',et_data=et_data)
#    save_zonal_prec_con(lc_type='shrub',et_data=et_data)
#    save_zonal_prec_con(lc_type='baresnow',et_data=et_data)
#    save_zonal_prec_con(lc_type='other',et_data=et_data)

    save_table(et_data='MODIS')
    save_table(et_data='ERA5')
