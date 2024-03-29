import xarray as xr
import pandas as pd
import numpy as np
from rasterstats import zonal_stats
from affine import Affine
from plot_prec_con_map import load_zonal_prec
from plot_et_prec_change import get_china_list
from process_et_data import get_tb_mask,load_prec_region

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
def save_zonal_prec_con(source_region='TP',lc_type='all',et_data='GLEAM_v3.5a',var='E'):
    shape_fn = '../data/shp/China_provinces_with_around_countries.shp'
    af = make_affine()
    re=get_region_list()
    if lc_type=='all':
        df = xr.open_dataset('../data/processed/utrack_climatology_prec_0.5_mon_%s_%s_%s.nc'%(source_region,et_data,var))
    else:
        df = xr.open_dataset('../data/processed/utrack_climatology_prec_0.5_mon_%s_%s_%s_%s.nc'%(source_region,lc_type,et_data,var))
    prec_temp = [calculate_zonal(shape_fn, df.prec[i].values, af, 'prec'+str(i+1)) for i in range(12)]
    
    df_result= pd.DataFrame(prec_temp).transpose()
    
    # Calculate year sum prec
    df_result.loc[:,'precYear']=df_result.iloc[:,0:12].sum(axis=1)

    if lc_type=='all':
        df_result.to_csv('../data/processed/prec_con_mon_%s_%s_%s_zonal.csv'%(source_region,et_data,var))
    else:
        df_result.to_csv('../data/processed/prec_con_mon_%s_%s_%s_%s_zonal.csv'%(source_region,lc_type,et_data,var))
    print('zonal results saved for source_region %s and lc_type %s, et_data %s, var %s'%(source_region,lc_type,et_data,var))

# save zonal resutls for prec change induced by et change
def save_zonal_prec_et_change(source_region='TP',var='E'):
    shape_fn = '../data/shp/China_provinces_with_around_countries.shp'
    
    df = xr.open_dataset('../data/processed/prec_change_by_%s_change_2000-2020_%s.nc'%(var,source_region))
    af = make_affine()
    re=get_region_list()
    prec_temp = [calculate_zonal(shape_fn, df.prec[i].values, af, 'prec'+str(i+1)) for i in range(12)]
    
    df_result= pd.DataFrame(prec_temp).transpose()
    
    # Calculate year sum prec
    df_result.loc[:,'precYear']=df_result.iloc[:,0:12].sum(axis=1)
    df_result.to_csv('../data/processed/prec_change_by_%s_mon_zonal.csv'%(var))
    print('zonal results saved for var %s'%var)

# zonal annual precipitation, used to calculate relative contribution 
def save_zonal_prec(source_region='TP',prec_data='ERA5'):
    shape_fn = '../data/shp/China_provinces_with_around_countries.shp'
    if prec_data=='CMFD':
        df = xr.open_dataset('../data/prec_CMFD_V0106_B-01_01mo_050deg_2008-2017_ymonmean_clean.nc')
    if prec_data=='ERA5':
        df=load_prec_region(prec_data=prec_data)
    af = make_affine()
    re=get_region_list()
    prec_temp = [calculate_zonal(shape_fn, df.prec[i].values, af, 'prec'+str(i+1)) for i in range(12)]
    
    df_result= pd.DataFrame(prec_temp).transpose()
    
    # Calculate year sum prec
    df_result.loc[:,'precYear']=df_result.iloc[:,0:12].sum(axis=1)
    df_result.to_csv('../data/processed/prec_mon_%s_%s_zonal.csv'%(source_region,prec_data))
    print('zonal results saved')

# calculate in-TP prec contribution
def intp_con(type='both',source_region='TP',scale='seasonal',lc_type='all',et_data='GLEAM_v3.5a',prec_data='ERA5',var='E'):
    # within-TP contribution
    if lc_type=='all':
        dpc = xr.open_dataset('../data/processed/utrack_climatology_prec_0.5_mon_%s_%s_%s.nc'%(source_region,et_data,var))
    else:
        dpc = xr.open_dataset('../data/processed/utrack_climatology_prec_0.5_mon_%s_%s_%s_%s.nc'%(source_region,lc_type,et_data,var))
    if prec_data=='ERA5':
        dp = load_prec_region()
    else:
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
def save_table(et_data='GLEAM_v3.5a', prec_data='ERA5',var='E'):
    china_list=get_china_list('china')
    
    ds_tp_abs = load_zonal_prec(type='absolute',time_scale='season',rank=1000,
                                et_data=et_data,prec_data=prec_data,var=var)
    ds_tp_rel = load_zonal_prec(type='relative',time_scale='season',rank=1000,
                                et_data=et_data,prec_data=prec_data,var=var)
    
    ds_lin = load_zonal_prec(rank=1000,source_region='lindibaohu',var=var)
    ds_cao = load_zonal_prec(rank=1000,source_region='caodibaohu',var=var)
    ds_shui = load_zonal_prec(rank=1000,source_region='shuituliushi',var=var)
    ds_sha = load_zonal_prec(rank=1000,source_region='shahuazhili',var=var)
    
    ds_forest = load_zonal_prec(rank=1000,lc_type='forest',var=var)
    ds_shrub = load_zonal_prec(rank=1000,lc_type='shrub',var=var)
    ds_grass = load_zonal_prec(rank=1000,lc_type='grass',var=var)
    ds_baresnow = load_zonal_prec(rank=1000,lc_type='baresnow',var=var)
    ds_other = load_zonal_prec(rank=1000,lc_type='other',var=var)
    
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


    # within-TP contribution, last element of temp_p is TP annual P
    [temp_pc,temp_p]=intp_con(et_data=et_data,prec_data=prec_data,var=var)
    
   # construct table 4 
    d = {"TP": np.append(temp_pc,np.array(temp_pc)/np.array(temp_p))}
    table4=pd.DataFrame(d,index=['MAM','JJA','SON','DJF','Annual','MAM_rel','JJA_rel','SON_rel','DJF_rel','Annual_rel']).transpose()
    
    # table 5: inTP contribution for subregion
    temp_lin=intp_con(type='absolute',source_region='lindibaohu',scale='year',var=var)[0]
    temp_cao=intp_con(type='absolute',source_region='caodibaohu',scale='year',var=var)[0]
    temp_shui=intp_con(type='absolute',source_region='shuituliushi',scale='year',var=var)[0]
    temp_sha=intp_con(type='absolute',source_region='shahuazhili',scale='year',var=var)[0]
    d = {"TP": np.array([temp_lin,temp_cao,temp_shui,temp_sha])}
    table5=pd.DataFrame(d,index=['lindibaohu','caodibaohu','shuituliushi','shahuazhili']).transpose()

    # table 6: inTP contribution for different lc groups 
    temp_forest=intp_con(type='absolute',lc_type='forest',scale='year',var=var)[0]
    temp_shrub=intp_con(type='absolute',lc_type='shrub',scale='year',var=var)[0]
    temp_grass=intp_con(type='absolute',lc_type='grass',scale='year',var=var)[0]
    temp_baresnow=intp_con(type='absolute',lc_type='baresnow',scale='year',var=var)[0]
    temp_other=intp_con(type='absolute',lc_type='other',scale='year',var=var)[0]
    d = {"TP": np.array([temp_forest,temp_shrub,temp_grass,temp_baresnow,temp_other])}
    table6=pd.DataFrame(d,index=['forest','shrub','grass','baresnow','other']).transpose()

   # table1.join(table2).join(table3).append(table4.join(table5).join(table6),sort=False) \
   #          .to_csv('../data/processed/summary_table_%s_%s.csv'%(et_data,prec_data))

    # summerzied table
    df=table1.join(table2).join(table3).append(table4.join(table5).join(table6),sort=False)

    # add relative change for land cover and eco region
    dpz=pd.read_csv('../data/processed/prec_mon_%s_%s_zonal.csv'%('TP',prec_data),
               index_col=0)
    df=df.join(dpz['precYear']).copy() # add precp column

    # Calculate ralative change
    for i in df.columns[10:19]:
        df.loc[:,i+'_rel'] = df.loc[:,i]/dpz['precYear']*100
        df.loc['TP':,i+'_rel'] = df.loc['TP',i]/temp_p[-1]*100 # relative for TP region

    # do rounding
#    df.iloc[:,5:10]=(df.iloc[:,5:10]*100)
    df.loc[:,['MAM_rel','JJA_rel','SON_rel','DJF_rel','Annual_rel']] \
            =df.loc[:,['MAM_rel','JJA_rel','SON_rel','DJF_rel','Annual_rel']]*100

    df[df>10]=df[df>10].round(0)
    df[(df>1)&(df<10)]=df[(df>1)&(df<10)].round(1)
    df[df<1]=df[df<1].round(2)
    # convert string format: abs(rel)
    for i in ['MAM','JJA','SON','DJF','Annual']:
        df.loc[:,i] = df.loc[:,i].astype(str) + '(' +df.loc[:,i+'_rel'].astype(str) +')'
    # convert string format for land cover and ecoregion
    for i in df.columns[10:19]:
        df.loc[:,i] = df.loc[:,i].astype(str) + '(' +df.loc[:,i+'_rel'].astype(str) +'%)'

    df.to_csv('../data/processed/summary_table_%s_%s_%s_1012.csv'%(et_data,prec_data,var))
    print('summary table saved; et data is %s prec data is %s, var is %s'%(et_data,prec_data,var))

if __name__=="__main__":
    var='Et'
#    save_zonal_prec_con(et_data='GLEAM_v3.5a',var=var)
#    save_zonal_prec(prec_data='ERA5')
#    save_zonal_prec_con(et_data='PML')
#    save_zonal_prec_et_change(var=var)

# save zonal for different TP subregions
    et_data='GLEAM_v3.5a'
#    save_zonal_prec_con(source_region='lindibaohu',et_data=et_data,var=var)
#    save_zonal_prec_con(source_region='caodibaohu',et_data=et_data,var=var)
#    save_zonal_prec_con(source_region='shahuazhili',et_data=et_data,var=var)
#    save_zonal_prec_con(source_region='shuituliushi',et_data=et_data,var=var)

# save for different ET data
#    save_zonal_prec_con(et_data='MODIS')
#    save_zonal_prec_con(et_data='ERA5')

# save zonal for different TP land cover groups 
#    save_zonal_prec_con(lc_type='forest',et_data=et_data,var=var)
#    save_zonal_prec_con(lc_type='grass',et_data=et_data,var=var)
#    save_zonal_prec_con(lc_type='shrub',et_data=et_data,var=var)
#    save_zonal_prec_con(lc_type='baresnow',et_data=et_data,var=var)
#    save_zonal_prec_con(lc_type='other',et_data=et_data,var=var)

#    save_table(et_data='MODIS')
    save_table(et_data='GLEAM_v3.5a',prec_data='ERA5',var=var)
#    save_table(et_data='PML',prec_data='ERA5')
