import xarray as xr
import pandas as pd
from rasterstats import zonal_stats
from affine import Affine

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
def save_zonal_prec(source_region='TP'):
    shape_fn = '../data/shp/China_provinces_with_around_countries.shp'
    df = xr.open_dataset('../data/processed/utrack_climatology_prec_0.5_mon_%s.nc'%source_region)
    af = make_affine()
    re=get_region_list()
    prec_temp = [calculate_zonal(shape_fn, df.prec[i].values, af, 'prec'+str(i+1)) for i in range(12)]
    
    df_result= pd.DataFrame(prec_temp).transpose()
    
    # Calculate year sum prec
    df_result.loc[:,'precYear']=df_result.iloc[:,0:12].sum(axis=1)
    df_result.to_csv('../data/processed/prec_mon_%s_zonal.csv'%source_region)
    print('zonal results saved')

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


if __name__=="__main__":
  #  save_zonal_prec()
    save_zonal_prec_et()

