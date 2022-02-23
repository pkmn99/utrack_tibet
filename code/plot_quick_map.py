import xarray as xr
import matplotlib.pyplot as plt
import cartopy.io.shapereader as shpreader
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.feature import ShapelyFeature

def load_data(region,time_scale='ymean'):
    d=xr.open_dataset('../data/processed/utrack_climatology_prec_0.5_mon_%s.nc'%region)
    return d

    # Use shaply interface https://stackoverflow.com/questions/20990381/how-to-add-custom-shapefile-to-map-using-cartopy
    # with statistical values

def plot_map(d, ax, plot_type='contourf',minmax=[]):
    # Load geographical data
    tb_shp=shpreader.Reader('../data/shp/DBATP_Polygon.shp')
#    china_shp=shpreader.Reader('../data/China_gis/行政边界/CN-sheng-A-WGS84.shp')
    china_shp=shpreader.Reader('../data/shp/China_provinces_with_around_countries.shp')
#    world_shp=shpreader.Reader('../../data/China_gis/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp')
#    tb=xr.open_dataset('../data/inputdata/DBATP_f19-touch.nc')

    tb_feature = ShapelyFeature(tb_shp.geometries(),
                                    ccrs.PlateCarree(), facecolor='none',edgecolor='g')
    china_feature = ShapelyFeature(china_shp.geometries(),
                                    ccrs.PlateCarree(), facecolor='none',edgecolor='k')
#    world_feature = ShapelyFeature(world_shp.geometries(),
#                                    ccrs.PlateCarree(), facecolor='none',edgecolor='grey', linewidth=0.5)
    ax.add_feature(tb_feature)
    ax.add_feature(china_feature)

    ax.set_extent([60, 140, 10, 55], ccrs.Geodetic())
    
    if plot_type=='contourf':
       d.plot.contourf(cmap='rainbow',
                       levels=[0,1,5,10,20,50,100,200,300,400,500,600],
                       vmax=800)
    if plot_type=='pcolormesh':
       d.plot.pcolormesh(ax=ax, cmap='rainbow')

 #   # Use hatch to represent significance
 #   cs = ax.contourf(p.lon,p.lat, 
 #                    (da1-da0).mean('time').where(dl.LANDFRAC_PFT.values>0.5).where(p<p_value),
 #                     colors='none',
 #                     hatches=['///'],
 #                     extend='lower',zorder=15
 #                     )

def make_plot(region,time_scale='ymean'):
    d=load_data(region,time_scale)

    # Plot change
    fig = plt.figure(figsize=[10,8])
    ax = fig.add_axes([0, 0, 1, 1], projection=ccrs.PlateCarree(),
                      frameon=False)

    ma=d['prec'].sum(dim='month')>1
    plot_map(d['prec'].sum(dim='month').where(ma), ax)
    ax.set_title(region)
    
    plt.savefig('../figure/fig_%s_%s.png'%(region,time_scale),dpi=300)
    print('figure saved succifully')

if __name__ == "__main__":
    make_plot('TP')
#    exp=['lindibaohu','caodibaohu','shuituliushi','shahuazhili']
#    for i in exp[2::]:
#        make_plot(i)
