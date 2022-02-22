import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl
import cartopy.io.shapereader as shpreader
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.feature import ShapelyFeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from plot_prec_con_map import plot_map,set_lat_lon,uneven_cmap,plot_subplot_label
from process_et_data import load_et_region, get_tb_mask

# get province list for China
def get_china_list(region):
    northwest=['Shaanxi','Gansu','Qinghai','Ningxia','Xinjiang']
    north=['Beijing','Tianjin','Hebei','Shanxi','Neimeng']
    northeast=['Liaoning','Jilin','Heilongjiang']
    east=['Shanghai','Jiangsu','Zhejiang','Anhui','Fujian','Jiangxi','Shandong','Taiwan']
    central_south=['Henan','Hubei','Hunan','Guangdong','Guangxi','Hainan','Hongkong']# ,'Macao'
    southwest=['Chongqing','Sichuan','Guizhou','Yunnan','Xizang']

    if region=='northwest':
        re=northwest
    if region=='north':
        re=north
    if region=='northeast':
        re=northeast
    if region=='east':
        re=east
    if region=='central_south':
        re=central_south
    if region=='southwest':
        re=southwest
    if region=='china':
        re=northwest+north+northeast+east+central_south+southwest
    return re

# plot map for TP extent
def plot_map_tb(d, ax, levels, minmax=[]):
    # Load geographical data
    tb_shp=shpreader.Reader('../data/shp/DBATP_Polygon.shp')
    tb_feature = ShapelyFeature(tb_shp.geometries(),
                                    ccrs.PlateCarree(), facecolor='none',edgecolor='g',linewidth=0.75)
    ax.add_feature(tb_feature)
    ax.set_extent([72, 105, 25, 40], ccrs.Geodetic())
    d.plot.contourf(cmap='bwr',
                    levels=levels,
                    add_colorbar=False,ax=ax)

def make_plot():
    # ET trends
    det=xr.open_dataset('../data/processed/Etrend_2000-2020_GLEAM_v3.5a_TP_mon.nc')
    
    # prec contribution
    dp = xr.open_dataset('../data/processed/prec_change_by_et_change_2000-2020_TP.nc')
    
    # zonal prec contribution
    ds=pd.read_csv('../data/processed/prec_change_by_et_mon_TP_zonal.csv')
    china_list=get_china_list('china')
    # only china
    ds_etp=ds.set_index('name').reindex(china_list)['precYear'].sort_values(ascending=False)

    # ET yearly data
    dety = load_et_region(scale='year')

    # mask for TP
    tb_mask=get_tb_mask('TP')
    # regional mean ET and trend over TP
    et_year=dety.where(tb_mask).mean(dim=['lat','lon'])
    et_trend=det.E_trend.where(tb_mask).mean(dim=['lat','lon'])
    
    # create levels for map
    levels1=[-6,-4,-2,-1,-0.5, -0.25, 0.25, 0.5,1,2,4,6]
    mycmap1,mynorm1=uneven_cmap(levels1,cmap='bwr') # for panel a
    
    levels2=[-50,-40,-20,-10,-5,-2,-1,-0.25,0.25,1,2,5,10,20,40,50]
    mycmap2,mynorm2=uneven_cmap(levels2,cmap='bwr') # for panel b
    
    #################### Panel a: ET regional trend
    fig = plt.figure(figsize=[10,8.5])
    ax1 = fig.add_axes([0.035, 0.475, 0.4, 0.4], projection=ccrs.PlateCarree(),
                                                frameon=False)
    plot_map_tb(det.E_trend.sum(dim='month'), ax1, levels1)
    set_lat_lon(ax1, range(75,105,10), range(25,40,10), label=True, pad=0.05, fontsize=10)
    ax1.set_title('ET trends',fontsize=12)
    
    # Add colorbar to big plot
    cbarbig1_pos = [ax1.get_position().x0, ax1.get_position().y0-0.03, ax1.get_position().width, 0.02]
    caxbig1 = fig.add_axes(cbarbig1_pos)
    
    cbbig1 = mpl.colorbar.ColorbarBase(ax=caxbig1, cmap=mycmap1, norm=mynorm1, orientation='horizontal',
                                      ticks=levels1)
    cbbig1.ax.set_xticklabels(levels1,fontsize=10)
    cbbig1.set_label('ET trend (mm/year)')

    ################### Panel b: ET time series and monthly trend
    # ET yearly change over TP
    ax2 = fig.add_axes([0.525, 0.525, 0.4, 0.3], frameon=True)
    et_year.to_pandas().plot(ax=ax2)
    ax2.set_ylabel('ET (mm)')
    ax2.set_xlabel('Year')
    ax2.set_ylim([260,325])
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.set_title('ET changes in TP',fontsize=12)

    # inset bar chart for ET monthly trend 
    ax2b = fig.add_axes([0.77, 0.55, 0.15, 0.125], frameon=True)
    et_trend.to_pandas().plot.bar(ax=ax2b)
    ax2b.set_ylabel('ET trend (mm/yr)',fontsize=8)
    ax2b.set_xlabel('')
    ax2b.spines['top'].set_visible(False)
    ax2b.spines['right'].set_visible(False)
    ax2b.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'],
                         fontsize=8,rotation=0)
    ax2b.set_yticklabels([0,0.1],fontsize=8)

    ################### Panel c:  ET induced prec change
    
    ax3 = fig.add_axes([0.035, 0.05, 0.4, 0.4], projection=ccrs.PlateCarree(),
                                                frameon=False)
    # due to the mismatch between color interval of contouf and functiion uneven_colors
    # uneven_colors return -1 colors with respect to levels; conturnf draws all levels, with default settings
    # use levels-1 to plot and levels for colorbar
    plot_map(dp.prec.sum(dim='month'), ax3, levels2[0:-1],cmap='bwr')
    set_lat_lon(ax3, range(70,140,20), range(10,51,20), label=True, pad=0.05, fontsize=10)
    
    ax3.set_title('ET-induced precipitation changes',fontsize=12)
    
    # Add colorbar to big plot
    cbarbig3_pos = [ax3.get_position().x0, ax3.get_position().y0-0.03, ax3.get_position().width, 0.02]
    caxbig3 = fig.add_axes(cbarbig3_pos)
    
    cbbig3 = mpl.colorbar.ColorbarBase(ax=caxbig3, cmap=mycmap2, norm=mynorm2, orientation='horizontal',
                                      ticks=levels2)
    cbbig3.ax.set_xticklabels(levels2,fontsize=10)
    cbbig3.set_label('Precipitation change (mm)')
    
    # ################ Panel d: ET induced prec change by region
    ax4 = fig.add_axes([0.525, 0.125, 0.4, 0.3],frameon=False)
    
    ds_etp.plot(kind='bar',ax=ax4, legend=False)
    ax4.set_ylabel('Precipitation change (mm)')
    ax4.set_xlabel('')
    ax4.set_title('ET-induced precipitation changes by province',fontsize=12)
    
    # panel label
    plot_subplot_label(ax1, 'a', left_offset=-0.1, upper_offset=0.15)
    plot_subplot_label(ax2, 'b', left_offset=-0.15,upper_offset=0)
    plot_subplot_label(ax3, 'c', left_offset=-0.1,upper_offset=0.1)
    plot_subplot_label(ax4, 'd', left_offset=-0.15,upper_offset=0)
    
    plt.savefig('../figure/figure_et_prec_change_0221.png',dpi=300,bbox_inches='tight')
    print('figure saved')

if __name__=="__main__":
   make_plot() 
