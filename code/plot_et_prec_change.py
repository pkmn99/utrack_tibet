import xarray as xr
import pandas as pd
import numpy as np
import seaborn as sns
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
    if region=='south':
        re=central_south
    if region=='southwest':
        re=southwest
    if region=='china':
        re=northwest+north+northeast+east+central_south+southwest
    return re

# plot map for TP extent
def plot_map_tb(d, ax, levels, minmax=[],region='TP',pr=ccrs.PlateCarree(),cmap='bwr'):
    # Load geographical data
    tb_shp=shpreader.Reader('../data/shp/DBATP_Polygon.shp')
    if region=='lindibaohu':
        tb_shp2=shpreader.Reader('../../data/Tibet/TP_ecoproject/%s.shp'%region)
    if region=='caodibaohu':
        tb_shp2=shpreader.Reader('../../data/Tibet/TP_ecoproject/%s1.shp'%region)
    if region=='shuituliushi':
        tb_shp2=shpreader.Reader('../../data/Tibet/TP_ecoproject/%s.shp'%region)
    if region=='shahuazhili':
        tb_shp2=shpreader.Reader('../../data/Tibet/TP_ecoproject/%s1.shp'%region)
    tb_feature = ShapelyFeature(tb_shp.geometries(), pr, facecolor='none')
    if region=='TP':
        ax.add_feature(tb_feature,edgecolor='k',linewidth=1)
        d.plot.contourf(cmap=cmap, levels=levels, add_colorbar=False,ax=ax)
    else:
        ax.add_feature(tb_feature,edgecolor='k',linewidth=1)
        tb_feature2 = ShapelyFeature(tb_shp2.geometries(), pr, facecolor='blue')
        ax.add_feature(tb_feature2,facecolor='blue',alpha=0.5,linewidth=0.75)
    ax.set_extent([72, 105, 25, 40], ccrs.Geodetic())
   # ax.set_extent([72, 105, 25, 40], crs=pr)

# load zonal changes in prec contribution due to ET
def load_zonal_prec_con_change():
    # zonal prec contribution
    ds=pd.read_csv('../data/processed/prec_change_by_et_mon_TP_zonal.csv')
    china_list=get_china_list('china')
    # only china
    ds_etp=ds.set_index('name').reindex(china_list)['precYear'].sort_values(ascending=False).to_frame()
    # attach region
    northwest=['Shaanxi','Gansu','Qinghai','Ningxia','Xinjiang']
    north=['Beijing','Tianjin','Hebei','Shanxi','Neimeng']
    northeast=['Liaoning','Jilin','Heilongjiang']
    east=['Shanghai','Jiangsu','Zhejiang','Anhui','Fujian','Jiangxi','Shandong','Taiwan']
    central_south=['Henan','Hubei','Hunan','Guangdong','Guangxi','Hainan','Hongkong']# ,'Macao'
    southwest=['Chongqing','Sichuan','Guizhou','Yunnan','Xizang']

    ds_etp.loc[northwest,'Region']='northwest'
    ds_etp.loc[north,'Region']='north'
    ds_etp.loc[northeast,'Region']='northeast'
    ds_etp.loc[east,'Region']='east'
    ds_etp.loc[central_south,'Region']='south'
    ds_etp.loc[southwest,'Region']='southwest'
    ds_etp.loc[ds_etp['Region'].isnull(),'Region']='international'
    return ds_etp

def make_plot():
    # ET trends
    det=xr.open_dataset('../data/processed/Etrend_2000-2020_GLEAM_v3.5a_TP_mon.nc')
    
    # prec contribution
    dp = xr.open_dataset('../data/processed/prec_change_by_et_change_2000-2020_TP.nc')
    
#    # zonal prec contribution
    ds_etp=load_zonal_prec_con_change()

    # ET yearly data
    dety = load_et_region(scale='year')

    # mask for TP
    tb_mask=get_tb_mask('TP')
    # regional mean ET and trend over TP
    et_year=dety.where(tb_mask).mean(dim=['lat','lon'])
    et_trend=det.E_trend.where(tb_mask).mean(dim=['lat','lon'])
    
    # create levels for map
    levels1=[-6,-4,-2,-1,-0.5, -0.25, 0.25, 0.5,1,2,4,6]
    mycmap1,mynorm1=uneven_cmap(levels1,cmap='bwr_r') # for panel a
    
    levels2=[-50,-40,-20,-10,-5,-2,-1,-0.25,0.25,1,2,5,10,20,40,50]
    mycmap2,mynorm2=uneven_cmap(levels2,cmap='bwr_r') # for panel b

    # order of Chian division 
    region_list=['northwest','southwest','north','south','east','northeast']
    
    #################### Panel a: ET regional trend
    fig = plt.figure(figsize=[10,8.5])
    ax1 = fig.add_axes([0.035, 0.475, 0.4, 0.4], projection=ccrs.PlateCarree(),
                                                frameon=False)
    plot_map_tb(det.E_trend.sum(dim='month'), ax1, levels1,cmap='bwr_r')
    set_lat_lon(ax1, range(75,105,10), range(25,40,10), label=True, pad=0.05, fontsize=10)
    ax1.set_title('ET trends in TP',fontsize=12)
    
    # Add colorbar to big plot
    cbarbig1_pos = [ax1.get_position().x0, ax1.get_position().y0-0.03, ax1.get_position().width, 0.02]
    caxbig1 = fig.add_axes(cbarbig1_pos)
    
    cbbig1 = mpl.colorbar.ColorbarBase(ax=caxbig1, cmap=mycmap1, norm=mynorm1, orientation='horizontal',
                                      ticks=levels1)
    cbbig1.ax.set_xticklabels(levels1,fontsize=10)
    cbbig1.set_label('ET trend (mm/year/year)')

    ################### Panel b: ET time series and monthly trend
    # ET yearly change over TP
    ax2 = fig.add_axes([0.525, 0.525, 0.4, 0.3], frameon=True)
    et_year.to_pandas().plot(ax=ax2)
    ax2.set_ylabel('ET (mm/year)')
    ax2.set_xlabel('Year')
    ax2.set_ylim([260,325])
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.set_title('ET trends in TP',fontsize=12)

    # inset bar chart for ET monthly trend 
    ax2b = fig.add_axes([0.77, 0.55, 0.15, 0.125], frameon=True)
    et_trend.to_pandas().plot.bar(ax=ax2b)
    ax2b.set_ylabel('ET trend (mm/mon/year)',fontsize=8)
    ax2b.set_xlabel('')
    ax2b.spines['top'].set_visible(False)
    ax2b.spines['right'].set_visible(False)
    ax2b.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'],
                         fontsize=8,rotation=0)
    ax2b.set_yticklabels([0,0.1],fontsize=8)

    ################### Panel c:  Map of ET induced prec change
    ax3 = fig.add_axes([0.035, 0.05, 0.4, 0.4], projection=ccrs.PlateCarree(),
                                                frameon=False)
    # due to the mismatch between color interval of contouf and functiion uneven_colors
    # uneven_colors return -1 colors with respect to levels; conturnf draws all levels, with default settings
    # use levels-1 to plot and levels for colorbar
    plot_map(dp.prec.sum(dim='month'), ax3, levels2[0:-1],lw=1,cmap='bwr_r')
    set_lat_lon(ax3, range(70,140,20), range(10,51,20), label=True, pad=0.05, fontsize=10)
    
    ax3.set_title('Changes in ET-contributed precipitation',fontsize=12)
    
    # Add colorbar to big plot
    cbarbig3_pos = [ax3.get_position().x0, ax3.get_position().y0-0.03, ax3.get_position().width, 0.02]
    caxbig3 = fig.add_axes(cbarbig3_pos)
    
    cbbig3 = mpl.colorbar.ColorbarBase(ax=caxbig3, cmap=mycmap2, norm=mynorm2, orientation='horizontal',
                                      ticks=levels2)
    cbbig3.ax.set_xticklabels(levels2,fontsize=10)
    cbbig3.set_label('Precipitation contribution changes (mm/year)')
    
    # ################ Panel d: ET induced prec change by region
    ax4 = fig.add_axes([0.525, 0.125, 0.4, 0.3],frameon=True)
    sns.barplot(x="name", y="precYear", hue="Region", hue_order=region_list,
                data=ds_etp.reset_index(), dodge=False, ax=ax4)
    ax4.set_xticklabels(ax4.get_xticklabels(),rotation=90)
    ax4.set_ylabel('Precipitation contribution changes (mm/year)')
    ax4.set_xlabel('')
    ax4.set_title('Changes in ET-contributed precipitation',fontsize=12)
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)

    # Inset to show China region division
    ax4inset = fig.add_axes([0.6, 0.25, 0.15, 0.125], projection=ccrs.PlateCarree(),
                                                       frame_on=False)
    ax4inset.outline_patch.set_visible(False) # Turn off borader
    # Load geographical data
    for i,r in enumerate(['西北','西南','华北','华南','East','东北']):
        china_shp=shpreader.Reader('/media/liyan/HDD/Project/data/China_gis/七大分区/%s.shp'%r)
        china_feature = ShapelyFeature(china_shp.geometries(), ccrs.PlateCarree(), facecolor=sns.color_palette()[i])
        ax4inset.add_feature(china_feature,edgecolor='w', linewidth=0.1)
    ax4inset.set_extent([70, 140, 10, 50],ccrs.Geodetic())

    # Add panel label
    plot_subplot_label(ax1, 'a', left_offset=-0.1, upper_offset=0.15)
    plot_subplot_label(ax2, 'b', left_offset=-0.15,upper_offset=0)
    plot_subplot_label(ax3, 'c', left_offset=-0.1,upper_offset=0.1)
    plot_subplot_label(ax4, 'd', left_offset=-0.15,upper_offset=0)
    
    plt.savefig('../figure/figure_et_prec_change_0728.png',dpi=300,bbox_inches='tight')
    print('figure saved')

if __name__=="__main__":
   make_plot() 
