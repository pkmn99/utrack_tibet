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
def load_zonal_prec_con_change(var='E'):
    # zonal prec contribution
    ds=pd.read_csv('../data/processed/prec_change_by_%s_mon_zonal.csv'%var)
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

def make_plot(var='E'):
    # ET trends
    dee=xr.open_dataset('../data/processed/%strend_2000-2020_GLEAM_v3.5a_TP_mon.nc'%'E') # ET
    det=xr.open_dataset('../data/processed/%strend_2000-2020_GLEAM_v3.5a_TP_mon.nc'%'Et') # T
    
    # prec contribution
    dpe = xr.open_dataset('../data/processed/prec_change_by_%s_change_2000-2020_TP.nc'%'E')
    dpt = xr.open_dataset('../data/processed/prec_change_by_%s_change_2000-2020_TP.nc'%'Et')
    
#    # zonal prec contribution
    ds_ep=load_zonal_prec_con_change(var='E')
    ds_tp=load_zonal_prec_con_change(var='Et')
    # merge to retain order of E contribution
    dsp=ds_ep.join(ds_tp['precYear'],rsuffix='1') # E is precYear, and T is precYear1

    # ET yearly data
    dey = load_et_region(scale='year',var='E')
    dty = load_et_region(scale='year',var='Et')

    # mask for TP
    tb_mask=get_tb_mask('TP')
    # regional mean ET time series and trend over TP
    e_year=dey.where(tb_mask).mean(dim=['lat','lon'])
    t_year=dty.where(tb_mask).mean(dim=['lat','lon'])
    e_trend=dee.E_trend.where(tb_mask).mean(dim=['lat','lon'])
    t_trend=det.E_trend.where(tb_mask).mean(dim=['lat','lon'])
    
    # create levels for map
    if var=='E':
       # levels1=[-6,-4,-2,-1,-0.5, -0.25, 0.25, 0.5,1,2,4,6]
        levels1=[-4,-3,-2,-1,-0.5, -0.1, 0.1, 0.5,1,2,3,4]
       # levels2=[-50,-40,-20,-10,-5,-2,-1,-0.25,0.25,1,2,5,10,20,40,50]
        levels2=[-40,-20,-10,-5,-2,-1,-0.1,0.1,1,2,5,10,20,40]
    if var=='Et':
        levels1=[-4,-3,-2,-1,-0.5, -0.1, 0.1, 0.5,1,2,3,4]
        levels2=[-40,-20,-10,-5,-2,-1,-0.1,0.1,1,2,5,10,20,40]

    mycmap1,mynorm1=uneven_cmap(levels1,cmap='bwr_r') # for panel a
    mycmap2,mynorm2=uneven_cmap(levels2,cmap='bwr_r') # for panel b

    # order of Chian division 
    region_list=['northwest','southwest','north','south','east','northeast']
    
    #################### Panel a: ET regional trend
    fig = plt.figure(figsize=[10,8.5])
    ax1 = fig.add_axes([0.035, 0.475, 0.4, 0.4], projection=ccrs.PlateCarree(),
                                                frameon=False)
    plot_map_tb(dee.E_trend.sum(dim='month').where(tb_mask), ax1, levels1,cmap='bwr_r')
    print('ET trend over TP is %f for 21 years'%dee.E_trend.sum(dim='month').where(tb_mask).mean().values*21)
    set_lat_lon(ax1, range(75,105,10), range(25,40,10), label=True, pad=0.05, fontsize=10)
#    ax1.set_title('ET trends in TP',fontsize=12)

    # Inset to show T trend
    ax1inset = fig.add_axes([ax1.get_position().x0, ax1.get_position().y0-0.0125, 0.12, 0.10], 
            projection=ccrs.PlateCarree(), frame_on=True)
    plot_map_tb(det.E_trend.sum(dim='month').where(tb_mask), ax1inset, levels1,cmap='bwr_r')
    print('T trend over TP is %f for 21 years'%det.E_trend.sum(dim='month').where(tb_mask).mean().values*21)
    ax1inset.text(0.075,0.15,'$T$', fontsize=14, transform=ax1inset.transAxes)
    ax1.text(0.05,0.4,'$ET$', fontsize=14, transform=ax1.transAxes)
    
    # Add colorbar to big plot
    cbarbig1_pos = [ax1.get_position().x0, ax1.get_position().y0-0.03, ax1.get_position().width, 0.02]
    caxbig1 = fig.add_axes(cbarbig1_pos)
    
    cbbig1 = mpl.colorbar.ColorbarBase(ax=caxbig1, cmap=mycmap1, norm=mynorm1, orientation='horizontal',
                                      ticks=levels1)
    cbbig1.ax.set_xticklabels(levels1,fontsize=10)
    cbbig1.set_label('ET trend (mm/yr/yr)')

    ################### Panel b: ET time series and monthly trend
    # ET yearly change over TP
    ax2 = fig.add_axes([0.525, 0.525, 0.4, 0.3], frameon=True)
    ax2x = ax2.twinx()
    # plot use numeric value as years 
    l1= ax2.plot(range(2000,2021),e_year.to_pandas())
    l2= ax2x.plot(range(2000,2021),t_year.to_pandas(),color='g')
    ax2.legend(l1+l2,['ET','T'],frameon=False)#,fontsize='xx-small')
    ax2.set_xlim([2000,2020])
    ax2.set_xticks(range(2000,2021,5))
    ax2x.tick_params(axis='y', colors='g')
    ax2x.set_ylabel('T (mm/yr)')

    ax2.set_ylabel('ET (mm/yr)')
    ax2.set_xlabel('Year')
    if var=="E":
        ax2.set_ylim([240,325])
        ax2x.set_ylim([145,175])
    if var=="Et":
        ax2.set_ylim([135,175])

#    ax2.set_xlim([2000,2020])

    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
#    ax2.set_title('ET trends in TP',fontsize=12)
    ax2x.spines['top'].set_visible(False)

    # inset bar chart for ET monthly trend 
    ax2b = fig.add_axes([0.8, 0.55, 0.125, 0.125], frameon=True)
    e_trend.to_pandas().plot.bar(ax=ax2b)
    ax2b.scatter(range(12),t_trend,marker='_',color='k',zorder=10, s=30,lw=1)
    l = ax2b.legend(['T','ET'])
    ax2b.legend(reversed(l.legendHandles), ['ET','T'],fontsize='xx-small',frameon=False)

#    ax2b.patch.set_facecolor('orange')
    ax2b.patch.set_alpha(0.0)

#    ax2b.set_ylabel('ET trend (mm/mon/yr)',fontsize=8)
    ax2b.set_xlabel('')
    ax2b.spines['top'].set_visible(False)
    ax2b.spines['right'].set_visible(False)
    ax2b.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'],
                         fontsize=8,rotation=0)
    ax2b.set_yticks([0,0.1])
    ax2b.set_yticklabels([0,0.1],fontsize=8)

    ################### Panel c:  Map of ET induced prec change
    ax3 = fig.add_axes([0.035, 0.05, 0.4, 0.4], projection=ccrs.PlateCarree(),
                                                frameon=False)
    # due to the mismatch between color interval of contouf and functiion uneven_colors
    # uneven_colors return -1 colors with respect to levels; conturnf draws all levels, with default settings
    # use levels-1 to plot and levels for colorbar
    plot_map(dpe.prec.sum(dim='month'), ax3, levels2[0:-1],lw=1,cmap='bwr_r')
    set_lat_lon(ax3, range(70,140,20), range(10,51,20), label=True, pad=0.05, fontsize=10)
    
#    ax3.set_title('Changes in ET-contributed precipitation',fontsize=12)

    # Inset to show T induced prec contribution change
    ax3inset = fig.add_axes([ax3.get_position().x0, ax3.get_position().y0, 0.15, 0.11],
                   projection=ccrs.PlateCarree(), frame_on=False)
    plot_map(dpt.prec.sum(dim='month'), ax3inset, levels2[0:-1],lw=1,cmap='bwr_r', extent=[70, 140, 15, 50])
    ax3inset.text(0.8,0.2,'$P_T$', fontsize=14, transform=ax3inset.transAxes)
    ax3.text(0.9,0.1,'$P_{ET}$', fontsize=14, transform=ax3.transAxes)
    
    # Add colorbar to big plot
    cbarbig3_pos = [ax3.get_position().x0, ax3.get_position().y0-0.03, ax3.get_position().width, 0.02]
    caxbig3 = fig.add_axes(cbarbig3_pos)
    
    cbbig3 = mpl.colorbar.ColorbarBase(ax=caxbig3, cmap=mycmap2, norm=mynorm2, orientation='horizontal',
                                      ticks=levels2)
    cbbig3.ax.set_xticklabels(levels2,fontsize=10)
    cbbig3.set_label('Precipitation contribution changes (mm/yr)')
    
    # ################ Panel d: ET induced prec change by region
    ax4 = fig.add_axes([0.525, 0.125, 0.4, 0.3],frameon=True)
    sns.barplot(x="name", y="precYear", hue="Region", hue_order=region_list,
                data=dsp.reset_index(), dodge=False, ax=ax4)
    sns.scatterplot(data=dsp.reset_index(), x="name", y="precYear1",
                    zorder=10,marker="_",color='k',linewidth=2, s=50,
                    ax=ax4)

    ax4.set_xticklabels(ax4.get_xticklabels(),rotation=90)
    ax4.set_ylabel('Precipitation contribution changes (mm/year)')
    ax4.set_xlabel('')

    ax4.set_xlim([-0.5,29.5])
   # ax4.set_ylim([0,18])
#    ax4.set_yticks(np.arange(0,18.1,3))
#    ax4.set_yticklabels(['%d'%i for i in np.arange(0,18.1,3)]) # remove digit

#    ax4.set_title('Changes in ET-contributed precipitation',fontsize=12)
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)

    ## a very clumsy way to recreate a combined legend  because dont know how to get line handel
    h1, l1 = ax4.get_legend_handles_labels()
    legend1=ax4.legend(h1,l1,title='PET')
    legend2 = ax4.legend('T',loc='lower right',title='PT')
    ax4.legend(legend1.get_patches() + legend2.get_lines(),l1 + ['$P_T$'], title='$P_{ET}$')

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
    
    plt.savefig('../figure/figure_et_prec_change_1003.png',dpi=300,bbox_inches='tight')
    print('figure saved')

if __name__=="__main__":
   make_plot(var='E') 
