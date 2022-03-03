import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import cartopy.crs as ccrs
import matplotlib as mpl

from process_et_data import subset_tb,get_tb_mask
from plot_prec_con_map import plot_map,uneven_cmap
from plot_et_prec_change import plot_map_tb

def make_plot():
    # load data
    dlc = xr.open_dataset('../data/processed/LC_fraction_clean.nc') # LC fration
    dlc_tb = subset_tb(dlc['LC_fraction'])
    tb_mask=get_tb_mask(scale='TP')
    dem = xr.open_dataset('../data/DEM_tb.nc') # DEM

    # combine LC group
    lc_list=[1,2,3,4,5]# forest
    eet_forest = dlc_tb.loc[lc_list].sum(dim='class')
    lc_list=[6,7,8,9]# shrub
    eet_shrub = dlc_tb.loc[lc_list].sum(dim='class')
    lc_list=[10]# grass
    eet_grass = dlc_tb.loc[lc_list].sum(dim='class')
    lc_list=[11,12,13,14]# other
    eet_other = dlc_tb.loc[lc_list].sum(dim='class')
    lc_list=[15,16]# bare and snow
    eet_baresnow = dlc_tb.loc[lc_list].sum(dim='class')

    levels=np.round(np.arange(0,1.01,0.2),1)
    levels2=np.arange(0,6501,500)

    fig = plt.figure(figsize=[12,12])
    gs = GridSpec(4, 4) # use gridspec but set ax1 mannualy
    ax2 = plt.subplot(gs[2,0],projection=ccrs.PlateCarree())
    ax3 = plt.subplot(gs[2,1],projection=ccrs.PlateCarree())
    ax4 = plt.subplot(gs[2,2],projection=ccrs.PlateCarree())
    ax5 = plt.subplot(gs[2,3],projection=ccrs.PlateCarree())
    
    ax6 = plt.subplot(gs[3,0],projection=ccrs.PlateCarree())
    ax7 = plt.subplot(gs[3,1],projection=ccrs.PlateCarree())
    ax8 = plt.subplot(gs[3,2],projection=ccrs.PlateCarree())
    ax9 = plt.subplot(gs[3,3],projection=ccrs.PlateCarree())
    
    # Panel: TP in China
    ax1 = fig.add_axes([ax2.get_position().x0, 0.45, 
                        ax5.get_position().x1 - ax2.get_position().x0, 0.5], 
                        projection=ccrs.PlateCarree(), frameon=True)
    
    plot_map(dem.DEM,ax1,levels2,cmap='terrain')
    
    # Add colorbar 
    cmap1,norm1=uneven_cmap(levels2,cmap='terrain')
    cbarbig1_pos = [ax1.get_position().x1 +.02, ax1.get_position().y0, 
                    0.02, ax1.get_position().height]
    caxbig1 = fig.add_axes(cbarbig1_pos)
    cbbig1 = mpl.colorbar.ColorbarBase(ax=caxbig1, cmap=cmap1, 
                                       norm=norm1, orientation='vertical',
                                       ticks=levels2)
    cbbig1.ax.set_yticklabels(levels2,fontsize=10)
    cbbig1.set_label('Elevation (m)')
    # # gs[3:5, 1:3], projection=ccrs.LambertConformal()
    
    # Panel for LC types 
    plot_map_tb(eet_forest.where(tb_mask), ax2,levels)
    plot_map_tb(eet_grass.where(tb_mask), ax3,levels)
    plot_map_tb(eet_shrub.where(tb_mask), ax4,levels)
    plot_map_tb(eet_baresnow.where(tb_mask), ax5,levels)
    
    # Panel for ecological projects
    plot_map_tb(np.nan, ax6,levels,region='lindibaohu')
    plot_map_tb(np.nan, ax7,levels,region='caodibaohu')
    plot_map_tb(np.nan, ax8,levels,region='shuituliushi')
    plot_map_tb(np.nan, ax9,levels,region='shahuazhili')
    
    ax2.set_title('Forest')
    ax3.set_title('Grass')
    ax4.set_title('Shrub')
    ax5.set_title('Bare and snow')
    
    ax6.set_title('Forest protection')
    ax7.set_title('Grass protection')
    ax8.set_title('Erosion control')
    ax9.set_title('Desertification control')
    
    # manual munipulation of axes position
    for a in [ax2,ax3,ax4,ax5]:
        a.set_position([a.get_position().x0,
                         a.get_position().y0-0.05,
                         a.get_position().width,
                         a.get_position().height])
    
    for a in [ax6,ax7,ax8,ax9]:
        a.set_position([a.get_position().x0,
                         a.get_position().y0,
                         a.get_position().width,
                         a.get_position().height])
        
    # Add colorbar panel LC
    cmap2,norm2=uneven_cmap(levels,cmap='bwr')
    cbar2_pos = [ax5.get_position().x1 +.02, ax5.get_position().y0, 
                 0.02, ax5.get_position().height]
    cax2 = fig.add_axes(cbar2_pos)
    cbbig2 = mpl.colorbar.ColorbarBase(ax=cax2, cmap=cmap2, 
                                       norm=norm2, orientation='vertical',
                                       ticks=levels)
    cbbig2.ax.set_yticklabels(levels,fontsize=10)
    cbbig2.set_label('Fraction')


    fig.text(0.5,0.41,"Ecosystem types",fontsize=14,ha='center')
    fig.text(0.5,0.26,"Ecological projects",fontsize=14,ha='center')
    
    plt.savefig('../figure/figure_study_area.png',bbox_inches='tight',dpi=300)
    print('figure saved')

if __name__=="__main__":
    make_plot()
