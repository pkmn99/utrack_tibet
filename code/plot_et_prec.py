import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
from plot_prec_con_map import plot_map,set_lat_lon,uneven_cmap,plot_subplot_label,get_china_mask
from process_et_data import load_et_region
from plot_et_prec_change import plot_map_tb

# figure show ET and P over TP
def make_plot():
    # annual mean ET and P over TP
    et = load_et_region()
#    dp = xr.open_dataset('../data/processed/prec_CMFD_V0106_B-01_01mo_050deg_2008-2017_ymonmean_clean.nc')
    dp = xr.open_dataset('../data/processed/prec_2008-2017_ERA5_ymonmean_360x720_clean.nc')
    cn_mask=get_china_mask() # China mask for ERA5

    levels1=np.arange(0,901,100) 
    mycmap1,mynorm1=uneven_cmap(levels1,cmap='bwr') # for panel a

    # define map projection
    pr=ccrs.PlateCarree()

    #################### Panel a: ET
    fig = plt.figure(figsize=[8, 5])
    
    ax1 = fig.add_axes([0.035, 0.475, 0.4, 0.4], projection=pr)
#    ax1 = fig.add_subplot(2, 1, 1, projection=pr)

    plot_map_tb(et.sum(dim='month'),ax=ax1, levels=levels1,pr=pr,cmap='bwr')

    set_lat_lon(ax1, range(75,105,10), range(25,40,10), label=True, pad=0.05, fontsize=10)
    ax1.set_title('ET',fontsize=12)

    ################### Panel b: annual Prec China
    ax2 = fig.add_axes([0.035, 0.05, 0.4, 0.4], projection=pr)
    plot_map(dp.prec.where(cn_mask).sum(dim='month',skipna=False),ax=ax2, levels=np.arange(0,901,100),lw=1,cmap='bwr')

    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    ax2.set_title('Precipitation',fontsize=12)
    set_lat_lon(ax2, range(70,140,20), range(10,51,20), label=True, pad=0.05, fontsize=10)


    # Add colorbar to big plot
    cbarbig1_pos = [ax2.get_position().x0, ax2.get_position().y0-0.03, ax2.get_position().width, 0.02]
    caxbig1 = fig.add_axes(cbarbig1_pos)
    
    cbbig1 = mpl.colorbar.ColorbarBase(ax=caxbig1, cmap=mycmap1, norm=mynorm1, orientation='horizontal',
                                      ticks=levels1)
    cbbig1.ax.set_xticklabels(levels1,fontsize=10)
    cbbig1.set_label('(mm/year)')
    
    ax1.set_position([ax2.get_position().x0,
                      ax1.get_position().y0,
                      ax2.get_position().width,
                      ax2.get_position().height])
    # panel label
    plot_subplot_label(ax1, 'a', left_offset=-0.1, upper_offset=0.1)
    plot_subplot_label(ax2, 'b', left_offset=-0.1,upper_offset=0.1)
    
    plt.subplots_adjust(hspace=0.1)
    plt.savefig('../figure/figure_et_prec0424.png',dpi=300,bbox_inches='tight')
    print('figure saved')

if __name__=="__main__":
   make_plot() 
