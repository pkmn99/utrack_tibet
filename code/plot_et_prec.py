import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
from plot_prec_con_map import plot_map,set_lat_lon,uneven_cmap,plot_subplot_label,get_china_mask
from process_et_data import load_et_region, get_tb_mask
from plot_et_prec_change import plot_map_tb

# figure show ET and P over TP
def make_plot():
    # annual mean ET and P over TP
    et = load_et_region()
    t = load_et_region(var='Et')
#    dp = xr.open_dataset('../data/processed/prec_CMFD_V0106_B-01_01mo_050deg_2008-2017_ymonmean_clean.nc')
    dp = xr.open_dataset('../data/processed/prec_2008-2017_ERA5_ymonmean_360x720_clean.nc')
    cn_mask=get_china_mask() # China mask for ERA5
    tb_mask=get_tb_mask(scale='TP')

    levels1=np.arange(0,901,100) 
    mycmap1,mynorm1=uneven_cmap(levels1,cmap='coolwarm_r') # for panel a

    # define map projection
    pr=ccrs.PlateCarree()

    #################### Panel a: ET
    fig = plt.figure(figsize=[8, 5])
    ax1 = fig.add_axes([0.135, 0.575, 0.3, 0.3], projection=pr)
    plot_map_tb(et.sum(dim='month').where(tb_mask),ax=ax1, levels=levels1,pr=pr,cmap='coolwarm_r')
    set_lat_lon(ax1, range(75,105,10), range(25,40,10), label=True, pad=0.05, fontsize=10)
    ax1.set_title('Evapotranspiration',fontsize=12)
    ax1.text(0.05,0.075,'ET: %dmm/yr'%et.sum(dim='month').where(tb_mask).mean().values, transform=ax1.transAxes,
            fontsize=10)

   ######################### Panel b: T
    ax2 = fig.add_axes([0.5, 0.575, 0.3, 0.3], projection=pr)
    plot_map_tb(t.sum(dim='month').where(tb_mask),ax=ax2, levels=levels1,pr=pr,cmap='coolwarm_r')
    set_lat_lon(ax2, range(75,105,10), range(25,40,10), label=True, pad=0.05, fontsize=10)
    ax2.set_title('Transpiration',fontsize=12)
    ax2.text(0.05,0.075,'T: %dmm/yr'%t.sum(dim='month').where(tb_mask).mean().values, transform=ax2.transAxes,fontsize=10)

    ################### Panel c: annual Prec China
    ax3 = fig.add_axes([0.265, 0.15, 0.4, 0.35], projection=pr)
    plot_map(dp.prec.where(cn_mask).sum(dim='month',skipna=False),ax=ax3, levels=np.arange(0,901,100),lw=1,cmap='coolwarm_r')

    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)

    ax3.set_title('Precipitation',fontsize=12)
    set_lat_lon(ax3, range(70,140,20), range(10,51,20), label=True, pad=0.05, fontsize=10)


    # Add colorbar to big plot
    cbarbig1_pos = [ax3.get_position().x0, ax3.get_position().y0-0.03, ax3.get_position().width, 0.02]
    caxbig1 = fig.add_axes(cbarbig1_pos)
    
    cbbig1 = mpl.colorbar.ColorbarBase(ax=caxbig1, cmap=mycmap1, norm=mynorm1, orientation='horizontal',
                                      ticks=levels1)
    cbbig1.ax.set_xticklabels(levels1,fontsize=10)
    cbbig1.set_label('(mm/yr)')
    
    # panel label
    plot_subplot_label(ax1, '(a)', left_offset=-0.11, upper_offset=0.1, fontsize=12)
    plot_subplot_label(ax2, '(b)', left_offset=-0.11,upper_offset=0.1, fontsize=12)
    plot_subplot_label(ax3, '(c)', left_offset=-0.14,upper_offset=0.1, fontsize=12)
    
    plt.subplots_adjust(hspace=0.1)
    plt.savefig('../figure/figure_et_prec1020.png',dpi=600,bbox_inches='tight')
    plt.savefig('../figure/figure_et_prec1020.pdf',bbox_inches='tight')
    print('figure saved')

if __name__=="__main__":
   make_plot() 
