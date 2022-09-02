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
    et_gleam = load_et_region()
   # et_modis = load_et_region(et_data="MODIS")
    et_pml = load_et_region(et_data="PML")
    et_era5 = load_et_region(et_data="ERA5")
    cn_mask=get_china_mask() # China mask for ERA5
    cmap='coolwarm_r'
    levels1=np.arange(0,901,100) 
    mycmap1,mynorm1=uneven_cmap(levels1,cmap=cmap)

    # define map projection
    pr=ccrs.PlateCarree()

    #################### Panel a: ET
    fig = plt.figure(figsize=[13, 6])
    ax1 = fig.add_subplot(3, 1, 1, projection=pr)

    plot_map_tb(et_gleam.sum(dim='month'),ax=ax1, levels=levels1,pr=pr,cmap=cmap)
    ax1.set_title('GLEAM ET',fontsize=12)

    ################### Panel b: ET
    ax2 = fig.add_subplot(3, 1, 2, projection=pr)
    plot_map_tb(et_pml.sum(dim='month'),ax=ax2, levels=levels1,pr=pr,cmap=cmap)
    ax2.set_title('PML ET',fontsize=12)
    ################### Panel c: ET
    ax3 = fig.add_subplot(3, 1, 3, projection=pr)
    plot_map_tb(et_era5.sum(dim='month'),ax=ax3, levels=levels1,pr=pr,cmap=cmap)
    ax3.set_title('ERA5 ET',fontsize=12)

    # Add colorbar to big plot
    cbarbig1_pos = [ax3.get_position().x0, ax3.get_position().y0-0.03, ax3.get_position().width, 0.02]
    caxbig1 = fig.add_axes(cbarbig1_pos)
    
    cbbig1 = mpl.colorbar.ColorbarBase(ax=caxbig1, cmap=mycmap1, norm=mynorm1, orientation='horizontal',
                                      ticks=levels1)
    cbbig1.ax.set_xticklabels(levels1,fontsize=10)
    cbbig1.set_label('(mm/year)')
    
    # panel label
    plot_subplot_label(ax1, 'a', left_offset=-0.1, upper_offset=0.1)
    plot_subplot_label(ax2, 'b', left_offset=-0.1,upper_offset=0.1)
    plot_subplot_label(ax3, 'c', left_offset=-0.1,upper_offset=0.1)
    
    plt.subplots_adjust(hspace=0.2)
    plt.savefig('../figure/figure_et_comparison0729.png',dpi=300,bbox_inches='tight')
    print('figure saved')

if __name__=="__main__":
   make_plot() 
