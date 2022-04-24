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
from process_et_data import load_prec_region

# # Figure1  abs and relative contribution

# Add lat lon to map figure
def set_lat_lon(ax, xtickrange, ytickrange, label=False,pad=0.05, fontsize=8,pr=ccrs.PlateCarree()):
    lon_formatter = LongitudeFormatter(zero_direction_label=True, degree_symbol='')
    lat_formatter = LatitudeFormatter(degree_symbol='')
    ax.set_yticks(ytickrange, crs=pr)
    ax.set_xticks(xtickrange, crs=pr)
    if label:
        ax.set_xticklabels(xtickrange,fontsize=fontsize)
        ax.set_yticklabels(ytickrange,fontsize=fontsize)
        ax.tick_params(axis='x', which='both', direction='out', bottom=False, top=True,labeltop=True,labelbottom=False, pad=pad)
        ax.tick_params(axis='y', which='both', direction='out', pad=pad)

    else:
        ax.tick_params(axis='x', which='both', direction='out', bottom=True, top=False, labeltop=False, labelleft=False, labelbottom=False)
        ax.tick_params(axis='y', which='both', direction='out', left=True, labelleft=False)

    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.set_ylabel('')
    ax.set_xlabel('')

def plot_map(d, ax, levels, minmax=[],cmap='rainbow',extent=[70, 140, 10, 50],pr=ccrs.PlateCarree(),
             lw=2):
    # Load geographical data
    tb_shp=shpreader.Reader('../data/shp/DBATP_Polygon.shp')
    china_shp=shpreader.Reader('../data/shp/China_provinces_with_around_countries.shp')

    tb_feature = ShapelyFeature(tb_shp.geometries(), pr, facecolor='none')
    china_feature = ShapelyFeature(china_shp.geometries(), pr, facecolor='none')

    ax.add_feature(china_feature,edgecolor='dimgrey', linewidth=0.5)
    ax.add_feature(tb_feature,edgecolor='k',linewidth=lw)
    ax.set_extent(extent,ccrs.Geodetic())

    d.plot.contourf(cmap=cmap, levels=levels, add_colorbar=False,ax=ax)

# create uneven color levels for maps
# https://stackoverflow.com/questions/61897393/unevenly-irregularly-spaced-data-for-colorbar-with-evenly-spaced-colors
def uneven_cmap(levels,cmap='rainbow'):
    cmap_rb = plt.get_cmap(cmap)
    colors = cmap_rb(np.linspace(0, 1, len(levels) - 1))
    mycmap, mynorm = mcolors.from_levels_and_colors(levels, colors)
    return mycmap,mynorm


# create four colors for plotting bar chart
def bin_color(df,bins):
    # Calculate sample number within each bin
    count = np.histogram(df, bins=bins)[0]
    # Creat color for each bin
    mycolor='r'*count[-1] + 'm'*count[-2] + 'y' *count[-3] + 'b'*count[-4]
    return mycolor

def plot_subplot_label(ax, txt, left_offset=-0.05, upper_offset=0.05):
    ax.text(left_offset, 1+upper_offset, txt, fontsize=14, transform=ax.transAxes, fontweight='bold')

def cal_season(ds,varname='prec'):
    ds.loc[:,'MAM']=ds.loc[:,slice('%s3'%varname,'%s5'%varname)].sum(axis=1)
    ds.loc[:,'JJA']=ds.loc[:,slice('%s6'%varname,'%s8'%varname)].sum(axis=1)
    ds.loc[:,'SON']=ds.loc[:,slice('%s9'%varname,'%s11'%varname)].sum(axis=1)
    ds.loc[:,'DJF']=ds.loc[:,[varname+'12',varname+'1',varname+'2']].sum(axis=1)
    return ds                        

# calculate et relative contribution to prec in different provinces
# return the top 30
# Manually edited Kashmir on two csv file to shorten the name
def load_zonal_prec(type='absolute',time_scale='year',rank=30, source_region='TP',lc_type='all',et_data='GLEAM_v3.5a',prec_data='ERA5'):
    if lc_type=='all':
        ds = pd.read_csv('../data/processed/prec_con_mon_%s_%s_zonal.csv'%(source_region,et_data))
    else:
        ds = pd.read_csv('../data/processed/prec_con_mon_%s_%s_%s_zonal.csv'%(source_region,lc_type,et_data))

    if (type=='absolute')&(time_scale=='year'):
        ds30=ds[['name','precYear']].set_index('name').sort_values(by='precYear',ascending=False)[:rank]
    if (type=='absolute')&(time_scale=='season'):
        ds=cal_season(ds)
        ds30=ds[['name','MAM','JJA','SON','DJF','precYear']].set_index('name').sort_values(by='precYear',ascending=False)[:rank]
    if (type=='relative')&(time_scale=='year'):
        dpz=pd.read_csv('../data/processed/prec_mon_%s_%s_zonal.csv'%(source_region,prec_data))
        ds30 = ((ds.set_index('name')['precYear']/dpz.set_index('name')['precYear']).replace(np.inf,np.nan).sort_values(ascending=False)).dropna().to_frame()[:rank]
    if (type=='relative')&(time_scale=='season'):
        dpz=pd.read_csv('../data/processed/prec_mon_%s_%s_zonal.csv'%(source_region,prec_data))
        ds=cal_season(ds)
        dpz=cal_season(dpz)
        ds30 = ((ds.set_index('name')[['MAM','JJA','SON','DJF','precYear']] \
                 /dpz.set_index('name')[['MAM','JJA','SON','DJF','precYear']]) \
                .replace(np.inf,np.nan).sort_values(by='precYear',ascending=False)).dropna()[:rank]
    return ds30

# calculate et relative contribution to prec 
def load_prec_conptc(prec_data='ERA5',et_data='GLEAM_v3.5a'):
    dep=xr.open_dataset('../data/processed/utrack_climatology_prec_0.5_mon_TP_%s.nc'%et_data)
    if prec_data=='CMFD':
        dp=xr.open_dataset('../data/prec_CMFD_V0106_B-01_01mo_050deg_2008-2017_ymonmean_clean.nc')
    if prec_data=='ERA5':
        dp=load_prec_region(prec_data=prec_data)
    prec_conpct = dep.prec.sum(dim='month')/dp.prec.sum(dim='month')
    return prec_conpct

# get China mask based CMFD data
def get_china_mask():
    p=load_prec_region(prec_data='CMFD')
    return p.prec.sum(dim='month')>0

def make_plot(prec_data='ERA5',et_data='GLEAM_v3.5a'):
    # load data
    dp = xr.open_dataset('../data/processed/utrack_climatology_prec_0.5_mon_TP_%s.nc'%et_data)# prec contribution
    dpct=load_prec_conptc(prec_data=prec_data,et_data=et_data) # relative prec contribution
    ds30=load_zonal_prec(type='absolute',et_data=et_data) # top 30 provincial prec contribution
    dsr = load_zonal_prec(type='relative',prec_data=prec_data,et_data=et_data) *100 # top 30 provincial prec relative contribution
    cn_mask=get_china_mask() 
    # create four colors for barchart
    mybincolor1=bin_color(ds30,[0, 10, 50, 100,1000])
    mybincolor2=bin_color(dsr,[0, 2, 5, 20, 1000])
    
    # create levels for maps
    if et_data=='GLEAM_v3.5a':
        levels1=[0,1,5,10,20,50,100,200,300,400,500]
    else:
        levels1=[0,1,5,10,20,50,100,200,300,400,500,600]
    if prec_data=='ERA5':
        levels2=[0,0.01,0.05,0.1,0.20,0.40,0.6,0.8]
    else:
        levels2=[0,0.01,0.05,0.1,0.20,0.50,0.8,1]
    
    # create uneven levels for cmap for maps
    mycmap1,mynorm1=uneven_cmap(levels1) # for panel a
    mycmap2,mynorm2=uneven_cmap(levels2) # for panel c
    
    
    ###########Panel A
    fig = plt.figure(figsize=[10,8])
    ax1 = fig.add_axes([0.075, 0.55, 0.4, 0.4], projection=ccrs.PlateCarree(),
                                         frameon=False)
    
    # only plot regions with prec contribution > 1 mm
    ma=dp['prec'].sum(dim='month')>1 
    
    plot_map(dp['prec'].sum(dim='month').where(ma), ax1, levels1, lw=1)
    set_lat_lon(ax1, range(70,140,20), range(10,51,20), label=True, pad=0.05, fontsize=10)
    
    # Add colorbar to big plot
    cbarbig1_pos = [ax1.get_position().x0, ax1.get_position().y0-0.03, ax1.get_position().width, 0.02]
    caxbig1 = fig.add_axes(cbarbig1_pos)
    
    cbbig1 = mpl.colorbar.ColorbarBase(ax=caxbig1, cmap=mycmap1, norm=mynorm1, orientation='horizontal',
                                      ticks=levels1)
    cbbig1.ax.set_yticklabels(levels1,fontsize=10)
    cbbig1.set_label('Precipitation contribution (mm)')
    
    ######################### Panel B: precipitation contribution in different provinces
    ax2 = fig.add_axes([0.575, 0.6, 0.4, 0.325],frameon=True)
    
    ds30.plot(kind='bar', color=mybincolor1, ax=ax2, legend=False, edgecolor='k',linewidth=0.75)
    ax2.set_ylabel('Precipitation contribution (mm)')
    ax2.set_xlabel('')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    
    
    ###################### Panel C: relative contribution map
    ax3 = fig.add_axes([0.075, 0.075, 0.4, 0.4], projection=ccrs.PlateCarree(),
                                         frameon=False)
    
    plot_map(dpct.where(cn_mask), ax3,levels2, lw=1)
    set_lat_lon(ax3, range(70,140,20), range(10,51,20), label=True, pad=0.05, fontsize=10)
    
    # Add colorbar to big plot
    cbarbig2_pos = [ax3.get_position().x0, ax3.get_position().y0-0.03, ax3.get_position().width, 0.02]
    caxbig2 = fig.add_axes(cbarbig2_pos)
    
    cbbig2 = mpl.colorbar.ColorbarBase(ax=caxbig2, cmap=mycmap2, norm=mynorm2, orientation='horizontal',
                                      ticks=levels2)
    cbbig2.ax.set_xticklabels((np.array(levels2)*100).astype(np.int),fontsize=10)
    cbbig2.set_label('Precipitation contribution (%)')
    
    ################### Panel D: relative precipitation contribution in different provinces
    ax4 = fig.add_axes([0.575, 0.125, 0.4, 0.325],frameon=True)
    
    dsr.plot(kind='bar', color=mybincolor2, ax=ax4, legend=False, edgecolor='k',linewidth=0.75)
    ax4.set_ylabel('Precipitation contribution (%)')
    ax4.set_xlabel('')
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)
    
    # add panel label
    plot_subplot_label(ax1, 'a', left_offset=-0.1, upper_offset=0.125)
    plot_subplot_label(ax2, 'b', left_offset=-0.05,upper_offset=0.05)
    plot_subplot_label(ax3, 'c', left_offset=-0.1,upper_offset=0.125)
    plot_subplot_label(ax4, 'd', left_offset=-0.05,upper_offset=0.05)
    
    plt.savefig('../figure/figure_prec_con_map_%s_%s_0321.png'%(prec_data,et_data),dpi=300)
    print('figure saved')
   
if __name__=="__main__":
    make_plot()
