import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import cartopy.io.shapereader as shpreader
import cartopy.crs as ccrs
from cartopy.feature import ShapelyFeature

# heatmap figure to show monthly precipitation contribution in different provinces

# return region rank
def region_rank(regionlist,df):
    return df.set_index('name').dropna().loc[regionlist].sort_values(by='precYear',ascending=False).index.values

def make_plot(et_data='GLEAM_v3.5a',var='E'):
    # load data
    depz=pd.read_csv('../data/processed/prec_con_mon_TP_%s_%s_zonal.csv'%(et_data,var))
    temp = depz.set_index('name').dropna().sort_values(by='precYear',ascending=False).iloc[:,0:12]
    
    # region list
    northwest=['Shaanxi','Gansu','Qinghai','Ningxia','Xinjiang']
    north=['Beijing','Tianjin','Hebei','Shanxi','Neimeng']
    northeast=['Liaoning','Jilin','Heilongjiang']
    east=['Shanghai','Jiangsu','Zhejiang','Anhui','Fujian','Jiangxi','Shandong','Taiwan']
    central_south=['Henan','Hubei','Hunan','Guangdong','Guangxi','Hainan','Hongkong']# ,'Macao'
    southwest=['Chongqing','Sichuan','Guizhou','Yunnan','Xizang']
    
    # order by precipitation contribution
    northwest_order=region_rank(northwest,depz)
    north_order=region_rank(north,depz)
    northeast_order=region_rank(northeast,depz)
    east_order=region_rank(east,depz)
    central_south_order=region_rank(central_south,depz)
    southwest_order=region_rank(southwest,depz)
    
    
    # Set up the matplotlib figure
    fig, axes = plt.subplots(3,2,figsize=(15, 12))
    
    # Generate a custom diverging colormap
   # cmap = sns.color_palette("Blues",n_colors=64)
    cmap = sns.color_palette("YlGnBu",n_colors=64)
    
    sns.heatmap(temp.loc[northwest_order],
                cmap=cmap, annot=True,
                linewidths=.5, cbar_kws={"extend":'both',"shrink": .75,"label":'$P_{ET}$ (mm/mon)'},
               # vmax=10,xticklabels=range(1,13),ax=axes[0,0]) #,fmt=".2g" two sig figure
                vmax=6,xticklabels=range(1,13),ax=axes[0,0],fmt=".2g")
    
    sns.heatmap(temp.loc[northeast_order], cmap=cmap, annot=True,
                linewidths=.5, cbar_kws={"extend":'both',"shrink": .75,"label":'$P_{ET}$ (mm/mon)'},
                vmax=1/2,xticklabels=range(1,13),ax=axes[0,1],fmt=".1f")
    
    sns.heatmap(temp.loc[north_order], cmap=cmap, annot=True,
                linewidths=.5, cbar_kws={"extend":'both',"shrink": .75,"label":'$P_{ET}$ (mm/mon)'},
               vmax=2/2,xticklabels=range(1,13),ax=axes[1,0],fmt=".1f")
    
    sns.heatmap(temp.loc[east_order], cmap=cmap, annot=True,
                linewidths=.5, cbar_kws={"extend":'both',"shrink": .75,"label":'$P_{ET}$ (mm/mon)'},
                vmax=2/2,xticklabels=range(1,13),ax=axes[1,1],fmt=".1f")
    
    sns.heatmap(temp.loc[southwest_order], cmap=cmap, annot=True,
                linewidths=.5, cbar_kws={"extend":'both',"shrink": .75,"label":'$P_{ET}$ (mm/mon)'},
                vmax=20/2,xticklabels=range(1,13),ax=axes[2,0])
    
    sns.heatmap(temp.loc[central_south_order], cmap=cmap, annot=True,
                linewidths=.5, cbar_kws={"extend":'both',"shrink": .75,"label":'$P_{ET}$ (mm/mon)'},
                vmax=2/2,xticklabels=range(1,13),ax=axes[2,1],fmt=".1f")
    
    # remove ylabel
    [a.set_ylabel('') for a in axes.flatten()]
    
    # add panel title
    region_color_order=[0,5,2,4,1,3] # order to match the map
   # [a.set_title(t,fontsize=12,fontweight='bold') for t,a in zip(['Northwest','Northeast','North','East','Southwest','South'],axes.flatten())]
    for i,t in enumerate(['Northwest','Northeast','North','East','Southwest','South']):
        axes.flatten()[i].set_title(t,fontsize=15,fontweight='bold',color=sns.color_palette()[region_color_order[i]])
    
    # disable tick label auto rotation
    axes[0,0].set_yticklabels(northwest_order,rotation=False)
    axes[1,0].set_yticklabels(north_order,rotation=False)
    axes[0,1].set_yticklabels(northeast_order,rotation=False)

    # set ytick color of different regions
   # for i,t in enumerate([northwest_order,northeast_order,north_order,east_order,southwest_order,central_south_order]):
   #     axes.flatten()[i].set_yticklabels(t,rotation=False,color=sns.color_palette()[region_color_order[i]])
    
    axes[2,0].set_xlabel('Month',fontsize=14)
    axes[2,1].set_xlabel('Month',fontsize=14)

    ############### Inset to show China region division
    ax2inset = fig.add_axes([0.4, 0.85, 0.15, 0.125], projection=ccrs.PlateCarree(),
                                                       frame_on=False)
    ax2inset.outline_patch.set_visible(False) # Turn off borader

    # Load geographical data
    for i,r in enumerate(['西北','西南','华北','华南','East','东北']):
        china_shp=shpreader.Reader('/media/liyan/HDD/Project/data/China_gis/七大分区/%s.shp'%r)
        china_feature = ShapelyFeature(china_shp.geometries(), ccrs.PlateCarree(), facecolor=sns.color_palette()[i])
        ax2inset.add_feature(china_feature,edgecolor='w', linewidth=0.1)
    ax2inset.set_extent([70, 140, 10, 50],ccrs.Geodetic())
    
    plt.subplots_adjust(hspace=0.3)
    
    plt.savefig('../figure/figure_prec_con_monthly_%s_1003.png'%var,dpi=300,bbox_inches='tight')
    print('figure saved')

if __name__=="__main__":
    make_plot(var='Et')
