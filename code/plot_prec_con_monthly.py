import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

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
    f, axes = plt.subplots(3,2,figsize=(15, 12))
    
    # Generate a custom diverging colormap
   # cmap = sns.color_palette("Blues",n_colors=64)
    cmap = sns.color_palette("YlGnBu",n_colors=64)
    
    sns.heatmap(temp.loc[northwest_order],
                cmap=cmap, annot=True,
                linewidths=.5, cbar_kws={"extend":'both',"shrink": .75,"label":'Precipitation contribution (mm/mon)'},
               # vmax=10,xticklabels=range(1,13),ax=axes[0,0]) #,fmt=".2g" two sig figure
                vmax=6,xticklabels=range(1,13),ax=axes[0,0],fmt=".2g")
    
    sns.heatmap(temp.loc[northeast_order], cmap=cmap, annot=True,
                linewidths=.5, cbar_kws={"extend":'both',"shrink": .75,"label":'Precipitation contribution (mm/mon)'},
                vmax=1/2,xticklabels=range(1,13),ax=axes[0,1],fmt=".1f")
    
    sns.heatmap(temp.loc[north_order], cmap=cmap, annot=True,
                linewidths=.5, cbar_kws={"extend":'both',"shrink": .75,"label":'Precipitation contribution (mm/mon)'},
               vmax=2/2,xticklabels=range(1,13),ax=axes[1,0],fmt=".1f")
    
    sns.heatmap(temp.loc[east_order], cmap=cmap, annot=True,
                linewidths=.5, cbar_kws={"extend":'both',"shrink": .75,"label":'Precipitation contribution (mm/mon)'},
                vmax=2/2,xticklabels=range(1,13),ax=axes[1,1],fmt=".1f")
    
    sns.heatmap(temp.loc[southwest_order], cmap=cmap, annot=True,
                linewidths=.5, cbar_kws={"extend":'both',"shrink": .75,"label":'Precipitation contribution (mm/mon)'},
                vmax=20/2,xticklabels=range(1,13),ax=axes[2,0])
    
    sns.heatmap(temp.loc[central_south_order], cmap=cmap, annot=True,
                linewidths=.5, cbar_kws={"extend":'both',"shrink": .75,"label":'Precipitation contribution (mm/mon)'},
                vmax=2/2,xticklabels=range(1,13),ax=axes[2,1],fmt=".1f")
    
    # remove ylabel
    [a.set_ylabel('') for a in axes.flatten()]
    
    # add panel title
    [a.set_title(t,fontsize=12,fontweight='bold') for t,a in zip(['Northwest','Northeast','North','East','Southwest','South'],axes.flatten())]
    
    # disable tick label auto rotation
    axes[0,0].set_yticklabels(northwest_order,rotation=False)
    axes[1,0].set_yticklabels(north_order,rotation=False)
    axes[0,1].set_yticklabels(northeast_order,rotation=False)
    
    axes[2,0].set_xlabel('Month')
    axes[2,1].set_xlabel('Month')
    
    plt.subplots_adjust(hspace=0.3)
    
    plt.savefig('../figure/figure_prec_con_monthly_0915.png',dpi=300,bbox_inches='tight')
    print('figure saved')

if __name__=="__main__":
    make_plot(var='Et')
