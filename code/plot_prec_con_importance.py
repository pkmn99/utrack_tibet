import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

# Plot fitting lines 
def plot_fitting(df,x_txt, y_txt, ax, order=2):
    p = np.poly1d(np.polyfit(df[x_txt], df[y_txt], order))
    xp = np.linspace(df[x_txt].min(), df[x_txt].max(), 100)
    ax.plot(xp, p(xp), '--',lw=1,color='b')

# Plot scatter value label and manually adjust label location to reduce overlap
def plot_name(df, x_txt, y_txt, ax, fontsize=8,alpha=1):
    for i in df.index.values:
        va='bottom' # default
        ha='center'
        if np.isin(i, ['Liaoning','Henan','Guangxi','Shandong']):
            ha='right'
        if np.isin(i, ['Henan','Hebei','Hongkong','Guangxi']):
            va='top'
        if np.isin(i, ['Guangdong','Hongkong']):
            va='top'
            ha='left'
        if np.isin(i, ['Beijing','Jiangsu','Hunan','Chongqing','Ningxia','Shanxi']):
            ha='left'
        if np.isin(i, ['Jiangxi','Yunnan']):
            ha='left'
            va='top'
        ax.text(df.loc[i,x_txt],df.loc[i,y_txt], i, fontsize=fontsize,alpha=alpha,va=va,ha=ha)

def make_plot():
    d = pd.read_csv('../data/processed/table_prec_con_importance.csv')
    # begin plot
    fig, axes = plt.subplots(1,1,figsize=(8,6))
    # set range for zoomed part which do not label value in the main plot
    range_y=(d['Importance']>3.5)&(d['Importance']<3.8)
    range_x=(d['Annual']>0)&(d['Annual']<20)
    
    d.plot.scatter(x='Annual',y='Importance',ax=axes,c='b',edgecolor='k')
    plot_name(d[~(range_x&range_y)].set_index('name'), 'Annual','Importance',axes,fontsize=10)
    plot_fitting(d, 'Annual','Importance',axes,order=1)
    axes.set_xlabel('Precipitation contributed by TP (mm/yr)',fontsize=12)
    axes.set_ylabel('Perceived importance score of TP',fontsize=12)
    plt.text(250,4.1,'r=%.2f'%d.corr()['Annual']['Importance'], fontsize=12,color='b')
    
    # Add zoom region box
    axes.add_patch(Rectangle((0,3.5), 20,0.3, facecolor='none',edgecolor='r',lw=1))
    
    pos2 = [0.55, 0.2, 0.3, 0.35] # [left, bottom, width, height] #
    ax2 = fig.add_axes(pos2)
    d.plot.scatter(x='Annual',y='Importance',ax=ax2,c='b',edgecolor='k')
    
    plot_name(d[range_x&range_y].set_index('name'), 'Annual','Importance',ax2,fontsize=8)
    plot_fitting(d, 'Annual','Importance',ax2,order=1)
    ax2.set_xlabel('')
    ax2.set_ylabel('')
    
    ax2.set_xticks(np.arange(0, 21, 5))
    ax2.set_yticks(np.arange(3.5, 3.81, 0.1))
    
    ax2.set_ylim([3.5,3.8])
    ax2.set_xlim([0,20])
    
    # change color of border 
    # https://stackoverflow.com/questions/7778954/elegantly-changing-the-color-of-a-plot-frame-in-matplotlib
    plt.setp(ax2.spines.values(), color='r')
    
   # plt.savefig('../figure/plot_prec_con_importance.pdf',dpi=300,bbox_inches='tight')
    plt.savefig('../figure/plot_prec_con_importance.png',dpi=300,bbox_inches='tight')
    print('figure saved')

if __name__=='__main__':
    make_plot()
