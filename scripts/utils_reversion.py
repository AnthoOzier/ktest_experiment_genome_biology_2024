import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
from sklearn.cluster import KMeans


def get_meta_from_df(df):
    meta = pd.DataFrame()
    meta['index'] = df.index
    meta.index = df.index
    meta['batch'] = meta['index'].apply(lambda x : x.split(sep='.')[0])
    meta['condition'] = meta['index'].apply(lambda x : x.split(sep='.')[1])
    return(meta)


def split_48HREV(self,comparaison,t,threshold=None,nobs=None,orientation='>',verbose=0):
    assert(any([a is not None for a in [threshold,nobs]]))
    
    self.set_test_data_info(samples=comparaison,condition='condition')
    self.projections(t)
    pop = self.select_observations_from_condition(threshold=threshold,t=t,orientation=orientation,nobs=nobs,sample='48HREV')
    ref_sample = comparaison[0] if comparaison[0] != '48HREV' else comparaison[1]
    
    criterium = f'threshold{threshold}' if threshold is not None else f'nobs{nobs}'
    condition = f'pop_{ref_sample}_t{t}_{criterium}'

    self.split_sample(pop,'48HREV',new_condition=condition,condition='condition',verbose=verbose)
    self.set_test_data_info(samples=comparaison,condition='condition')
    
    return(pop,condition)

def custom_kde(data,
               normalize=True,
               coef=1,
               minmax=None,
               bw_method=.2,
               yshift=0,
               xshift=0
            ):
    if minmax is not None:
        min,max = minmax
        min,max = min - .1*(max-min),max +.1*(max-min)
    else:
        min,max = data.min(),data.max()
        min,max = min - .1*(max-min),max +.1*(max-min)

    
    x = np.linspace(min,max,200)
    density = gaussian_kde(data,bw_method=bw_method)
    y = density(x)
    if normalize : 
        y = y/np.max(y) 
    y = coef*y+yshift
    x = x+xshift
    return(x,y)


def custom_violin(data,
               fig=None,
               ax=None,
               orientation='vertical',
               alpha=.5,
               label=None,
               color=None, 
               lw=2,
               minmax=None,
               bw_method=.2,
               yshift=0,
               xshift=0
               ):
    if fig is None:
        fig,ax = plt.subplots(figsize=(12,6))
    
    for coef in [-1,1]:
        if coef == 1:
            label = None
        fig,ax,c = custom_rug(data,
                   fig=fig,
                   ax=ax,
                   orientation=orientation,
                   alpha=alpha,
                   label=label,
                   color=color, 
                   lw=lw,
                   normalize=True,
                   coef=coef,
                   minmax=minmax,
                   bw_method=bw_method,
                   yshift=yshift,
                   xshift=xshift
                   )
        if color is None: 
            color=c
    return(fig,ax,color)

def custom_rug(data,
               fig=None,
               ax=None,
               orientation='vertical',
               alpha=.5,
               label=None,
               color=None, 
               lw=2,
               normalize=True,
               coef=1,
               minmax=None,
               bw_method=.2,
               yshift=0,
               xshift=0
               ):
    if fig is None:
        fig,ax = plt.subplots(ncols=1,figsize=(12,6))

    x,y = custom_kde(data=data,
               normalize=normalize,
               coef=coef,
               minmax=minmax,
               bw_method=bw_method,
               yshift=yshift,
               xshift=xshift)

    if orientation == 'vertical':
        ax.plot(x,y,color=color,lw=lw,)
    else:
        ax.plot(y,x,color=color,lw=lw)

    if color is None:
        color = ax._children[-1]._color


    if orientation == 'vertical':
        ax.fill_between(x,y,y2=yshift,color=color,label=label,alpha=alpha)
    else:
        ax.fill_betweenx(x,y,x2=yshift,color=color,label=label,alpha=alpha)
        

    return(fig,ax,color)


def double_histogram(data,
                     fig=None,
                     ax=None,
                     alpha=.5,
                     orientation='vertical',
                     label=None,
                     color=None,
                     edgecolor=None,
                     coef_bins=3,
                     yshift=0
                     ):
    
    if fig is None:
        fig,ax = plt.subplots(ncols=1,figsize=(12,6))
    bins = coef_bins*int(np.floor(np.sqrt(len(data))))

    ax.hist(data,
        density=True,
        histtype='bar',
        bins=bins,
        alpha=alpha,
        orientation=orientation,
        label=label,color=color)
    
    if edgecolor is None:
        edgecolor=ax._children[-1]._facecolor
    if color is None:
        color = ax._children[-1]._facecolor

    ax.hist(data,
            density=True,
            histtype='step',
            bins=bins,
            lw=3,
            orientation=orientation,
            edgecolor=edgecolor)    

    return(fig,ax,color)

def custom_histogram(data,
                     fig=None,
                     ax=None,
                     orientation='vertical',
                     alpha=.5,
                     label=None,
                     color=None,
                     edgecolor=None,
                     lw=2,
                     coef_bins=3,
                     means=True,
                     hist_type='kde',
                     kde_bw=.2,
                     minmax=None,
                     xshift=0,
                     yshift=0,
                     normalize=True
                     ):
    """
    Parameters
    ----------
        type (default = 'kde') : in ['kde','violin','hist'] 
    """

    if hist_type=='kde':
        nz_data = data[data!=0]
        if len(nz_data)>0 and len(nz_data.unique())>1:
            fig,ax,color=custom_rug(data=data,
               fig=fig,
               ax=ax,
               orientation=orientation,
               alpha=alpha,
               label=label,
               color=color, 
               lw=lw,
               minmax=minmax,
               bw_method=kde_bw,
               xshift=xshift,
               yshift=yshift,
               normalize=normalize
               )

    if hist_type == 'violin':
        if not len(data[data==0])==len(data):
            fig,ax,color=custom_violin(data=data,
               fig=fig,
               ax=ax,
               orientation=orientation,
               alpha=alpha,
               label=label,
               color=color, 
               lw=lw,
               minmax=minmax,
               bw_method=kde_bw,
               xshift=xshift,
               yshift=yshift,
               )


    if hist_type == 'hist':
        fig,ax,color= double_histogram(data=data,
                     fig=fig,
                     ax=ax,
                     alpha=alpha,
                     orientation=orientation,
                     label=label,
                     color=color,
                     edgecolor=edgecolor,
                     coef_bins=coef_bins
                     )

    if means == 'line': 
        if orientation =='vertical':
            ax.axvline(data.mean(),c=color,lw=1.5)
        else:
            ax.axhline(data.mean(),c=color,lw=1.5)
    elif means:
        if orientation == 'vertical':
            ax.scatter(data.mean()+xshift,yshift,marker='+',s=100,c=color)
        else:
            ax.scatter(yshift,data.mean()+xshift,marker='+',s=100,color=color)            
    return(fig,ax)


def points_in_boxplot(df,ax,colors=None,vert=False):

    # ajouter les points 
    vals, ys = [], [] 
    for i,c in enumerate(df.columns):
        vals.append(df[c])
        ys.append(np.random.normal(i+1, 0.04, len(df)))
    ngroup = len(vals)
    clevels = np.linspace(0., 1., ngroup)
    colors = list(colors.values()) if colors is not None else [None]*len(ys)
    for x, val, color in zip(ys, vals, colors):
        if vert :
            ax.scatter(x,val, c=color, alpha=1)  
        else:
            ax.scatter(val,x,c=color,alpha=1)

def filled_boxplot(df,ax,colors=None,alpha=.5,vert=False):
    bp = df.boxplot(ax=ax,return_type='both',patch_artist=True,vert=vert)

    for i,patch in enumerate(bp[1]['boxes']):
        if colors is not None:
            color=list(colors.values())[i]
            patch.set(facecolor=color,edgecolor=color)

        patch.set(alpha=alpha,
                 fill=True,
                 linewidth=0)

def contours_boxplot(df,ax,colors=None,lw=3,vert=False):
    bp = df.boxplot(ax=ax,return_type='both',patch_artist=True,vert=vert)
    for i,patch in enumerate(bp[1]['boxes']):
        if colors is not None:
            color=list(colors.values())[i]
            patch.set(edgecolor=color)

        patch.set(alpha=1,
                 fill=False,
                 linewidth=lw)
    if colors is not None:
        for i,cap in enumerate(bp[1]['caps']):
            cap.set(color=list(colors.values())[i//2],linewidth=lw)

        for i,whisker in enumerate(bp[1]['whiskers']):
            whisker.set(color=list(colors.values())[i//2],linewidth=lw)
        
        
def custom_boxplot(df,colors=None,alpha=.5,lw=3,scatter=True,fig=None,ax=None,vert=True):
    if fig is None:
        fig,ax = plt.subplots(figsize=(20,7))
    
    if colors is not None:
        colors = {s:c for s,c in colors.items() if s in df} 
        df = df[list(colors.keys())]
    
    filled_boxplot(df=df,colors=colors,alpha=alpha,ax=ax,vert=vert)
    contours_boxplot(df=df,colors=colors,lw=lw,ax=ax,vert=vert)
    if scatter: 
        points_in_boxplot(df=df,colors=colors,ax=ax,vert=vert)
    return(fig,ax)



def cluster_genes_from_condition(self,condition,n_clusters,colors = {'0H':'xkcd:bright blue','24H':'xkcd:leaf green','48HDIFF':'crimson','48HREV':'xkcd:pale purple',
         '48HREV_1':'xkcd:aqua green','48HREV_2':'xkcd:light magenta'},ylim=(-2,2),vert=True,
         verbose=0,fig=None,axes=None):
    if fig is None:
        fig,axes = plt.subplots(ncols = n_clusters,figsize=(8*n_clusters,6))
    dfm = self.get_dataframe_of_means(condition=condition,samples='all',
                                      )
    kmeans = KMeans(
         init="random",
         n_clusters=n_clusters,
         n_init=10,
         max_iter=300,
         random_state=4)
    kmeans.fit(dfm)
    
    dfm = dfm.rename(columns={c:f'{c}_{condition}' for c in dfm.columns},inplace=False)
    self.update_var_from_dataframe(dfm)
    var = self.var

    col_cluster = f'kmeans_{condition}_nc{n_clusters}'
    var[col_cluster]=kmeans.labels_
    var[col_cluster] = var[col_cluster].astype('category')
    
    for cat in range(n_clusters):
        ax = axes[cat]
        df = var[var[col_cluster]==cat].copy()
        cols = [c for c in df.columns if condition in c ]
        df = df.loc[:,cols]
        df = df.rename(columns={c:c.replace(f'_{condition}','') for c in df.columns},inplace=False)
        custom_boxplot(df,colors,fig=fig,ax=ax,vert=vert)
        if verbose>0:
            if verbose >1 or len(df) <200:
                print(f'cluster {cat} ({len(df)}) :',",".join(df.index.tolist()))
        ax.set_title(f'cluster {cat} ({len(df)} genes) ',fontsize=30)
        if vert:
            ax.set_ylim(ylim)
        else:
            ax.set_xlim(ylim)
    return(col_cluster,fig,axes)