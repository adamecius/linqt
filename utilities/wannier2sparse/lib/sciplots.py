import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator

mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.rc('text', usetex=True)


from cycler import cycler
default_cycler = (cycler(color=['k', '#a90308', '#056eee','#5fa052']) +
                  cycler(linestyle=['-', '--', '-.', ':'])+
                  cycler(linewidth=[3, 2,2,3]))
mpl.rcParams['axes.prop_cycle'] = default_cycler
mpl.rcParams.update({'font.size': 18})


def load_data( fnames, **kwargs ):
    """
    
    If sahr
    """
    fnames= [fnames] if isinstance(fnames, str ) else fnames;
    
    Xs,Ys = [],[];
    for f in fnames:
        X,Y = np.loadtxt(f,unpack=True);
        Xs.append(X);
        Ys.append(Y);
    

    #If xrange == (xmin,xmax)
    #Then the Xs and Ys will be filtered based on this value
    key="xrange"
    if( key in kwargs  ):
        xmin,xmax = kwargs[key];
        for i,X in enumerate(Xs):
            idx = (X>=xmin)*(X<=xmax)
            Xs[i]=Xs[i][idx];
            Ys[i]=Ys[i][idx];
              
    
    #If sharex == True, Xs returns a single shared X array
    key="sharex"
    if( key in kwargs and kwargs[key]):
        X0 = Xs[0]; 
        if all( [np.array_equal(X0,X) for X in Xs] ):
            Xs=X0;
        else:
            print("The x-axis are not equivalent, cannot be combined")

    

    return Xs,np.array(Ys);

def plot(x,y, ax=None,**kwargs):

    if ax is None:
        ax = plt.gca();  
        
    ax.margins(x = 0.00) # 5% padding in all directions

    ts =1
    key="textscale";
    if( key in kwargs  ):  
        ts = kwargs[key];
        
    
    label=""
    key="label";
    if( key in kwargs  ):  
        label = kwargs[key];

    
    labelfs =int(18*ts);

    key="xlabel";
    if( key in kwargs  ):  
        xlabel = kwargs[key];
        ax.set_xlabel( r"$ "+xlabel+" $", fontsize=labelfs );

    key="ylabel";
    if( key in kwargs  ):  
        ylabel = kwargs[key];
        ax.set_ylabel( r"$ "+ylabel+" $", fontsize=labelfs);
        

    ax.plot(x,y, label=label);

    tickfs =int(14*ts);

    ax.xaxis.set_tick_params(labelsize=tickfs)
    ax.yaxis.set_tick_params(labelsize=tickfs)
    return ax;

def add_legends( ax, legends ):
    
    lines = [ line for line in ax.lines ]
    for i,legend in enumerate(legends):
        lines[i].set_label(r"$"+legend+" $")
        ax.legend(fontsize=int(18), frameon=False,labelspacing=0.2, handlelength=1,borderpad=0.2,handletextpad=0.2);
        
        
        
def plot4fig( Xs, Ys, **kwargs ):
    #Since this will plot 4 figure, it is assume Ys has len 4
    num_figs = 4;
    assert len(Ys) == num_figs, "The dataset submited through Ys should have length 4"
    Ys = np.array(Ys);
    
    
    key = "shareX";
    sharex = kwargs[key] if key in kwargs else False;
    if sharex:
        skey = "Xlabels";
        if skey in kwargs:
            labels = kwargs[skey];
            labels= np.broadcast_to( labels, (num_figs, len(labels)) );
        Xs    = np.broadcast_to( Xs    , (num_figs, len(Xs)) )

    key = "fig";
    fig, axs = None, None;
    if key in kwargs:
        fig = kwargs[key];
        axs = fig.axes;
    else:
        fig, axs = plt.subplots( 2, 2, dpi=300 );
        axs = axs.flatten();

    key = "Ylabels";
    if key in kwargs:
        for (ax,label) in zip(axs,kwargs[key]):
            ax.set_ylabel(label);

    key = "Xlabels";
    if key in kwargs:
        for (ax,label) in zip(axs,kwargs[key]):
            ax.set_xlabel(label);

    for ax,X,Ylist in zip(axs,Xs,Ys):
        Ylist = [Ylist,] if len(Ylist.shape) == 1 else Ylist;
        for Y in Ylist:
            ax.plot(X,Y);
        ax.margins(x = 0.00) # 5% padding in all directions
    
    plt.tight_layout(pad=0.5, h_pad=None, w_pad=None);
    return fig;

def plot_ZwithArrows(XX,YY,ZZ,UU,VV, Zlabel="", mask = None,npoints=(3,3), ax=None):
    if mask is not None:
        UU,VV,ZZ = [ np.ma.array(XX, mask=mask) for XX in (UU,VV,ZZ) ];

    # Interpolate to regularly-spaced quad grid.
    x,y = XX[:,0].flatten(),YY[0].flatten();
    fzi,fui,fvi =[ RegularGridInterpolator((x,y),  V) for V in [ZZ,UU,VV] ];

    #numdivisions
    kgrid= np.meshgrid(*[ np.linspace(x.min(),x.max(),n) for n,x in zip(npoints,[XX,YY]) ],indexing='ij');
    xi,yi= kgrid;
    kgrid= np.transpose(kgrid).swapaxes(0,1);
    
    zi,ui,vi =[ fx(kgrid) for fx in (fzi,fui,fvi)];

    if mask is not None:
        UU,VV,ZZ = [ np.ma.array(XX, mask=mask(ZZ)) for XX in (UU,VV,ZZ) ];
        ui,vi,zi = [ np.ma.array(xi, mask=mask(zi)) for xi in (ui,vi,zi) ];

    plotcb = False;
    if ax is  None:
        ax = plt.gca();
        plotcb = True;

    cos = ax.contour (XX,YY, ZZ,levels=[0], cmap="bone")
    cs = ax.contourf(XX,YY, ZZ,levels=100, cmap="bone")
    Qv = ax.quiver(xi,yi,ui,vi, color="C1",  scale=10.0, width=0.021/2.0);
    for c in cs.collections:
        c.set_edgecolor("face")
    for c in cos.collections:
        c.set_edgecolor("white")
    if plotcb:
        cbar = plt.gcf().colorbar(cs, shrink=0.9,label=Zlabel);
    return 0;

