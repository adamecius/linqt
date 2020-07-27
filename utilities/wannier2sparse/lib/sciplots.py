import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.rc('text', usetex=True)


from cycler import cycler
default_cycler = (cycler(color=['k', '#a90308', '#056eee','#5fa052']) +
                  cycler(linestyle=['-', '--', '-.', ':'])+
                  cycler(linewidth=[3, 2,2,3]))
mpl.rcParams['axes.prop_cycle'] = default_cycler



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
        ax.legend(fontsize=int(18), frameon=False,labelspacing=0.2, handlelength=1,borderpad=0.2,handletextpad=0.2)
