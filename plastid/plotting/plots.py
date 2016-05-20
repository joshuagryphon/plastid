#!/usr/bin/env python
"""Sommon common plots that are not directly implemented in :mod:`matplotlib`,
as well as some specific plots used by :mod:`plastid`. Demos of these appear
in :doc:`/examples/z_plotting`.

General plots
-------------
    :func:`stacked_bar`
        Create a stacked bar graph

    :func:`kde_plot`
        Plot a kernel density estimate (continuous histogram) of data

    :func:`scatterhist_x`, :func:`scatterhist_y`, and :func:`scatterhist_xy`
        Create scatter plots with kernel density estimates of the marginal 
        distributions along side them

    :func:`profile_heatmap`
        Plot a heatmap, with a columnwise median (or other profile summarizing
        the heatmap) in an upper panel above it, with aligned coordinates

    :func:`triangle_plot`
        Plot data lying on the plane x + y + z = k (e.g. codon phasing)
        in a homogeneous 2D representation


Plots for genomics
------------------
    :func:`ma_plot`
        Plot fold changes between `x` and `y` (:math:`\log_{2} (y/x)`) as
        a function of the mean of x and y (:math:`0.5*(x+y)`).

    :func:`phase_plot`
        For ribosome profiling. Plot sub-codon phasing of ribosome-protected
        foorpints stratified by read length, as well as the fraction of total
        reads represented by each length.

"""
import copy
import numpy
import scipy.stats
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.artist import Artist
from plastid.plotting.colors import lighten, darken, process_black
from plastid.plotting.plotutils import get_fig_axes, split_axes, clean_invalid, get_kde



#==============================================================================
# Default keyword arguments for plots or subplots
# and helper functions
#==============================================================================

plastid_default_scatter = {
    "marker"     : "o",
    "alpha"      : 0.7,
    "facecolor"  : "none",
    "s"          : 8,
    "rasterized" : True,
}
"""Default parameters for scatter plots"""


def get_color_cycle(ax):
    """Get color cycle iterator from multiple versions of matplotlib axes

    Parameters
    ----------
    ax : :class:`matplotlib.axes.Axes`


    Returns
    -------
    iterator
        Iterator over colors, passable to matplotlib `color` keyword
    """
    if hasattr(ax,"_get_lines"):
        if hasattr(ax._get_lines,"prop_cycler"):
            return (X["color"] for X in ax._get_lines.prop_cycler)
        elif hasattr(ax._get_lines,"color_cycle"):
            return ax._get_lines.color_cycle
    else:
        raise AssertionError("get_color_cycle: Could not find color cycle")




#==============================================================================
# Stacked bar
#==============================================================================

def stacked_bar(data,axes=None,labels=None,lighten_by=0.1,cmap=None,**kwargs):
    """Create a stacked bar graph
    
    Parameters
    ----------
    data : :class:`numpy.ndarray`
        Array of data, in which each row is a stack, each column a value in that stack.

    axes : :class:`matplotlib.axes.Axes` or `None`, optional
        Axes in which to place plot. If `None`, a new figure is generated.
        (Default: `None`)
        
    labels : list, optional
        Labels for each stack. If `None`, stacks are labeled sequentially by number.
        (Default: `None`)
    
    lighten_by : float, optional
        Amount by which to lighten sequential blocks in each stack. (Default: 0.10)
    
    cmap : :class:`matplotlib.colors.Colormap`, optional
        Colormap from which to generate bar colors. If supplied, will override
        any `color` attribute in `**kwargs`. (Default: `None`)
        
    **kwargs : keyword arguments
        Other keyword arguments to pass to :func:`matplotlib.pyplot.bar`

    Returns
    -------
    :class:`matplotlib.figure.Figure`
        Parent figure of axes
    
    :class:`matplotlib.axes.Axes`
        Axes containing plot
    """
    fig, ax = get_fig_axes(axes)
    rows, cols = data.shape
    labels = labels if labels is not None else range(rows)
    defaults = [("align","center"),
                ("width",0.8)]

    if cmap is not None:
        kwargs["color"] = cmap(numpy.linspace(0,1.0,num=10))
    elif kwargs.get("color",None) is None:
        kwargs["color"] = [next(get_color_cycle(ax)) for _ in range(rows)]
        
    x = numpy.arange(rows) + 0.5
    xaxis = ax.xaxis
    xaxis.set_ticks(x)
    xaxis.set_ticklabels(labels)
    bottoms = numpy.zeros(rows)
   
    for k,v in defaults:
        if k not in kwargs:
            kwargs[k] = v
    
    for i in range(cols):
        color = kwargs["color"]
        if i > 0:
            kwargs["color"] = lighten(color,amt=lighten_by)

        heights = data[:,i]
        plt.bar(x,heights,bottom=bottoms,**kwargs)
        heights.shape
        bottoms += heights
    
    ax.set_xlim(-0.5,rows+0.5)

    return fig, ax



#==============================================================================
# Kernel density estimate
#==============================================================================

def kde_plot(data,axes=None,color=None,label=None,alpha=0.7,vert=False,
            log=False,base=10,points=500,bw_method="scott"):
    """Plot a kernel density estimate of `data` on `axes`.

    Parameters
    ----------
    data : :class:`numpy.ndarray`
        Array of data

    axes : :class:`matplotlib.axes.Axes` or `None`, optional
        Axes in which to place plot. If `None`, a new figure is generated.
        (Default: `None`)
        
    color : matplotlib colorspec, optional
        Color to use for plotting (Default: use next in matplotlibrc)

    label : str, optional
        Name of data series (used for legend; default: `None`)

    alpha : float, optional
        Amount of alpha transparency to use (Default: 0.7)

    vert : bool, optional
        If true, plot kde vertically

    log : bool, optional
        If `True`, `data` is log-transformed before the kde is estimated.
        Data are converted back to non-log space afterwards.

    base : 2, 10, or :obj:`numpy.e`, optional
        If `log` is `True`, this serves as the base of the log space.
        If `log` is `False`, this is ignored. (Default: 2)

    points : int
        Number of points over which to evaluate kde. (Default: 100)

    bw_method : str
        Bandwith estimation method. See documentation for
        :obj:`scipy.stats.gaussian_kde`. (Default: "scott")    


    Returns
    -------
    :class:`matplotlib.figure.Figure`
        Parent figure of axes
    
    :class:`matplotlib.axes.Axes`
        Axes containing plot
    """
    fig, axes = get_fig_axes(axes)

    if color is None:
        color = next(get_color_cycle(axes))

    a, b = get_kde(data,log=log,base=base,points=points,bw_method=bw_method)

    fbargs = { "alpha" : alpha,
               "facecolor" : lighten(color),
               "edgecolor" : color
             }
    if label is not None:
        fbargs["label"] = label

    if vert == True:
        axes.fill_betweenx(a,b,0,**fbargs)
        axes.plot(b,a,color=color,alpha=alpha,label=label) # this is a bit of a hack to get labels to print; fill_between doesn't work with legends
        if log == True:
            axes.semilogy()
    else:
        axes.fill_between(a,b,0,**fbargs)
        axes.plot(a,b,color=color,alpha=alpha,label=label)
        if log == True:
            axes.semilogx()

    return fig, axes


#==============================================================================
# Triangle plot
#==============================================================================

rotate = numpy.array([[0,1,0],
                      [0,0,1],
                      [1,0,0]])

_triA = numpy.array([[1,  0, 0],
                     [-1,-1, 1]])
_triT = numpy.array([[0.5,         1 ],
                     [0.5*(3**0.5),0]])
_triTA = _triT.dot(_triA)

_triverts = numpy.array([[1.0,0.0],
             [0.0,1.0],
             [0.0,0.0]])


def trianglize(data):
    """Convert points from triangular space to Cartesian space for plotting.
    Called internally by :func:`triangle_plot`.

    Parameters
    ----------
    data : class:`numpy.ndarray`
        Mx2 or Mx3 list or array of points in triangular space, where
        the first column is the first coordinate, the second column
        the second, and the third, the third.
    
    Returns
    -------
    :class:`numpy.ndarray`
        Corresponding Nx2 array of points in Cartesian space, for plotting
        on a standard axis
    """

    if data.shape[1] == 2:
        data = _triTA.dot(numpy.hstack([data,numpy.ones((data.shape[0],1))]).T).T
    else:
        data = _triT.dot(data[:,(0,2)].T).T 

    return data

def triangle_plot(data,axes=None,fn="scatter",vertex_labels=None,grid=None,clip=True,do_setup=True,**kwargs):
    """Plot data lying in a plane x + y + z = k in a homogenous triangular space.

    Parameters
    ----------
    data : :class:`numpy.ndarray`
        Mx2 or Mx3 list or array of points in triangular space, where
        the first column is the first coordinate, the second column
        the second, and the third, the third.

    axes : :class:`matplotlib.axes.Axes` or `None`, optional
        Axes in which to place plot. If `None`, a new figure is generated.
        (Default: `None`)
       
    fn : str, optional
        Name of plotting function. Must correspond to an attribute of a
        :class:`matplotlib.axes.Axes` (e.g. 'scatter', 'plot', 'hexbin'et c.),
        that is be able to take an Nx2 :class:`numpy.ndarray` in Cartesian space
        as input (e.g. 'plot', 'scatter', 'hexbin'; Default: 'scatter').

    vertex_labels : list or None, optional
        Labels for vertex. If `None`, vertices aren't labeled. (Default: `None`)

    grid : :class:`numpy.ndarray` or None, optional
        If not `None`, draw gridlines at intervals specified in `grid`,
        as long as the grid coordinate is > 0.33333 (center of triangle)
        and <= 1.0 (border).

    clip : bool, optional
        If `True` clipping masks corresponding to the triangle boundaries
        will be applied to all plot elements (Default: `True`)
        
    do_setup : bool, optional
        If `True`, the plot area will be prepared. A triangle will be drawn,
        gridlines drawn, et c. Specify `False` if plotting a second dataset
        ontop of an already-prepared axes (Default: `True`)

    **kwargs : keyword arguments
        Other keyword arguments to pass to function specified by `fn`.

    Returns
    -------
    :class:`matplotlib.figure.Figure`
        Parent figure of axes
    
    :class:`matplotlib.axes.Axes`
        Axes containing plot

    """
    fig, ax = get_fig_axes(axes)
    mplrc = matplotlib.rcParams

    if do_setup == True:
        triverts = trianglize(_triverts)
        tripatch =  matplotlib.patches.Polygon(triverts,
                                               closed=True,
                                               facecolor=mplrc["axes.facecolor"],
                                               edgecolor=mplrc["axes.edgecolor"],
                                               linewidth=mplrc["axes.linewidth"],
                                               zorder=-10
                                               )
        ax.add_patch(tripatch)

        # format axes
        ax.set_xlim((0,1))
        ax.set_ylim((0,_triverts[:,1].max()))
        ax.set_frame_on(False)
        ax.set_xticks([])
        ax.set_yticks([])

        # label vertices
        if vertex_labels is not None:
            l1,l2,l3 = vertex_labels

            tkwargs = { "fig"   : fig,
                        "units" : "points"
                    }
            p1trans = matplotlib.transforms.offset_copy(ax.transData,x=0,  y=8,  **tkwargs)
            p2trans = matplotlib.transforms.offset_copy(ax.transData,x=-10,y=-12,**tkwargs)
            p3trans = matplotlib.transforms.offset_copy(ax.transData,x=10, y=-12,**tkwargs)
            ax.text(triverts[0,0],triverts[0,1],l1,transform=p1trans)
            ax.text(triverts[1,0],triverts[1,1],l2,transform=p2trans)
            ax.text(triverts[2,0],triverts[2,1],l3,transform=p3trans)

        # add gridlines
        grid_kwargs = { K.replace("grid.","") : V for (K,V) in mplrc.items() if K.startswith("grid") }
        if grid is not None:
            grid = numpy.array(grid)
            remainders = (1.0 - grid)/2
            for i, r in zip(grid,remainders):
                if i >= 1.0/3:
                    points = [numpy.array([i,r,r])]
                    for _ in range(3):
                        points.append(rotate.dot(points[-1]))

                    points = numpy.array(points)
                    points = trianglize(points[:,[0,2]])

                    myline = matplotlib.lines.Line2D(points[:,0],
                                                     points[:,1],
                                                     **grid_kwargs)
                    ax.add_line(myline)
    

    # scale data
    data = trianglize(data)

    # plot data
    artists = []
    fn = getattr(ax,fn)
    res = fn(*zip(*data),**kwargs)
    if isinstance(res,Artist):
        artists.append(res)
    elif isinstance(res,list):
        artists.extend([X for X in res if isinstance(X,Artist)])

    # clip
    if clip == True:
        for artist in artists:
            artist.set_clip_path(tripatch)
            artist.set_clip_on(True)

    return fig, ax


#==============================================================================
# Heatmaps with profiles on top
#==============================================================================

def sort_max_position(data):
    """Generate indices that sort rows in `data` by column in which
    the row's maximal value is attained

    Parameters
    ----------
    data : :class:`numpy.ndarray`

    Returns
    -------
    :class:`numpy.ndarray`
        Indices of rows that sort data by max position
    """
    maxvals = numpy.nanmax(data,1)
    maxidx  = numpy.zeros(len(maxvals))
    for i,maxval in enumerate(maxvals):
        maxidx[i] = (data[i,:] == maxval).argmax()

    return numpy.argsort(maxidx)

_heatmap_defaults = {
    "aspect"        :  "auto",
    "origin"        : "upper",
    "interpolation" : "none",
}

def profile_heatmap(data,profile=None,x=None,axes=None,sort_fn=sort_max_position,
                    cmap=None,nancolor="#666666",im_args={},plot_args={}):
    """Create a dual-paned plot in which `profile` is displayed in a top
    panel, above a heatmap showing the intensities of each row of `data`,
    optionally sorted top-to-bottom by `sort_fn`.

    Parameters
    ----------
    data : :class:`numpy.ndarray`
        Array of data, in which each row is an individual aligned vector of data
        for region of interest, and each column a position in that vector

    profile : :class:`numpy.ndarray` or None
        Reduced profile of data, often a column-wise median. If not
        supplied, it will be calculated.

    x : :class:`numpy.ndarray`
        Array of values for X-axis

    axes : :class:`matplotlib.axes.Axes` or `None`, optional
        Axes in which to place plot. If `None`, a new figure is generated.
        (Default: `None`)
       
    sort_fn : function, optional
        Sort rows in `data` by this function before plotting. Function must
        return a :class:`numpy.ndarray` of indices corresponding to rows in `data`
        (Default: sort by ascending argmax of each row)

    cmap : :class:`~matplotlib.colors.Colormap`, optional
        Colormap to use in heatmap. It not `None`, overrides any value
        in `im_args`. (Default: `None`) 

    nancolor : str or matplotlib colorspec
        Color used for plotting `nan` and other illegal or masked values
        
    im_args : dict
        Keyword arguments to pass to :func:`matplotlib.pyplot.imshow`

    plot_args : dict
        Keyword arguments to pass to :func:`matplotlib.pyplot.plot`
        for plotting metagene average


    Returns
    -------
    :class:`matplotlib.figure.Figure`
        Parent figure of axes
    
    dict
        Dictionary of :class:`matplotlib.axes.Axes`. "top" refers to the 
        panel containing the summary profile. "main" refers to the heatmap
        of individual values
    """
    fig, ax = get_fig_axes(axes)
    axes = split_axes(ax,top_height=0.2)

    if sort_fn is None:
        sort_indices = numpy.arange(data.shape[0])
    else:
        sort_indices = sort_fn(data)

    if x is None:
        x = numpy.arange(0,data.shape[1])

    if profile is None:
        profile = numpy.nanmedian(data,axis=0)
    
    im_args     = copy.deepcopy(im_args)
    
    # populate with defaults
    for k,v in _heatmap_defaults.items():
        if k not in im_args:
            im_args[k] = v

    if "extent" not in im_args:            
        im_args["extent"] = [x.min(),x.max(),0,data.shape[0]]#,0]
    if "vmin" not in im_args:
        im_args["vmin"] = numpy.nanmin(data)
    if "vmax" not in im_args:
        im_args["vmax"] = numpy.nanmax(data)

    if cmap is not None:
        im_args["cmap"] = cmap
    elif "cmap" in im_args:
        cmap = matplotlib.cm.get_cmap(im_args["cmap"])
    else:
        cmap = matplotlib.cm.get_cmap()
    
    cmap.set_bad(nancolor,1.0)
    
    axes["top"].plot(x,profile,**plot_args)
    axes["top"].set_ylim(0,profile.max())
    axes["top"].set_xlim(x.min(),x.max())
    #axes["top"].set_yticks([])
    axes["top"].set_yticks([0,profile.max()])
    axes["top"].xaxis.tick_bottom()
    axes["top"].grid(True,which="both")

    axes["main"].xaxis.tick_bottom()
    axes["main"].imshow(data[sort_indices,:],**im_args)

    return fig, axes


#==============================================================================
# Scatter plots with marginal distributions
#==============================================================================


def _scatterhist_help(axes=None,
                      top_height=0,left_width=0,right_width=0,bottom_height=0,
                      ):
    """Create a scatter plot with the marginal distribution for `x`
    plotted in a separate pain as a kernel density estimate.

    Parameters
    ----------
    color : matplotlib colorspec, or None, optional
        Color to plot data series

    axes : :class:`matplotlib.axes.Axes`, dict, or `None`, optional
        If a :class:`matplotlib.axes.Axes`, an axes in which to place plot.
        This axes will be split into relevant panels.

        If a dict, this is expected to be a group of pre-split axes.

        If `None`, a new figure is generated, and axes are split.
        (Default: `None`)
      
    top_height, left_width, right_width, bottom_height : float, optional
        If not `None`, a panel on the corresponding side of the `ax` will
        be created, using whatever fraction is specified (e.g. 0.1 to use
        10% of total height).    


    Returns
    -------
    :class:`matplotlib.figure.Figure`
        Parent figure of axes
    
    dict of :class:`matplotlib.axes.Axes`
        axes containing plot space
    """
    if axes is None:
        fig = plt.figure()
        axes = plt.gca()
    if isinstance(axes,matplotlib.axes.Axes):
        fig = axes.figure
        axes = split_axes(axes,top_height=top_height,left_width=left_width,
                          right_width=right_width,bottom_height=bottom_height)
    elif isinstance(axes,dict):
        fig = axes["main"].figure
        if left_width > 0:
            assert "left" in axes
        if right_width > 0:
            assert "right" in axes
        if bottom_height > 0:
            assert "bottom" in axes
        if top_height > 0:
            assert "top" in axes

    return fig, axes


def scatterhist_x(x,y,color=None,axes=None,label=None,
                  top_height=0.2,mask_invalid=True,
                  log="",
                  min_x=-numpy.inf,min_y=-numpy.inf,max_x=numpy.inf,max_y=numpy.inf,
                  scargs=plastid_default_scatter,bw_method="scott",
                  kdalpha=0.7):
    """Produce a scatter plot with a kernel density estimate of the marginal `x` distribution

    Parameters
    ----------
    x, y : :class:`numpy.ndarray` or list
        Pair arrays or lists of corresponding numbers

    color : matplotlib colorspec, optional
        Color to use in plot

    label : str, or None
        If not None, a label for plotting

    axes : :class:`matplotlib.axes.Axes`, dict, or `None`, optional
        If a :class:`matplotlib.axes.Axes`, an axes in which to place plot.
        This axes will be split into relevant panels.

        If a dict, this is expected to be a group of pre-split axes.

        If `None`, a new figure is generated, and axes are split.
        (Default: `None`)

    top_height : float, optional
        fraction of `axes` height to use in top panel containing
        marginal distribution (Default: 0.2)

    mask_invalid : bool, optional
        If `True` mask out any `nan`s or `inf`s, as these mess up kernel density
        estimates and histograms in matplotlib, even if in masked arrays

    min_x, min_y, max_x, max_y : number, optional
        If supplied, set values below `min_x` to `min_x`, values larger
        than `max_x` to `max_x` and so for `min_y` and `max_y`

    log : str, "", "x", "xy", or "xy", optional
        Plot these axes on a log scale (Default: "" for no log axes)

    scargs : Keyword arguments, optional
        Arguments to pass to :func:`~matplotlib.pyplot.scatter`
        (Default: :obj:`plastid_default_scatter`). We highly recommend
        setting `rasterized` to `True`!

    kdalpha : float, optional
        Alpha level (transparency) for marginal distributions (Default: 0.7)

    bw_method : str
        Bandwith estimation method. See documentation for
        :obj:`scipy.stats.gaussian_kde`. (Default: `"scott"`)


    Returns
    -------
    :class:`matplotlib.figure.Figure`
        Figure

    dict
        Dictionary of axes. `'orig'` refers to `ax`. The central panel is `'main'`.
        Other panels will be mapped to `'top'`, `'left'` et c, if they are created.
    """
    fig, axes = _scatterhist_help(axes=axes,top_height=top_height)
    xlog = False

    if color is None:
        color = next(get_color_cycle(axes["main"]))
    
    if mask_invalid == True:
        x, y = clean_invalid(x,y,min_x=min_x,max_x=max_x,min_y=min_y,max_y=max_y)

    if label is not None:
        scargs = copy.deepcopy(scargs)
        scargs["label"] = label

    if "x" in log:
        axes["main"].semilogx()
        axes["top"].semilogx()
        xlog = True
        xmask = x > 0
    else:
        xmask = numpy.tile(True,x.shape)

    if "y" in log:
        axes["main"].semilogy()


    axes["main"].scatter(x,y,edgecolor=color,**scargs)

    # kernel densities
    kargs = { "color" : color,
              "alpha" : kdalpha,
              "bw_method" : bw_method,
             }

    if xmask.sum() > 0:
        kde_plot(x[xmask],log=xlog,axes=axes["top"],**kargs)

    return fig, axes


def scatterhist_y(x,y,color=None,axes=None,label=None,
                  right_width=0.2,mask_invalid=True,log="xy",
                  min_x=-numpy.inf,min_y=-numpy.inf,max_x=numpy.inf,max_y=numpy.inf,
                  scargs=plastid_default_scatter,bw_method="scott",
                  kdalpha=0.7):
    """Produce a scatter plot with a kernel density estimate of the marginal `y` distribution

    Parameters
    ----------
    x, y : :class:`numpy.ndarray` or list
        Pair arrays or lists of corresponding numbers

    color : matplotlib colorspec, optional
        Color to use in plot

    label : str, or None
        If not None, a label for plotting

    axes : :class:`matplotlib.axes.Axes`, dict, or `None`, optional
        If a :class:`matplotlib.axes.Axes`, an axes in which to place plot.
        This axes will be split into relevant panels.

        If a dict, this is expected to be a group of pre-split axes.

        If `None`, a new figure is generated, and axes are split.
        (Default: `None`)

    right_width : float, optional
        fraction of `axes` width to use in right panel containing
        marginal distribution (Default: 0.2)

    mask_invalid : bool, optional
        If `True` mask out any `nan`s or `inf`s, as these mess up kernel density
        estimates and histograms in matplotlib, even if in masked arrays

    min_x, min_y, max_x, max_y : number, optional
        If supplied, set values below `min_x` to `min_x`, values larger
        than `max_x` to `max_x` and so for `min_y` and `max_y`

    log : str, "", "x", "xy", or "xy", optional
        Plot these axes on a log scale (Default: "" for no log axes)

    scargs : Keyword arguments, optional
        Arguments to pass to :func:`~matplotlib.pyplot.scatter`
        (Default: :obj:`plastid_default_scatter`). We highly recommend
        setting `rasterized` to `True`!

    kdalpha : float, optional
        Alpha level (transparency) for marginal distributions (Default: 0.7)

    bw_method : str
        Bandwith estimation method. See documentation for
        :obj:`scipy.stats.gaussian_kde`. (Default: "scott")

    Returns
    -------
    :class:`matplotlib.figure.Figure`
        Figure

    dict
        Dictionary of axes. `'orig'` refers to `ax`. The central panel is `'main'`.
        Other panels will be mapped to `'top'`, `'left'` et c, if they are created.
    """
    fig, axes = _scatterhist_help(axes=axes,right_width=right_width)
    ylog = False

    if color is None:
        color = next(get_color_cycle(axes["main"]))

    if label is not None:
        scargs = copy.deepcopy(scargs)
        scargs["label"] = label
    
    if mask_invalid == True:
        x, y = clean_invalid(x,y,min_x=min_x,max_x=max_x,min_y=min_y,max_y=max_y)

    if "x" in log:
        axes["main"].semilogx()

    if "y" in log:
        axes["main"].semilogy()
        axes["right"].semilogy()
        ylog = True
        ymask = y > 0
    else:
        ymask = numpy.tile(True,y.shape)

    axes["main"].scatter(x,y,edgecolor=color,**scargs)

    # kernel density
    kargs = { "color" : color,
              "alpha" : kdalpha,
              "bw_method" : bw_method,
             }
             
    if ymask.sum() > 0:
        kde_plot(y[ymask],log=ylog,axes=axes["right"],vert=True,**kargs)

    return fig, axes


def scatterhist_xy(x,y,color=None,axes=None,label=None,
                   top_height=0.2,right_width=0.2,mask_invalid=True,log="xy",
                   min_x=-numpy.inf,min_y=-numpy.inf,max_x=numpy.inf,max_y=numpy.inf,
                   scargs=plastid_default_scatter,
                   kdalpha=0.7,bw_method="scott"):
    """Produce a scatter plot with kernel density estimate of the marginal `x` and `y` distributions

    Parameters
    ----------
    x, y : :class:`numpy.ndarray` or list
        Pair arrays or lists of corresponding numbers

    color : matplotlib colorspec, optional
        Color to use in plot

    label : str, or None
        If not None, a label for plotting    top_height : float, optional
        fraction of `axes` height to use in top panel containing
        marginal distribution (Default: 0.2)
        
    right_width : float, optional
        fraction of `axes` width to use in right panel containing
        marginal distribution (Default: 0.2)

    top_height : float, optional
        fraction of `axes` height to use in top panel containing
        marginal distribution (Default: 0.2)

    axes : :class:`matplotlib.axes.Axes`, dict, or `None`, optional
        If a :class:`matplotlib.axes.Axes`, an axes in which to place plot.
        This axes will be split into relevant panels.

        If a dict, this is expected to be a group of pre-split axes.

        If `None`, a new figure is generated, and axes are split.
        (Default: `None`)

    mask_invalid : bool, optional
        If `True` mask out any `nan`s or `inf`s, as these mess up kernel density
        estimates and histograms in matplotlib, even if in masked arrays

    min_x, min_y, max_x, max_y : number, optional
        If supplied, set values below `min_x` to `min_x`, values larger
        than `max_x` to `max_x` and so for `min_y` and `max-y`

    log : str, "", "x", "xy", or "xy", optional
        Plot these axes on a log scale (Default: "" for no log axes)

    scargs : Keyword arguments, optional
        Arguments to pass to :func:`~matplotlib.pyplot.scatter`
        (Default: :obj:`plastid_default_scatter`). We highly recommend
        setting `rasterized` to `True`!

    kdalpha : float, optional
        Alpha level (transparency) for marginal distributions (Default: 0.7)

    bw_method : str
        Bandwith estimation method. See documentation for
        :obj:`scipy.stats.gaussian_kde`. (Default: "scott")

    Returns
    -------
    :class:`matplotlib.figure.Figure`
        Figure

    dict
        Dictionary of axes. `'orig'` refers to `ax`. The central panel is `'main'`.
        Other panels will be mapped to `'top'`, `'left`' et c, if they are created.
    """
    fig, axes = _scatterhist_help(axes=axes,top_height=top_height,right_width=right_width)
    xlog = ylog = False

    if color is None:
        color = next(get_color_cycle(axes["main"]))
    
    if mask_invalid == True:
        x, y = clean_invalid(x,y,min_x=min_x,max_x=max_x,min_y=min_y,max_y=max_y)

    if label is not None:
        scargs = copy.deepcopy(scargs)
        scargs["label"] = label

    if "x" in log:
        axes["main"].semilogx()
        axes["top"].semilogx()
        xlog = True
        xmask = x > 0
    else:
        xmask = numpy.tile(True,x.shape)

    if "y" in log:
        axes["main"].semilogy()
        axes["right"].semilogy()
        ylog = True
        ymask = y > 0
    else:
        ymask = numpy.tile(True,y.shape)
        
    axes["main"].scatter(x,y,edgecolor=color,**scargs)

    # kernel densities
    kargs = { "color" : color,
              "alpha" : kdalpha,
              "bw_method" : bw_method,
             }
    
    if ymask.sum() > 0:
        kde_plot(y[ymask],log=ylog,axes=axes["right"],vert=True,**kargs)

    if xmask.sum() > 0:
        kde_plot(x[xmask],log=xlog,axes=axes["top"],**kargs)

    return fig, axes



#==============================================================================
# Plots specific for genomics
#==============================================================================

def ma_plot(x,y,axes=None,color=None,label=None,xlabel=None,ylabel=None,title=None,
            right_width=0.2,log="xy",
            min_x=-numpy.inf,max_x=numpy.inf,min_y=-numpy.inf,max_y=numpy.inf,
            scargs=plastid_default_scatter,mask_invalid=True,
            kdalpha=0.7):
    """Plot fold changes (:math:`\log_{2} (y/x)`) as a function of the mean of x and y (:math:`0.5*(x+y)`).

    Parameters
    ----------
    x, y : :class:`numpy.ndarray` or list
        Pair arrays or lists of corresponding numbers

    color : matplotlib colorspec, optional
        Color to use in plot

    label : str or None, optional
        If not `None`, a label for plotting

    xlabel : str or None, optional
        If not `None`, an x-axis label

    ylabel : str or None, optional
        If not `None`, a y-axis label
        
    right_width : float, optional
        fraction of `axes` width to use in right panel containing
        marginal distribution (Default: 0.2)

    axes : :class:`matplotlib.axes.Axes`, dict, or `None`, optional
        If a :class:`matplotlib.axes.Axes`, an axes in which to place plot.
        This axes will be split into relevant panels.

        If a dict, this is expected to be a group of pre-split axes.

        If `None`, a new figure is generated, and axes are split.
        (Default: `None`)

    mask_invalid : bool, optional
        If `True` mask out any `nan`s or `inf`s, as these mess up kernel density
        estimates and histograms in matplotlib, even if in masked arrays

    min_x, min_y, max_x, max_y : number, optional
        If supplied, set values below `min_x` to `min_x`, values larger
        than `max_x` to `max_x` and so for `min_y` and `max_y`

    log : str, "", "x", "xy", or "xy", optional
        Plot these axes on a log scale (Default: "xy")

    scargs : Keyword arguments, optional
        Arguments to pass to :func:`~matplotlib.pyplot.scatter`
        (Default: :obj:`plastid_default_scatter` ). Recommend: set `rasterized`
        to `True`

    kdalpha : float, optional
        Alpha level (transparency) for marginal distributions (Default: 0.7)


    Returns
    -------
    :class:`matplotlib.figure.Figure`
        Figure

    dict
        Dictionary of axes. `'orig'` refers to `ax`. The central panel is `'main'`.
        Other panels will be mapped to `'top'`, `'left'` et c, if they are created.
    """

    do_setup = axes is None
    logs = numpy.ma.masked_invalid(numpy.log2(y/x))
    imask = ~logs.mask

    ratio = y/x
    mean = 0.5*(x+y)

    fig, axdict = scatterhist_y(mean[imask],ratio[imask],axes=axes,
                                min_x=min_x,max_x=max_x,
                                min_y=min_y,max_y=max_y,
                                log=log,right_width=right_width,
                                color=color,mask_invalid=mask_invalid,
                                label=label,kdalpha=kdalpha)

    if do_setup == True:
        axdict["main"].axhline(1,color=process_black,zorder=-5,linewidth=1)
        axdict["right"].axhline(1,color=process_black,zorder=-5,linewidth=1)
        axdict["right"].xaxis.set_ticklabels([])

        if ylabel is None:
            ylabel = "log2 fold change"

        axdict["main"].set_ylabel(ylabel)
        if xlabel is not None:
            axdict["main"].set_xlabel(xlabel)

        if title is not None:
            plt.suptitle(title)

    return fig, axdict


def phase_plot(counts,labels=None,cmap=None,color=None,lighten_by=0.2,fig={},line={},bar={}):
    """Phasing plot for ribosome profiling

    Creates a two-panel plot:

      - the top panel is a line graph indicating the fraction of reads
        as a function of read length

      - the bottom panel is a stacked bar graph, showing the fraction
        of reads in each codon position for each read length, with
        codon position 2 stacked above position 1 stacked above position 0


    Parameters
    ----------
    counts : :class:`numpy.ndarray`
        Nx3 array of raw counts, where each row represents a read length,
        and each column a codon phase

    labels : list, optional
        Labels for each stack. If `None`, stacks are labeled sequentially by number.
        (Default: `None`)
    
    lighten_by : float, optional
        Amount by which to lighten sequential blocks in each stack. (Default: 0.10)
    
    cmap : :class:`matplotlib.colors.Colormap`, optional
        Colormap from which to generate bar colors. If supplied, will override
        any `color` attribute in `**kwargs`.

    color : matplotlib colorspec, or list of these
        Colors to use in plot. Overridden if `cmap` is supplied.

    fig : dict
        Keyword arguments to :func:`matplotlib.pylot.figure`
        
    line : dict
        Keyword arguments to pass to :func:`matplotlib.pyplot.plot`
        in top panel

    bar : dict
        Keyword arguments to pass to :func:`matplotlib.pyplot.bar`
        in bottom panel


    Returns
    -------
    :class:`matplotlib.figure.Figure`
        Figure

    tuple
        Tuple of :class:`matplotlib.axes.Axes`; the first corresponding 
        to the line graph (top panel), the second, the bar graph (bottom).
    """
    fig, (ax1,ax2) = plt.subplots(nrows=2,ncols=1,sharex=True,**fig)

    totals = counts.sum(1)
    phases = (counts.astype(float).T/totals).T
    
    stacked_bar(phases,axes=ax2,labels=labels,lighten_by=lighten_by,cmap=cmap,color=color,**bar)
    ax2.set_xlabel("Read length (nt)")
    ax2.set_ylabel("Fraction in each phase")
    x = numpy.arange(len(totals)) + 0.5

    if "color" not in line:
        line["color"] = process_black
    
    ax1.plot(x,(totals.astype(float))/totals.sum(),**line)

    ax1.set_ylabel("Fraction of reads")
    
    return fig, (ax1,ax2)


