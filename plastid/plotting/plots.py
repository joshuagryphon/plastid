#!/usr/bin/env python
"""Sommon common plots that are not directly implemented in :mod:`matplotlib`,
as well as some specific plots used by :mod:`plastid`. For example:

    :func:`stacked_bar`
        Create a stacked bar graph

    :func:`scatterhist_x`, :func:`scatterhist_y`, and :func:`scatterhist_xy`
        Create scatter plots with kernel density estimates of the marginal 
        distributions along side them

    :func:`profile_heatmap`
        Plot a heatmap, with a columnwise median (or other profile summarizing
        the heatmap) in an upper panel above it, with aligned coordinates

    :func:`triangle_plot`
        Plot data lying on the plane x + y + z = k (e.g. codon phasing)
        in a homogeneous 2D representation

In addition, there are several utility functions:

    :func:`split_axes`
        Split a :class:`matplotlib.axes.Axes` into one or more panels,
        with tied x and y axes. Also hides overlapping tick labels

    :func:`remove_invalid`
        Remove pairs of values from two arrays of data if either 
        element in the pair is `nan` or `inf`. Used to prepare data
        for histogram or violin plots; as even masked `nan`s and `inf`s
        are not handled by these functions.

"""
import copy
import numpy
import matplotlib
import matplotlib.pyplot as plt
from plastid.plotting.colors import lighten, darken, process_black


def _plot_helper(axes):
    """Helper function - performs setup for plotting funcs below.

    Parameters
    ----------
    axes : :class:`matplotlib.axes.Axes` or `None`
        Axes in which to place plot. If `None`, a new figure is generated.

    Returns
    -------
    :class:`matplotlib.figure.Figure`
        Parent figure of axes
    
    :class:`matplotlib.axes.Axes`
        Axes containing plot
    """
    if axes is None:
        fig = plt.figure()
        ax = plt.gca()
    else:
        ax = axes
        fig = ax.figure

    return fig, ax

def split_axes(ax,top_height=0,left_width=0,right_width=0,bottom_height=0,main_ax_kwargs={},
                other_ax_kwargs={}):
    """Split the spaces taken by one axes into one or more panes. The original axes is made invisible.

    Parameters
    ----------
    ax : :class:`matplotlib.axes.Axes`
        Axes to split

    top_height, left_width, right_width, bottom_height : float, optional
        If not `None`, a panel on the corresponding side of the `ax` will be created, using
        whatever fraction is specified (e.g. 0.1 -> 10% of total height).

    main_ax_kwargs : dict
        Dictionary of keyword arguments for central panes, passed
        to :meth:`matplotlib.figure.Figure.add_axes`

    other_ax_kwargs : dict
        Dictionary of keyword arguments for peripheral panes, passed
        to :meth:`matplotlib.figure.Figure.add_axes`


    Returns
    -------
    dict
        Dictionary of axes. `'orig'` refers to `ax`. The central panel is `'main'`.
        Other panels will be mapped to `'top'`, `'left`' et c, if they are created.
    """
    fig = ax.figure
    ax.set_visible(False)
    rects = {}
    axes = { "orig" : ax }
    mplrc = matplotlib.rcParams

    buf_left  = mplrc["figure.subplot.left"]
    buf_bot   = mplrc["figure.subplot.bottom"]
    buf_right = 1.0 - mplrc["figure.subplot.right"] 
    buf_top   = 1.0 - mplrc["figure.subplot.top"]

    hscale = 1.0 - buf_top - buf_bot
    wscale = 1.0 - buf_left - buf_right

    main_height = (1.0 - bottom_height - top_height)*hscale
    main_width  = (1.0 - left_width - right_width)*wscale

    bottom_height *= hscale
    top_height    *= hscale
    left_width    *= wscale
    right_width   *= wscale

    main_left  = buf_left + left_width
    main_right = main_left + main_width
    main_bot   = buf_bot + bottom_height
    main_top   = main_bot + main_height


    # each rect is left,bottom,width,height
    rects = {}
    rects["main"] = [main_left,
                     main_bot,
                     main_width,
                     main_height]
    if left_width > 0:
        rects["left"] = [buf_left,
                         main_bot,
                         left_width,
                         main_height]
    if right_width > 0:
        rects["right"] = [main_right,
                          main_bot,
                          right_width,
                          main_height]
    if bottom_height > 0:
        rects["bottom"] = [main_left,
                           main_bot,
                           main_width,
                           bottom_height]
    if top_height > 0:
        rects["top"] = [main_left,
                        main_top,
                        main_width,
                        top_height]

    axes["main"] = fig.add_axes(rects["main"],zorder=100,**main_ax_kwargs)

    for axes_name, rect in rects.items():
        if axes_name == "main":
            pass
        else:
            ax_kwargs = other_ax_kwargs
            if axes_name in ("right","left"):
                ax_kwargs["sharey"] = axes["main"]
                if "sharex" in ax_kwargs:
                    ax_kwargs.pop("sharex")
            if axes_name in ("top","bottom"):
                ax_kwargs["sharex"] = axes["main"]
                if "sharey" in ax_kwargs:
                    ax_kwargs.pop("sharey")

            axes[axes_name] = fig.add_axes(rect,zorder=50,**ax_kwargs)

    if "top" in axes:
        axes["top"].xaxis.tick_top()
        # prevent tick collisions
        axes["top"].yaxis.get_ticklabels()[0].set_visible(False)

    if "bottom" in axes:
        axes["bottom"].xaxis.tick_bottom()
        for t in axes["main"].xaxis.get_ticklabels():
            t.set_visible(False)

    if "left" in axes:
        axes["left"].yaxis.tick_left()
        axes["left"].xaxis.get_ticklabels()[-1].set_visible(False)
        for t in axes["main"].yaxis.get_ticklabels():
            t.set_visible(False)

    if "right" in axes:
        axes["right"].yaxis.tick_right()
        axes["right"].xaxis.get_ticklabels()[0].set_visible(False)
        

    return axes

def remove_invalid(x,y,min_x=-numpy.inf,min_y=-numpy.inf,max_x=numpy.inf,max_y=numpy.inf):
    """Remove corresponding values from x and y when one or both of those is nan or inf,
    and optionally truncate values to minima and maxima

    Parameters
    ----------
    x, y : :class:`numpy.ndarray` or list
        Pair arrays or lists of corresponding numbers

    min_x, min_y, max_x, max_y : number, optional
        If supplied, set values below `min_x` to `min_x`, values larger
        than `max_x` to `max_x` and so for `min_y` and `max-y`

    Returns
    -------
    :class:`numpy.ndarray
        A shortened version of `x`, excluding invalid values

    :class:`numpy.ndarray
        A shortened version of `y`, excluding invalid values
    """
    x = numpy.array(x).astype(float)
    y = numpy.array(y).astype(float)
    newmask = numpy.isinf(x) | numpy.isnan(x) | numpy.isinf(y) | numpy.isnan(y) 
    x = x[~newmask]
    y = y[~newmask]
    x[x < min_x] = min_x
    x[x > max_x] = max_x
    y[y < min_y] = min_y
    y[y > max_y] = max_y

    return x,y 


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
        any `color` attribute in `**kwargs`.
        
    **kwargs : keyword arguments
        Other keyword arguments to pass to :func:`matplotlib.pyplot.bar`

    Returns
    -------
    :class:`matplotlib.figure.Figure`
        Parent figure of axes
    
    :class:`matplotlib.axes.Axes`
        Axes containing plot
    """
    fig, ax = _plot_helper(axes)
    rows, cols = data.shape
    labels = labels if labels is not None else range(rows)
    defaults = [("align","center"),
                ("width",0.8)]

    if cmap is not None:
        kwargs["color"] = cmap(numpy.linspace(0,1.0,num=10))
    elif kwargs.get("color",None) is None:
        kwargs["color"] = [next(ax._get_lines.color_cycle) for _ in range(rows)]
        
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
    """Convert points from triangular space to Cartesian space for plotting

    Parameters
    ----------
    data : class:`numpy.ndarray`
        Nx2 or Nx3 list or array of points in triangular space, where
        the first column is the first coordinate, the second column
        the second, and the third, if supplied, the third
    
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
        Array of data, in which each row is a stack, each column a value in that stack.

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
    fig, ax = _plot_helper(axes)
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

def profile_heatmap(data,profile=None,x=None,axes=None,sort_fn=sort_max_position,
                    im_args={},plot_args={}):
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
        Sort rows in `data` by this function before plotting
        (Default: sort by ascending argmax of each row)

    im_args : dict
        Keyword arguments to pass to :func:`matplotlib.pyplot.imshow`

    plot_args : dict
        Keyword arugments to pass to :func:`matplotlib.pyplot.plot`
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
    fig, ax = _plot_helper(axes)
    axes = split_axes(ax,top_height=0.2)

    if sort_fn is None:
        sort_indices = numpy.arange(data.shape[0])
    else:
        sort_indices = sort_fn(data)

    if x is None:
        x = numpy.arange(0,data.shape[1])

    if profile is None:
        profile = numpy.nanmedian(data,axes=0)
    
    plot_kwargs = copy.deepcopy(plot_args)
    im_args     = copy.deepcopy(im_args)
    im_args["aspect"] = "auto"
    im_args["extent"] = [x.min(),x.max(),0,data.shape[0]]#,0]
    im_args["origin"] = "upper"

    axes["top"].plot(x,profile,**plot_args)
    axes["top"].set_ylim(0,profile.max())
    axes["top"].set_xlim(x.min(),x.max())
    axes["top"].set_yticks([])
    axes["top"].xaxis.tick_bottom()
    axes["top"].set_frame_on(False)

    axes["main"].xaxis.tick_bottom()
    axes["main"].imshow(data[sort_indices,:],
                        vmin=numpy.nanmin(data),
                        vmax=numpy.nanmax(data),
                        **im_args)

    return fig, axes


#==============================================================================
# Scatter plots with marginal distributions
#==============================================================================


plastid_default_scatter = {
    "marker"    : "o",
    "alpha"     : 0.7,
    "facecolor" : "none",
    "s"         : 8,
}
"""Default parameters for scatter plots"""


plastid_default_marginalplot = {
    "showextrema" : False,
}


def _scatterhist_help(axes=None,
                      top_height=0,left_width=0,right_width=0,bottom_height=0,
                      ):
    """Create a scatter plot with the marginal distribution for `x`
    plotted in a separate pain as a kernel density estimate.

    Parameters
    ----------
    x : list or :class:`numpy.ndarray`
        x values

    y : list or :class:`numpy.ndarray`
        y values

    color : matplotlib colorspec, or None, optional
        Color to plot data series

    axes : :class:`matplotlib.axes.Axes`, dict, or `None`, optional
        If a :class:`matplotlib.axes.Axes`, an axes in which to place plot.
        This axes will be split into relevant panels.

        If a dict, this is expected to be a group of pre-split axes.

        If `None`, a new figure is generated, and axes are split.
        (Default: `None`)
       
    scargs : dict
        Keyword arguments to pass to :func:`~matplotlib.pyplot.scatter`
        Default: :obj:`plastid_default_scatter`)

    args : dict
        Keyword arguments to pass to :func:`~matplotlib.pyplot.violinplot`,
        which draws the marginal distributions


    Returns
    -------
    :class:`matplotlib.figure.Figure`
        Parent figure of axes
    
    dict of :class:`matplotlib.axes.Axes`
        axes containing plot space
    """
    if axes is None:
        do_setup = True
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


def scatterhist_x(x,y,color=None,axes=None,
                  top_height=0.2,mask_invalid=True,
                  log="",
                  min_x=-numpy.inf,min_y=-numpy.inf,max_x=numpy.inf,max_y=numpy.inf,
                  scargs=plastid_default_scatter,
                  vargs=plastid_default_marginalplot,valpha=0.7):
    """Produce a scatter plot with a kernel density estimate of the marginal `x` distribution

    Parameters
    ----------
    x, y : :class:`numpy.ndarray` or list
        Pair arrays or lists of corresponding numbers

    color : matplotlib colorspec, optional
        Color to use in plot

    top_height : float, optional
        fraction of `axes` height to use in top panel containing
        marginal distribution (Default: 0.2)

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
        (Default: :obj:`plastid_default_scatter`)

    vargs : Keyword arguments, optional
        Arguments to pass to :func:`~matplotlib.pyplot.violinplot`, which draws
        the marginal distributions (Default : :obj:`plastid_default_marginalplot`)

    valpha : float, optional
        Alpha level (transparency) for marginal distributions (Default: 0.7)


    Returns
    -------
    :class:`matplotlib.figure.Figure`
        Figure

    dict
        Dictionary of axes. `'orig'` refers to `ax`. The central panel is `'main'`.
        Other panels will be mapped to `'top'`, `'left`' et c, if they are created.
    """
    fig, axes = _scatterhist_help(axes=axes,top_height=top_height)

    if "x" in log:
        axes["main"].semilogx()
        axes["top"].semilogx()

    if "y" in log:
        axes["main"].semilogy()

    if color is None:
        color = next(axes["main"]._get_lines.color_cycle)
    
    if mask_invalid == True:
       x, y = remove_invalid(x,y,min_x=min_x,max_x=max_x,min_y=min_y,max_y=max_y)

    axes["main"].scatter(x,y,edgecolor=color,**scargs)
    vdict = axes["top"].violinplot(x,positions=[0],vert=False,**vargs)

    violins = vdict["bodies"][0]
    violins.set_facecolor(color)
    violins.set_edgecolor(darken(color))
    violins.set_alpha(valpha)
    axes["top"].set_ylim(0.0,0.5)
    axes["top"].yaxis.set_ticklabels([])

    return fig, axes

def scatterhist_y(x,y,color=None,axes=None,
                  right_width=0.2,mask_invalid=True,log="xy",
                  min_x=-numpy.inf,min_y=-numpy.inf,max_x=numpy.inf,max_y=numpy.inf,
                  scargs=plastid_default_scatter,
                  vargs=plastid_default_marginalplot,valpha=0.7):
    """Produce a scatter plot with a kernel density estimate of the marginal `y` distribution

    Parameters
    ----------
    x, y : :class:`numpy.ndarray` or list
        Pair arrays or lists of corresponding numbers

    color : matplotlib colorspec, optional
        Color to use in plot

    right_width : float, optional
        fraction of `axes` width to use in right panel containing
        marginal distribution (Default: 0.2)

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
        (Default: :obj:`plastid_default_scatter`)

    vargs : Keyword arguments, optional
        Arguments to pass to :func:`~matplotlib.pyplot.violinplot`, which draws
        the marginal distributions (Default : :obj:`plastid_default_marginalplot`)

    valpha : float, optional
        Alpha level (transparency) for marginal distributions (Default: 0.7)


    Returns
    -------
    :class:`matplotlib.figure.Figure`
        Figure

    dict
        Dictionary of axes. `'orig'` refers to `ax`. The central panel is `'main'`.
        Other panels will be mapped to `'top'`, `'left`' et c, if they are created.
    """
    fig, axes = _scatterhist_help(axes=axes,right_width=right_width)

    if "x" in log:
        axes["main"].semilogx()

    if "y" in log:
        axes["main"].semilogy()
        axes["right"].semilogy()
        
    if color is None:
        color = next(axes["main"]._get_lines.color_cycle)
    
    if mask_invalid == True:
        x, y = remove_invalid(x,y,min_x=min_x,max_x=max_x,min_y=min_y,max_y=max_y)

    axes["main"].scatter(x,y,edgecolor=color,**scargs)
    vdict = axes["right"].violinplot(y,positions=[0],vert=True,**vargs)

    violins = vdict["bodies"][0]
    violins.set_facecolor(color)
    violins.set_edgecolor(darken(color))
    violins.set_alpha(valpha)
    axes["right"].set_xlim(0,0.5)


    return fig, axes

def scatterhist_xy(x,y,color=None,axes=None,
                   top_height=0.2,right_width=0.2,mask_invalid=True,log="xy",
                   min_x=-numpy.inf,min_y=-numpy.inf,max_x=numpy.inf,max_y=numpy.inf,
                   scargs=plastid_default_scatter,
                   vargs=plastid_default_marginalplot,valpha=0.7):
    """Produce a scatter plot with kernel density estimate of the marginal `x` and `y` distributions

    Parameters
    ----------
    x, y : :class:`numpy.ndarray` or list
        Pair arrays or lists of corresponding numbers

    color : matplotlib colorspec, optional
        Color to use in plot

    top_height : float, optional
        fraction of `axes` height to use in top panel containing
        marginal distribution (Default: 0.2)
        
    right_width : float, optional
        fraction of `axes` width to use in right panel containing
        marginal distribution (Default: 0.2)

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
        (Default: :obj:`plastid_default_scatter`)

    vargs : Keyword arguments, optional
        Arguments to pass to :func:`~matplotlib.pyplot.violinplot`, which draws
        the marginal distributions (Default : :obj:`plastid_default_marginalplot`)

    valpha : float, optional
        Alpha level (transparency) for marginal distributions (Default: 0.7)


    Returns
    -------
    :class:`matplotlib.figure.Figure`
        Figure

    dict
        Dictionary of axes. `'orig'` refers to `ax`. The central panel is `'main'`.
        Other panels will be mapped to `'top'`, `'left`' et c, if they are created.
    """
    fig, axes = _scatterhist_help(axes=axes,top_height=top_height,right_width=right_width)


    if "x" in log:
        axes["main"].semilogx()
        axes["top"].semilogx()

    if "y" in log:
        axes["main"].semilogy()
        axes["right"].semilogy()

    if color is None:
        color = next(axes["main"]._get_lines.color_cycle)
    
    if mask_invalid == True:
        x, y = remove_invalid(x,y,min_x=min_x,max_x=max_x,min_y=min_y,max_y=max_y)

    axes["main"].scatter(x,y,edgecolor=color,**scargs)
    vdictx = axes["right"].violinplot(y,positions=[0],vert=True,**vargs)
    vdicty = axes["top"].violinplot(x,positions=[0],vert=False,**vargs)

    vix = vdictx["bodies"][0]
    vix.set_facecolor(color)
    vix.set_edgecolor(darken(color))
    vix.set_alpha(valpha)

    viy = vdicty["bodies"][0]
    viy.set_facecolor(color)
    viy.set_edgecolor(darken(color))
    viy.set_alpha(valpha)

    axes["right"].set_xlim(0,0.5)
    axes["top"].set_ylim(0,0.5)


    return fig, axes




#==============================================================================
# Plots specific for genomics
#==============================================================================

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


