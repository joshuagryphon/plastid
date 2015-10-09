#!/usr/bin/env python
"""Sommon common plots that are not directly implemented in :mod:`matplotlib`,
as well as some specific plots used by :mod:`plastid`. For example:

    :func:`stacked_bar`
        Create a stacked bar graph.

    :func:`triangle_plot`
        Plot data lying on the plane x + y + z = k (e.g. codon phasing)
        in a homogeneous 2D representation

"""
import numpy
import matplotlib.pyplot as plt
from plastid.plotting.colors import lighten, process_black


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
    elif "color" not in kwargs:
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

def triangle_plot(data,axes=None,fn="scatter",vertex_labels=None,grid=None,clip=True,**kwargs):
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

    # set up canvas
    artists = []
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
# 
#==============================================================================



def phase_plot(counts,labels=None,cmap=None,lighten_by=0.2,line={},bar={}):
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
        
    line : dict
        Other keyword arguments to pass to :func:`matplotlib.pyplot.plot`
        in top panel

    bar : dict
        Other keyword arguments to pass to :func:`matplotlib.pyplot.bar`
        in bottom panel


    Returns
    -------
    :class:`matplotlib.figure.Figure`
        Figure

    tuple
        Tuple of :class:`matplotlib.axes.Axes`; the first corresponding 
        to the line graph (top panel), the second, the bar graph (bottom).
    """
    fig, (ax1,ax2) = plt.subplots(nrows=2,ncols=1,sharex=True,squeeze=True)

    totals = counts.sum(1)
    phases = (counts.astype(float).T/totals).T
    
    stacked_bar(phases,axes=ax2,labels=labels,lighten_by=lighten_by,cmap=cmap,**bar)
    ax2.set_xlabel("Read length (nt)")
    ax2.set_ylabel("Fraction in each phase")
    x = numpy.arange(len(totals)) + 0.5

    if "color" not in line:
        line["color"] = process_black
    
    ax1.plot(x,(totals.astype(float))/totals.sum(),**line)

    ax1.set_ylabel("Fraction of reads")
    
    return fig, (ax1,ax2)
