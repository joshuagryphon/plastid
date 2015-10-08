#!/usr/bin/env python
"""Sommon common plots that are not directly implemented in :mod:`matplotlib`,
as well as some specific plots used by :mod:`plastid`. For example:

    :func:`stacked_bar`
        Create a stacked bar graph.

"""
import numpy
import matplotlib.pyplot as plt
from plastid.plotting.colors import lighten


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


def stacked_bar(data,axes=None,labels=None,lighten_by=0.1,cmap=Spectral,**kwargs):
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
    
    cmap : colormap, optional
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
