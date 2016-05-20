#!/usr/bin/env python
"""This module contains utility functions for plotting:

    :func:`get_fig_axes`
        Helper function to check arguments of plotting functions. Retrieve figure
        and axes from `axes`. If `axes` is None, create figure and axes, and
        return those.

    :func:`split_axes`
        Split a :class:`matplotlib.axes.Axes` into one or more panels, with tied
        x and y axes. Also hides overlapping tick labels

    :func:`clean_invalid`
        Remove pairs of values from two arrays of data if either  element in the
        pair is `nan` or `inf`. Used to prepare data for histogram or violin
        plots; as even masked `nan` and `inf` values are not handled by these
        functions.

    :func:`get_kde`
        Evaluate a kernel density estimate over observed data in linear or
        log-transformed space (e.g. for making violin plots in log space,
        but having kernels appropriately scaled).
"""
import numpy
import scipy.stats
import matplotlib
import matplotlib.pyplot as plt


def get_fig_axes(axes=None):
    """Retrieve figure and axes from `axes`. If `axes` is None, both.
    
    Used as a helper function for replotting atop existing axes, by functions
    defined in :mod:`plastid.plotting.plots`.


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
    """Split the spaces taken by one axes into one or more panes, setting the original axes invisible.

    Parameters
    ----------
    ax : :class:`matplotlib.axes.Axes`
        Axes to split

    top_height, left_width, right_width, bottom_height : float, optional
        If not `None`, a panel on the corresponding side of the `ax` will
        be created, using whatever fraction is specified (e.g. 0.1 to use
        10% of total height).

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
        axes["main"].xaxis.tick_bottom()
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
        axes["main"].yaxis.tick_left()
        axes["right"].xaxis.get_ticklabels()[0].set_visible(False)
        

    return axes

def clean_invalid(x,y,min_x=-numpy.inf,min_y=-numpy.inf,max_x=numpy.inf,max_y=numpy.inf):
    """Remove corresponding values from x and y when one or both of those is `nan` or `inf`,
    and optionally truncate values to minima and maxima

    Parameters
    ----------
    x, y : :class:`numpy.ndarray` or list
        Pair arrays or lists of corresponding numbers

    min_x, min_y, max_x, max_y : number, optional
        If supplied, set values below `min_x` to `min_x`, values larger
        than `max_x` to `max_x` and so for `min_y` and `max_y`

    Returns
    -------
    :class:`numpy.ndarray`
        A shortened version of `x`, excluding invalid values

    :class:`numpy.ndarray`
        A shortened version of `y`, excluding invalid values
    """
    x = numpy.array(x).astype(float)
    y = numpy.array(y).astype(float)

    x[x < min_x] = min_x
    x[x > max_x] = max_x
    y[y < min_y] = min_y
    y[y > max_y] = max_y
    
    newmask = numpy.isinf(x) | numpy.isnan(x) | numpy.isinf(y) | numpy.isnan(y) 
    x = x[~newmask]
    y = y[~newmask]


    return x,y 

def get_kde(data,log=False,base=2,points=100,bw_method="scott"):
    """Estimate a kernel density (kde) over `data`

    Parameters
    ----------
    data : :class:`numpy.ndarray`
        Data to build kde over

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
    :class:`numpy.ndarray`
        Points over which kde is evaluated (x-values), in non-log space

    :class:`numpy.ndarray`
        Value of kde (y-values), in non-log space
    """
    if log == True:
        if base == 2:
            func = numpy.log2
        elif base == 10:
            func = numpy.log10
        elif base == numpy.e:
            func = numpy.log
        else:
            raise ValueError("kde: Base must be 2, 10, or numpy.e")

        data = func(data)
        domain = func(numpy.logspace(data.min(),data.max(),base=base,num=points))
    else:
        domain = numpy.linspace(data.min(),data.max(),points)

    kde = scipy.stats.gaussian_kde(data,bw_method=bw_method)

    curve = kde.evaluate(domain)

    if log == True:
        domain = base**domain

    return domain, curve

