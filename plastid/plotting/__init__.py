#!/usr/bin/env python
"""This subpackage contains plotting utilities to build some compound plots
not found in matplotlib. For demos, see :doc:`/examples/z_plotting`.


    ======================================    =========================================================
    **Submodule**                             **Contents**
    --------------------------------------    ---------------------------------------------------------
    :mod:`plastid.plotting.colors`            Utilities for manipulating, converting, darkening,
                                              and lightening colors

    :mod:`plastid.plotting.plots`             Plots, such as heatmaps with summary profiles,
                                              scatter plots with kernel density estimates of 
                                              marginal distributions, MA plots, et c

    :mod:`plastid.plotting.plotutils`         Utility functions for preprocessing data or manipulating
                                              axes
    ======================================    =========================================================
"""
from plastid.plotting.plots import scatterhist_x, scatterhist_y, scatterhist_xy, \
                                   ma_plot, phase_plot, \
                                   profile_heatmap, \
                                   triangle_plot, \
                                   stacked_bar, \
                                   kde_plot
