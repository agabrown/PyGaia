__all__ = ['order_points_for_sky_plot', 'plot_coordinate_transformation_on_sky']

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np

from pygaia.astrometry.coordinates import CoordinateTransformation


def order_points_for_sky_plot(x, y):
    """
    This code takes care of ordering the points (x,y), calculated for a sky map parallel or merdian, such
    that the drawing code can start at one end of the curve and end at the other (so no artifacts due to
    connecting the disjoint ends occur).

    Parameters
    ----------

    x - Set of x coordinates
    y - Set of y coordinates

    Returns
    -------

    x, y: Order set of (x,y) points
    """
    xroll = np.roll(x, 1)
    yroll = np.roll(y, 1)
    distance = (xroll - x) ** 2 + (yroll - y) ** 2
    indexmax = np.argmax(distance)
    return np.roll(x, -indexmax), np.roll(y, -indexmax)


def plot_coordinate_transformation_on_sky(transformation, outfile=None, no_title=False, no_labels=False,
                                          return_plot_object=False, lc=plt.cm.get_cmap('tab10').colors[0],
                                          tc=plt.cm.get_cmap('tab10').colors[1],
                                          lonpos=True):
    """
    Produce a sky-plot in a given coordinate system with the meridians and paralles for another
    coordinate system overlayed. The coordinate systems are specified through the
    pygaia.coordinates.Transformations enum. For example for Transformations.GAL2ECL the sky plot will be
    in Ecliptic coordinates with the Galactic coordinate grid overlayed.

    Keywords
    --------

    transformation - The coordinate transformation for which to make the plot (e.g.,
                     Transformations.GAL2ECL).
    outfile        - Save plot to this output file (default is to plot on screen). Make sure an extension
                     (.pdf, .png, etc) is included.
    noTitle        - If true do not include the plot title.
    noLabels       - If true do not include plot labels.
    returnPlotObject - If true return the matplotlib object used for plotting. Further plot elements can
                       then be added.
    lc             - Colour for gridlines.
    tc             - Colour for text labels.
    lonpos         - If true use longitude labels between 0 and 360 degrees.
    """
    ct = CoordinateTransformation(transformation)

    parallels = np.arange(-80.0, 90.0, 10.0)
    meridians = np.arange(0.0, 375.0, 15.0)
    meridian_max = np.deg2rad(85.0)

    default_proj = ccrs.PlateCarree()
    addtolabel = 0
    if lonpos:
        addtolabel = 360

    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Mollweide())
    ax.invert_xaxis()

    for thetaDeg in parallels:
        phi = np.linspace(-np.pi, np.pi, 1001)
        theta = np.zeros_like(phi) + np.deg2rad(thetaDeg)
        phirot, thetarot = ct.transform_sky_coordinates(phi, theta)
        phirot[(phirot > np.pi)] = phirot[(phirot > np.pi)] - 2 * np.pi
        x, y = np.rad2deg(phirot), np.rad2deg(thetarot)

        indices = (phirot >= 0.0)
        xplot = x[indices]
        yplot = y[indices]
        if any(indices):
            xplot, yplot = order_points_for_sky_plot(xplot, yplot)
        ax.plot(xplot, yplot, '-', color=lc, transform=default_proj)

        indices = (phirot < 0.0)
        xplot = x[indices]
        yplot = y[indices]
        if any(indices):
            xplot, yplot = order_points_for_sky_plot(xplot, yplot)
        ax.plot(xplot, yplot, '-', color=lc, transform=default_proj)

    for phiDeg in meridians:
        theta = np.linspace(-meridian_max, meridian_max, 1001)
        phi = np.zeros_like(theta) + np.deg2rad(phiDeg)
        phirot, thetarot = ct.transform_sky_coordinates(phi, theta)
        phirot[(phirot > np.pi)] = phirot[(phirot > np.pi)] - 2 * np.pi
        x, y = np.rad2deg(phirot), np.rad2deg(thetarot)

        indices = (phirot >= 0.0)
        xplot = x[indices]
        yplot = y[indices]
        if any(indices):
            xplot, yplot = order_points_for_sky_plot(xplot, yplot)
        ax.plot(xplot, yplot, '-', color=lc, transform=default_proj)

        indices = (phirot < 0.0)
        xplot = x[indices]
        yplot = y[indices]
        if any(indices):
            xplot, yplot = order_points_for_sky_plot(xplot, yplot)
        ax.plot(xplot, yplot, '-', color=lc, transform=default_proj)

    if not no_title:
        plt.title("Sky projection in " + ct.transformationStrings[1] + " coordinates with the corresponding " +
                  ct.transformationStrings[0] + " grid overlayed")

    if not no_labels:
        for theta in np.arange(-60, 90, 30):
            phirot, thetarot = ct.transform_sky_coordinates(0.0, np.deg2rad(theta))
            x, y = (np.rad2deg(phirot), np.rad2deg(thetarot))
            ax.text(x, y, "${0}$".format(theta), fontsize=16, va='bottom', ha='center', color=tc,
                    transform=default_proj)
        for phi in np.arange(-150, 0, 30):
            phirot, thetarot = ct.transform_sky_coordinates(np.deg2rad(phi), 0.0)
            x, y = (np.rad2deg(phirot), np.rad2deg(thetarot))
            ax.text(x, y, "${0}$".format(phi + addtolabel), fontsize=16, va='bottom', ha='center', color=tc,
                    transform=default_proj)
        for phi in np.arange(30, 210, 30):
            phirot, thetarot = ct.transform_sky_coordinates(np.deg2rad(phi), 0.0)
            x, y = (np.rad2deg(phirot), np.rad2deg(thetarot))
            ax.text(x, y, "${0}$".format(phi), fontsize=16, va='bottom', ha='center', color=tc,
                    transform=default_proj)

    if outfile is not None:
        plt.savefig(outfile)
    elif return_plot_object:
        return plt.gca()
    else:
        plt.show()
