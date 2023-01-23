"""
Function for visualizing coordinate transformations on the sky.
"""
import matplotlib.pyplot as plt
import numpy as np

from pygaia.astrometry.coordinates import CoordinateTransformation

__all__ = ["plot_coordinate_transformation_on_sky"]


def _order_points_for_sky_plot(x, y):
    """
    This code takes care of ordering the points (x,y), calculated for a sky map parallel
    or merdian, such that the drawing code can start at one end of the curve and end at
    the other (so no artifacts due to connecting the disjoint ends occur).

    Parameters
    ----------
    x : float array
        Set of x coordinates
    y : float array
        Set of y coordinates

    Returns
    -------
    x, y : float arrays
        Ordered set of (x,y) points
    """
    xroll = np.roll(x, 1)
    yroll = np.roll(y, 1)
    distance = (xroll - x) ** 2 + (yroll - y) ** 2
    indexmax = np.argmax(distance)
    return np.roll(x, -indexmax), np.roll(y, -indexmax)


def plot_coordinate_transformation_on_sky(
    transformation,
    fig,
    outfile=None,
    notitle=False,
    nolabels=False,
    lc=plt.cm.get_cmap("tab10").colors[0],
    tc=plt.cm.get_cmap("tab10").colors[1],
    lonpos=True,
    skyproj="hammer",
):
    """
    Produce a sky-plot in a given coordinate system with the meridians and parallels for
    another coordinate system overlayed. The coordinate systems are specified through
    the pygaia.coordinates.Transformations enum. For example for Transformations.GAL2ECL
    the sky plot will be in Ecliptic coordinates with the Galactic coordinate grid
    overlayed.

    Parameters
    --------
    transformation : pygaia.astrometry.coordinates.Transformations instance
        The coordinate transformation for which to make the plot (e.g.,
        Transformations.GAL2ECL).
    fig : matplotlib.figure.Figure
        Empty Figure instance in which to create the plot.
    outfile: str
        Save plot to this output file. Make sure an extension (.pdf, .png, etc) is included.
    notitle : boolean
        If true do not include the plot title.
    nolabels : boolean
        If true do not include plot labels.
    inax : matplotlib.axes.Axes
        If provided use the input Axes instance for plotting.
    lc : Matplotlib color specification
        Colour for gridlines.
    tc : Matplotlib color specification
        Colour for text labels.
    lonpos : boolean
        If true use longitude labels between 0 and 360 degrees.
    skyproj : str
        The geographic projection to use. Must be one of "aitoff", "hammer", or "mollweide". Default is "hammer".

    Returns
    -------

    Nothing.
    """
    if not skyproj.lower() in ["hammer", "mollweide", "aitoff"]:
        raise ValueError("Only Aitoff, Hammer, or Mollweide projections are supported")
    ct = CoordinateTransformation(transformation)

    parallels = np.arange(-80.0, 90.0, 10.0)
    meridians = np.arange(0.0, 375.0, 15.0)
    meridian_max = np.deg2rad(85.0)

    addtolabel = 0
    if lonpos:
        addtolabel = 360

    ax = fig.add_subplot(projection=skyproj)
    ax.set_longitude_grid(360)
    ax.set_latitude_grid(180)

    for thetaDeg in parallels:
        if np.mod(thetaDeg, 30) == 0:
            a = 1
        else:
            a = 0.3
        phi = np.linspace(-np.pi, np.pi, 1001)
        theta = np.zeros_like(phi) + np.deg2rad(thetaDeg)
        phirot, thetarot = ct.transform_sky_coordinates(phi, theta)
        phirot[(phirot > np.pi)] = phirot[(phirot > np.pi)] - 2 * np.pi

        indices = phirot >= 0.0
        xplot = phirot[indices]
        yplot = thetarot[indices]
        if any(indices):
            xplot, yplot = _order_points_for_sky_plot(xplot, yplot)
            ax.plot(-xplot, yplot, "-", color=lc, alpha=a)

        indices = phirot < 0.0
        xplot = phirot[indices]
        yplot = thetarot[indices]
        if any(indices):
            xplot, yplot = _order_points_for_sky_plot(xplot, yplot)
            ax.plot(-xplot, yplot, "-", color=lc, alpha=a)

    for phiDeg in meridians:
        if np.mod(phiDeg, 30) == 0:
            a = 1
        else:
            a = 0.3
        theta = np.linspace(-meridian_max, meridian_max, 1001)
        phi = np.zeros_like(theta) + np.deg2rad(phiDeg)
        phirot, thetarot = ct.transform_sky_coordinates(phi, theta)
        phirot[(phirot > np.pi)] = phirot[(phirot > np.pi)] - 2 * np.pi

        indices = phirot >= 0.0
        xplot = phirot[indices]
        yplot = thetarot[indices]
        if any(indices):
            xplot, yplot = _order_points_for_sky_plot(xplot, yplot)
            ax.plot(-xplot, yplot, "-", color=lc, alpha=a)

        indices = phirot < 0.0
        xplot = phirot[indices]
        yplot = thetarot[indices]
        if any(indices):
            xplot, yplot = _order_points_for_sky_plot(xplot, yplot)
            ax.plot(-xplot, yplot, "-", color=lc, alpha=a)

    if not notitle:
        plt.title(
            f"{skyproj.capitalize()} projection in "
            + ct.target_coordinates()
            + " coordinates with the corresponding "
            + ct.start_coordinates()
            + " grid overlayed"
        )

    if not nolabels:
        for theta in np.arange(-60, 90, 30):
            phirot, thetarot = ct.transform_sky_coordinates(0.0, np.deg2rad(theta))
            phirot[(phirot > np.pi)] = phirot[(phirot > np.pi)] - 2 * np.pi
            ax.text(
                -phirot,
                thetarot,
                "${0}$".format(theta),
                fontsize=16,
                va="bottom",
                ha="center",
                color=tc,
            )
        for phi in np.arange(-150, 0, 30):
            phirot, thetarot = ct.transform_sky_coordinates(np.deg2rad(phi), 0.0)
            phirot[(phirot > np.pi)] = phirot[(phirot > np.pi)] - 2 * np.pi
            ax.text(
                -phirot,
                thetarot,
                "${0}$".format(phi + addtolabel),
                fontsize=16,
                va="bottom",
                ha="center",
                color=tc,
            )
        for phi in np.arange(30, 210, 30):
            phirot, thetarot = ct.transform_sky_coordinates(np.deg2rad(phi), 0.0)
            phirot[(phirot > np.pi)] = phirot[(phirot > np.pi)] - 2 * np.pi
            ax.text(
                -phirot,
                thetarot,
                "${0}$".format(phi),
                fontsize=16,
                va="bottom",
                ha="center",
                color=tc,
            )

    if outfile is not None:
        plt.savefig(outfile)
