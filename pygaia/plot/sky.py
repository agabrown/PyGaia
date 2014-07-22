__all__ = ['plotCoordinateTransformationOnSky']

from pygaia.astrometry.coordinates import Transformations, CoordinateTransformation
from pygaia.utils import degreesToRadians, radiansToDegrees
from numpy import arange, argmax, roll, any, linspace, pi, zeros_like

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

def _orderGridlinePoints(x, y):
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
  xroll=roll(x,1)
  yroll=roll(y,1)
  distance=(xroll-x)**2+(yroll-y)**2
  indexmax=argmax(distance)
  return roll(x,-indexmax), roll(y,-indexmax)

def plotCoordinateTransformationOnSky(transformation, outfile=None, myProjection='hammer',
    noTitle=False, noLabels=False, returnPlotObject=False):
    """
    Produce a sky-plot in a given coordinate system with the meridians and paralles for another
    coordinate system overlayed. The coordinate systems are specified through the
    pygaia.coordinates.Transformations enum. For example for Transformations.GAL2ECL the sky plot will be
    in Ecliptic coordinates with the Galactic coordinate grid overlayed.

    Keywords
    --------

    transformation - The coordinate transformation for which to make the plot (e.g.,
                     Transformations.GAL2ECL)
    outfile        - Save plot to this output file (default is to plot on screen). Make sure an extension
                     (.pdf, .png, etc) is included.
    myProjection   - Use this map projection (default is 'hammer', see basemap documentation)
    noTitle        - If true do not include the plot title.
    noLabels       - If true do not include plot labels.
    returnPlotObject - If true return the matplotlib object used for plotting. Further plot elements can
                       then be added.
    """
    ct = CoordinateTransformation(transformation)

    parallels=arange(-80.0,90.0,10.0)
    meridians=arange(0.0,375.0,15.0)
    meridianMax=degreesToRadians(85.0)
    parallelsMax=degreesToRadians(179.0)

    fig=plt.figure(figsize=(12,6))
    basemapInstance=Basemap(projection=myProjection,lon_0=0, celestial=True)
    basemapInstance.drawmapboundary()

    for thetaDeg in parallels:
      phi=linspace(-pi,pi,1001)
      theta=zeros_like(phi)+degreesToRadians(thetaDeg)
      phirot, thetarot = ct.transformSkyCoordinates(phi, theta)
      x ,y = basemapInstance(radiansToDegrees(phirot), radiansToDegrees(thetarot))
      
      indices=(phirot>=0.0)
      xplot=x[indices]
      yplot=y[indices]
      if any(indices):
        xplot, yplot = _orderGridlinePoints(xplot, yplot)
      plt.plot(xplot,yplot,'b-')
     
      indices=(phirot<0.0)
      xplot=x[indices]
      yplot=y[indices]
      if any(indices):
        xplot, yplot = _orderGridlinePoints(xplot, yplot)
      plt.plot(xplot,yplot,'b-')

    for phiDeg in meridians:
      theta=linspace(-meridianMax,meridianMax,1001)
      phi=zeros_like(theta)+degreesToRadians(phiDeg)
      phirot, thetarot = ct.transformSkyCoordinates(phi, theta)
      x ,y = basemapInstance(radiansToDegrees(phirot), radiansToDegrees(thetarot))

      indices=(phirot>=0.0)
      xplot=x[indices]
      yplot=y[indices]
      if any(indices):
        xplot, yplot = _orderGridlinePoints(xplot, yplot)
      plt.plot(xplot,yplot,'b-')
     
      indices=(phirot<0.0)
      xplot=x[indices]
      yplot=y[indices]
      if any(indices):
        xplot, yplot = _orderGridlinePoints(xplot, yplot)
      plt.plot(xplot,yplot,'b-')

    if (not noTitle):
      plt.title("Sky projection in " + ct.transformationStrings[1] + " coordinates with the corresponding " + ct.transformationStrings[0] + " grid overlayed")

    if (not noLabels):
      for theta in arange(-60,90,30):
        phirot, thetarot=ct.transformSkyCoordinates(0.0,degreesToRadians(theta))
        x, y = basemapInstance(radiansToDegrees(phirot), radiansToDegrees(thetarot))
        plt.text(x,y,"${0}$".format(theta),fontsize=16,va='bottom',ha='center',color='r')
      for phi in arange(-180,0,30):
        phirot, thetarot=ct.transformSkyCoordinates(degreesToRadians(phi),0.0)
        x, y = basemapInstance(radiansToDegrees(phirot), radiansToDegrees(thetarot))
        plt.text(x,y,"${0}$".format(phi),fontsize=16,va='bottom',ha='center',color='r')
      for phi in arange(30,180,30):
        phirot, thetarot=ct.transformSkyCoordinates(degreesToRadians(phi),0.0)
        x, y = basemapInstance(radiansToDegrees(phirot), radiansToDegrees(thetarot))
        plt.text(x,y,"${0}$".format(phi),fontsize=16,va='bottom',ha='center',color='r')

    if (outfile != None):
      plt.savefig(outfile)
    elif (returnPlotObject):
      return plt.gca(), basemapInstance
    else:
      plt.show()
