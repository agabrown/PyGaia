Plotting utilities
==================

Provides plotting utilities built on top of `Matplotlib <https://matplotlib.org/>`_.

Coordinate transformations
--------------------------

Detailed API: :py:mod:`pygaia.plot.sky`

This package provides a function to visualize the coordinate transformations defined in :py:mod:`pygaia.astrometry.coordinates` on the sky. The plots reproduce figures 3.1.2 to 3.1.7 in the `Hipparcos and Tycho Catalogues Volume 1 <https://www.cosmos.esa.int/documents/532822/552851/vol1_all.pdf/99adf6e3-6893-4824-8fc2-8d3c9cbba2b5>`_.

An example plot:

.. plot::
    
    import matplotlib.pyplot as plt
    from pygaia.plot.sky import plot_coordinate_transformation_on_sky
    from pygaia.astrometry.coordinates import Transformations

    fig = plt.figure(figsize=(10,5))
    plot_coordinate_transformation_on_sky(Transformations.GAL2ICRS, fig)
    plt.show()
