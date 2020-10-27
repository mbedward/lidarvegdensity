# NEWS lidarvegdensity

Version 0.1.8 2020-10-29

  * Fixed problem where matrix of ground elevations was being read incorrectly due
    to different X-Y index order of result arrays and the elevation raster. The
    conventional X-Y index order used by the raster package is now applied to all
    inputs and outputs.

  * Function get_prop_reflected now checks for scan points outside the vertical
    limits of the null reference model. Also added arguments to specify whether
    the function should ignore scan points outside the horizontal or vertical
    range, or stop with an error. Default behaviour is to ignore points outside
    the horizontal range and stop with an error for points outside the vertical
    range.

  * Function get_prop_reflected now accepts either a RasterLayer object or a
    matrix (in raster cell order) for the ground argument.

  * The named list returned by function get_max_returns now has a class attribute
    'TLSNullModel'.

  * The named list returned by function get_prop_reflected now has a class attribute
    'TLSResult'.

  * Added function res_to_brick to create a RasterBrick object from one of the
    arrays in the result object returned by get_prop_reflected. By default, the
    function creates a brick for the 'preflected' array. A brick object is useful
    for plotting and GIS-style calculations.
