XBeach Matlab Toolbox
=====================
    
The XBeach Matlab toolboxes accomodates several frequantly used pre-
and post-processing scripts for the XBeach model. Be aware that this
toolbox is not exhaustive. XBeach has many more functionalities than
exposed by this toolbox. Also note that this toolbox is *not* required
to run XBeach. You can run XBeach by providing the right input text
files in any way that suits your situation best

xb_analysis_
    Perform analysis on XBeach results
xb_gui_
    Simple gui to setup XBeach models (deprecated)
xb_io_
    Read and write XBeach files
xb_lib_
    General functions
xb_modelsetup_
    Facilitate model set-up, including grid optimization
xb_nesting_
    Routines to nest an XBeach model in another model (Delft3D or WW3)
xb_run_
    Routines to run XBeach models from Matlab, also support remote running on UNIX clusters
xb_testbed_
    Routines to facilitate XBeach testing in the testbed/skillbed
xb_visualise_
    Facilitate model output visualisation

.. _sec-xb_analysis:

xb_analysis
-----------

.. _xb_get_activeprofile:

xb_get_activeprofile
~~~~~~~~~~~~~~~~~~~~

.. code:: text

   XB_GET_ACTIVEPROFILE  Determines the range of the profile that is active
   
   Determines the x-range of the profile that is active by locating the
   area in which a certain percentage of the bed level change is located.
   By default, this percentage is 90%.
   
   Syntax:
   [x ix] = xb_get_activeprofile(xb, varargin)
   
   Input:
   xb        = XBeach output structure
   varargin  = dzfrac:     Fraction of total profile change that should be
                           included in the range returned
   
   Output:
   x         = Array with minimum and maximum x-locations with in between
               the requested percentage of bed level change
   xi        = Array with the corresponding indices in the x-grid
   
   Example:
   x = xb_get_activeprofile(xb)
   [x xi] = xb_get_activeprofile(xb)
   [x xi] = xb_get_activeprofile(xb, 'dzfrac', .8)
   
.. seealso:: xb_get_profile_


.. _xb_get_coastline:

xb_get_coastline
~~~~~~~~~~~~~~~~

.. code:: text

   XB_GET_COASTLINE  Determines coastline from 2D grid
   
   Determines coastline based on 2D grid by first determining the
   orientation of the grid and than finding the first grid cell that
   exceeds a certain elevation (default 0).
   
   TODO: add interpolation option
   
   varargout = xb_get_coastline(varargin)
   
   Syntax:
   
   Input:
   x           = x-coordinates of bathymetric grid
   y           = y-coordinates of bathymetric grid
   z           = elevations in bathymetric grid
   varargin    = level:        Level that needs to be exceeded
                 interpolate:  Boolean flag to determine whether result
                               should be interpolated to obtain a better
                               apporximation
   
   Output:
   xc          = x-coordinates of coastline
   yc          = y-coordinates of coastline
   dim         = cross-shore dimension (1 or 2)
   dir         = landward direction (negative or positive)
   idx         = logical matrix of the size of z indicating whether a cell
                 is on the coastline or not (xc = x(idx) and yc = y(idx)
                 if x and y are matrices of size of z)
   
   Example:
   xb_get_coastline
   
.. seealso:: xb_grid_orientation_


.. _xb_get_hydro:

xb_get_hydro
~~~~~~~~~~~~

.. code:: text

   XB_GET_HYDRO  Compute hydrodynamic parameters from XBeach output structure
   
   Compute hydrodynamic parameters like RMS wave heights over a
   cross-section split in low and high freqnecy waves. The same is done
   for orbital velocities and mean velocities. Also the water level setup
   is computed, if possible. The results are stored in an XBeach
   hydrodynamics structure and can be plotted with xb_plot_hydro.
   
   Syntax:
   xbo = xb_get_hydro(xb, varargin)
   
   Input:
   xb        = XBeach output structure
   varargin  = Trep:   repesentative wave period
   
   Output:
   xbo       = XBeach hydrodynamics structure
   
   Example:
   xbo = xb_get_hydro(xb)
   xbo = xb_get_hydro(xb, 'Trep', 12)
   
.. seealso:: xb_plot_hydro_, xb_get_morpho_, xb_get_spectrum_


.. _xb_get_morpho:

xb_get_morpho
~~~~~~~~~~~~~

.. code:: text

   XB_GET_MORPHO  Compute morphological parameters from XBeach output structure
   
   Compute morphological parameters like bed level change, erosion and
   sedimentation volumes and retreat distances from XBeach output
   structure. The results are stored in an XBeach morphology structure and
   can be plotted with xb_plot_morpho.
   
   Syntax:
   xbo = xb_get_morpho(xb, varargin)
   
   Input:
   xb        = XBeach output structure
   varargin  = level:  assumed storm surge level
   
   Output:
   xbo       = XBeach morphology structure
   
   Example:
   xbo = xb_get_morpho(xb)
   xbo = xb_get_morpho(xb, 'level', 0)
   
.. seealso:: xb_plot_morpho_, xb_get_hydro_, xb_get_spectrum_


.. _xb_get_profilespecs:

xb_get_profilespecs
~~~~~~~~~~~~~~~~~~~

.. code:: text

   XB_GET_PROFILESPECS  Determines a variety of profile characteristics from a dune profile
   
   More detailed description goes here.
   
   Syntax:
   varargout = xb_get_profilespecs(varargin)
   
   Input:
   varargin  =
   
   Output:
   varargout =
   
   Example:
   xb_get_profilespecs
   
   See also

.. _xb_get_sedbal:

xb_get_sedbal
~~~~~~~~~~~~~

.. code:: text

   XB_GET_SEDBAL  Computes sediment balance from XBeach output structure
   
   Computes total sedimentation, erosion and transports over domain
   borders and determines total sediment budget continuity.
   
   Syntax:
   xbo = xb_get_sedbal(xb, varargin)
   
   Input:
   xb        = XBeach output structure
   varargin  = t:          time at which balance should be computed
                           (approximately)
               margin:     grid margin
               porosity:   porosity of bed
               morfac:     morphological factor between transports and bed
   
   Output:
   xbo       = XBeach sediment balance structure
   
   Example:
   xbo = xb_get_sedbal(xb);
   xs_show(xbo);
   
.. seealso:: xb_read_output_


.. _xb_get_sedero:

xb_get_sedero
~~~~~~~~~~~~~

.. code:: text

   XB_GET_SEDERO  Compute sedimentation and erosion from profiles
   
   Compute sedimentation and erosion from profile development in time. A
   single x-axis and multiple z-axes are provided. Crossings with the
   initial profile are computed and areas of erosion and sedimentation
   distinguished. Based on a given surge level, the erosion volume and
   retreat distance above surge level are computed.
   
   Syntax:
   [sed ero dz R Q P] = xb_get_sedero(x,z)
   
   Input:
   x         = x-axis vector
   z         = z-axes matrix with time in first dimension
   varargin  = level:      maximum surge level
   
   Output:
   sed       = total sedimentation volume
   ero       = total erosion volume
   dz        = profile change
   R         = first profile crossing above surge level
   Q         = last profile crossing below surge level
   P         = one but last profile crossing below surge level
   
   Example:
   [sed ero] = xb_get_sedero(x, z, 'level, 5)
   
.. seealso:: xb_get_morpho_


.. _xb_get_sedtrans:

xb_get_sedtrans
~~~~~~~~~~~~~~~

.. code:: text

   XB_GET_SEDTRANS  Compute sediment transport parameters from XBeach output structure
   
   Compute sediment transport parameters like sediment concentrations and
   transport volumes from XBeach output structure. The results are stored
   in an XBeach sedimenttransport structure and can be plotted with
   xb_plot_sedtrans.
   
   Syntax:
   xbo = xb_get_sedtrans(xb, varargin)
   
   Input:
   xb        = XBeach output structure
   varargin  = Trep:   representative wave period
               rho:    sediment density
               por:    porosity
   
   Output:
   xbo       = XBeach sedimenttransport structure
   
   Example:
   xbo = xb_get_sedtrans(xb)
   
.. seealso:: xb_plot_sedtrans_, xb_get_hydro_, xb_get_morpho_, xb_get_spectrum_


.. _xb_get_spectrum:

xb_get_spectrum
~~~~~~~~~~~~~~~

.. code:: text

   XB_GET_SPECTRUM  Computes a spectrum from a timeseries
   
   Computes a spectrum from a timeseries. The result is stored in an
   XBeach spectrum structure and can be plotted using the xb_plot_spectrum
   function.
   
   FUNCTION IS AN ADAPTED VERSION OF R.T. MCCALL'S MAKESPECTRUM FUNCTION
               WITH MODIFICATIONS FROM HIS SPECTRUMSIMPLE FUNCTION
   
   Syntax:
   xbo = xb_get_spectrum(ts, varargin)
   
   Input:
   ts        = Timeseries in columns
   varargin  = sfreq:          sample frequency
               fsplit:         split frequency between high and low
                               frequency waves
               fcutoff:        cut-off frequency for high frequency waves
               detrend:        boolean to determine whether timeseries
                               should be linearly detrended before
                               computation
               filterlength:   smoothing window
   
   Output:
   xbo       = XBeach spectrum structure
   
   Example:
   xbo = xb_get_spectrum(ts)
   
.. seealso:: xb_plot_spectrum_, xb_get_hydro_, xb_get_morpho_


.. _xb_get_transect:

xb_get_transect
~~~~~~~~~~~~~~~

.. code:: text

   XB_GET_TRANSECT  Squeezes an XBeach output structure to a single transect
   
   Squeezes an XBeach output structure to a single transect
   
   Syntax:
   xb = xb_get_transect(xb, varargin)
   
   Input:
   xb        = XBeach output structure
   varargin  = transect:   transect number
               dim:        dimension that should be squeezed
   
   Output:
   xb        = squeezed XBeach output structure
   
   Example:
   xb = xb_get_transect(xb)
   xb = xb_get_transect(xb, 'transect', 10)

.. _xb_skill:

xb_skill
~~~~~~~~

.. code:: text

   XB_SKILL  Computes a variety skill scores
   
   Computes a variety skill scores: R^2, Sci, Relative bias, Brier Skill
   Score. Special feature: within the XBeach testbed, the results are
   stored to be able to show the development of the different skill scores
   in time.
   
   Syntax:
   [r2 sci relbias bss] = xb_skill(measured, computed, varargin)
   
   Input:
   measured  = Measured data where the first column contains independent
               values and the second column contains dependent values
   computed  = Computed data where the first column contains independent
               values and the second column contains dependent values
   initial   = Initial data where the first column contains independent 
               values and the second column contains dependent values
   varargin  = var:    Name of the variable that is supplied
   
   Output:
   r2        = R^2 skill score
   sci       = Sci skill score
   relbias   = Relative bias
   bss       = Brier Skill Score
   
   Example:
   [r2 sci relbias bss] = xb_skill(measured, computed)
   [r2 sci relbias bss] = xb_skill(measured, computed, 'var', 'zb')
   
.. seealso:: xb_plot_skill_


.. _sec-xb_gui:

xb_gui
------

.. _xb_gui:

xb_gui
~~~~~~

.. code:: text

   XB_GUI  One line description goes here.
   
   More detailed description goes here.
   
   Syntax:
   varargout = xb_gui(varargin)
   
   Input:
   varargin  =
   
   Output:
   varargout =
   
   Example:
   xb_gui
   
   See also 

.. _xb_gui_crop:

xb_gui_crop
~~~~~~~~~~~

.. code:: text

   XB_GUI_CROP  Rotate and crop bathymetry
   
   Aligns the coastline to the y-axis and enables the used to visually rop
   the bathymetry. The result is the rotated, cropped bathymetry.
   
   Syntax:
   varargout = xb_gui_crop(x,y,z,varargin)
   
   Input:
   x         = matrix with x-values
   y         = matrix with y-values
   z         = matrix with z-values
   varargin  = none
   
   Output:
   varargout = rotated and cropped versions of the x, y and z matrices
   
   Example:
   [xc yc zc] = xb_gui_crop(x,y,z)
   
.. seealso:: xb_gui_mergebathy_


.. _xb_gui_dragselect:

xb_gui_dragselect
~~~~~~~~~~~~~~~~~

.. code:: text

   XB_GUI_DRAGSELECT  Creates an drag and select interface
   
   Creates an drag and select interface for maps. Can fire a function
   after selection is complete. The function should take 5 arguments:
   default object and event arguments, the drag-and-select object, the x
   and y position of the selection.
   
   Syntax:
   xb_gui_dragselect(obj, varargin)
   
   Input:
   obj       = Object handle which should be drag-and-select enabled
   varargin  = cursor:     Enables crosshair cursor
               select:     Enables area selection
               fcn:        Function handle to be fired after selection
   
   Output:
   none
   
   Example:
   xb_gui_dragselect(axes, 'select', true, 'cursor', false, 'fcn', @drawrectangle)
   
.. seealso:: xb_gui_normconditions_, xb_gui_mergebathy_


.. _xb_gui_mergebathy:

xb_gui_mergebathy
~~~~~~~~~~~~~~~~~

.. code:: text

   XB_GUI_MERGEBATHY  Merge bathymetries
   
   Select JARKUS transects and Vaklodingen from a map, possible add some
   ArcGIS files and XBeach bathymetries to this selection and generate a
   merged bathymetry from these sources.
   
   Syntax:
   varargout = xb_gui_mergebathy
   
   Input:
   none
   
   Output:
   varargout = x:  x-coordinates
               y:  y-coordinates
               z:  z-coordinates
   
   Example:
   [x y z] = xb_gui_mergebathy
   
.. seealso:: xb_gui_normconditions_, xb_gui_dragselect_


.. _xb_gui_normconditions:

xb_gui_normconditions
~~~~~~~~~~~~~~~~~~~~~

.. code:: text

   XB_GUI_NORMCONDITIONS  Returns normative conditions based on location and frequency of exceedence
   
   Returns normative conditions along the Dutch coast based on location
   and frequency of exceedence according to probabilistic derivation of
   hydraulic boundary conditions for VTV2006.
   The location can be selected on a map or defined by a string that is
   interpreted by Google Maps.
   
   Syntax:
   varargout = xb_gui_modelsetup_hydro_norm
   
   Input:
   none
   
   Output:
   varargout = h:      Storm surge level
               Hs:     Significant wave height
               Tp:     Peak wave period
               x:      x-coordinate of picked location (RD)
               y:      y-coordinate of picked location (RD)
               freq:   Selected frequency of exceedance
   
   Example:
   h = xb_gui_modelsetup_hydro_norm
   [h Hs Tp] = xb_gui_modelsetup_hydro_norm
   [h Hs Tp x y freq] = xb_gui_modelsetup_hydro_norm
   
.. seealso:: xb_gui_dragselect_, xb_gui_mergebathy_


.. _sec-xb_io:

xb_io
-----

.. _xb_dat2nc:

xb_dat2nc
~~~~~~~~~

.. code:: text

   XB_DAT2NC  One line description goes here.
   
   More detailed description goes here.
   
   Syntax:
   varargout = xb_dat2nc(varargin)
   
   Input:
   varargin  =
   
   Output:
   varargout =
   
   Example:
   xb_dat2nc
   
   See also

.. _xb_get_dt:

xb_get_dt
~~~~~~~~~

.. code:: text

   XB_GET_DT  Read series of (average) dt values from XBlog.txt file
   
   Function to extract dt values from XBlog.txt file by means of regular
   expression.
   
   Syntax:
   dt = xb_average_dt(XBlog, varargin)
   
   Input:
   XBlog  = filename including path of XBlog.txt
   
   Output:
   dt = array of dt values
   
   Example:
   xb_average_dt
   
   See also

.. _xb_get_output:

xb_get_output
~~~~~~~~~~~~~

.. code:: text

   XB_GET_OUTPUT  Reads the spaceparams.tmpl file from the XBeach source code into a struct
   
   Reads the spaceparans.tmpl file from the XBeach source code into a
   struct. The file contains information on the possible output variables,
   their dimensions, name and description.
   
   Syntax:
   Input:
   fpath   = Path to spaceparams.tmpl
   
   Output:
   
   Output:
   output  = Structure array containing data from spaceparams.tmpl
   
   output = xb_get_output
   output = xb_get_output('spaceparams.tmpl')
   
   Example:
.. seealso:: xb_get_params_


.. _xb_get_params:

xb_get_params
~~~~~~~~~~~~~

.. code:: text

   XB_GET_PARAMS  Reads XBeach parameter types and defaults from XBeach source code
   
   Function to read XBeach params types and defaults from XBeach source
   code (params.F90). Useful to link to latest trunk update.
   
   Syntax:
   [params params_array] = xb_get_params(xbdir)
   
   Input:
   xbdir           = Directory in which XBeach source code can be found.
                     If not given, an attempt is made to use a default
                     path.
   
   Output:
   params          = structure array with listing of every parameter in
                     XBeach, including type, name, units, comment,
                     parameter type, default, minimum recommended and
                     maximum recommended values data.
   params_array    = array-version of params
   
   Example:
   [params params_array] = xb_get_params(xbdir)
   
.. seealso:: xb_read_params_, xb_write_params_


.. _xb_get_varinfo:

xb_get_varinfo
~~~~~~~~~~~~~~

.. code:: text

   XB_GET_VARINFO  One line description goes here.
   
   More detailed description goes here.
   
   Syntax:
   varargout = xb_get_varinfo(varargin)
   
   Input:
   varargin  =
   
   Output:
   varargout =
   
   Example:
   xb_get_varinfo
   
   See also

.. _xb_get_vars:

xb_get_vars
~~~~~~~~~~~

.. code:: text

   XB_GET_VARS  One line description goes here.
   
   More detailed description goes here.
   
   Syntax:
   varargout = xb_get_vars(varargin)
   
   Input:
   varargin  =
   
   Output:
   varargout =
   
   Example:
   xb_get_vars
   
   See also 

.. _xb_get_wavefiletype:

xb_get_wavefiletype
~~~~~~~~~~~~~~~~~~~

.. code:: text

   XB_GET_WAVEFILETYPE  Determines the type of wave definition file for XBeach input
   
   Analyzes the contents of a wave definition file for XBeach input and
   returns a string specifying the type of wave definition files.
   Currently, the following types can be returned: unknown, filelist,
   jonswap, jonswap_mtx, vardens, bcflist
   
   Syntax:
   type = xb_get_wavefiletype(filename)
   
   Input:
   filename  = filename of wave definition file to be analyzed
   
   Output:
   type      = string specifying the wave definition filetype
   types     = wave definition filetypes available
   counts    = matching scores of each filetype
   
   Example:
   type = xb_get_wavefiletype(filename)
   
.. seealso:: xb_read_params_, xb_read_waves_


.. _xb_read_bathy:

xb_read_bathy
~~~~~~~~~~~~~

.. code:: text

   XB_READ_BATHY  Read xbeach bathymetry files
   
   Routine to read xbeach bathymetry files.
   
   Syntax:
   xb = xb_read_bathy('xfile', <filename>, yfile, <filename>, depfile, <filename>, nefile, <filename>)
   
   Input:
   varargin    = xfile:    file name of x-coordinates file (cross-shore)
                 yfile:    file name of y-coordinates file (alongshore)
                 depfile:  file name of bathymetry file
                 ne_layer: file name of non erodible layer file
   
   Output:
   xb          = XBeach structure array
   
   Example:
   xb = xb_read_bathy('xfile', xfile, 'yfile', yfile)
   
.. seealso:: xb_write_bathy_, xb_read_input_


.. _xb_read_bcffile:

xb_read_bcffile
~~~~~~~~~~~~~~~

.. code:: text

   XB_READ_BCFFILE  Reads a bcf file generated by XBeach
   
   Reads a wave field realisation generated by XBeach.
   
   Syntax:
   data = xb_read_bcffile(filename, varargin)
   
   Input:
   filename  = Path to bcf file
   varargin  = none
   
   Output:
   data      = Matrix with wave field data
   
   Example:
   xb_read_bcffile('E001.bcf')
   xb_read_bcffile('Q001.bcf')
   
.. seealso:: xb_read_bcflist_


.. _xb_read_bcflist:

xb_read_bcflist
~~~~~~~~~~~~~~~

.. code:: text

   XB_READ_BCFLIST  Reads bcflist files generated by XBeach
   
   Reads ebcflist or qbcflist files generated by XBeach. The files contain
   references to other files containing realized wave and flux fields. The
   referred files are read as well. The result is retruned in the form of
   a XBeach structure.
   
   Syntax:
   xb = xb_read_bcflist(filename, varargin)
   
   Input:
   filename    = filename of the ebcflist or qbcflist file to be read
   varargin    = range:    unity-based numerical range of files to be read
                           (e.g. [4 5], 4 or [1 10])
   
   Output:
   xb          = XBeach structure array
   
   Example:
   xb = xb_read_bcflist(filename)
   
.. seealso:: xb_read_waves_


.. _xb_read_dat:

xb_read_dat
~~~~~~~~~~~

.. code:: text

   XB_READ_DAT  Reads DAT formatted output files from XBeach
   
   Reads DAT formatted output files from XBeach in the form of an XBeach
   structure. Specific variables can be requested in the varargin by means
   of an exact match, dos-like filtering or regular expressions (see
   strfilter)
   
   Syntax:
   xb = xb_read_dat(fname, varargin)
   
   Input:
   fname       = directory name that contains the dat files.
   varargin    = vars:     variable filters
                 start:    Start positions for reading in each dimension,
                           first item is zero
                 length:   Number of data items to be read in each
                           dimension, negative is unlimited
                 stride:   Stride to be used in each dimension
                 index:    Cell array with indices to read in each
                           dimension (overwrites start/length/stride)
                 dims:     Force the use of certain dimensions in
                           xb_dat_read. These dimensions are used for all
                           requested variables!
   
   Output:
   xb          = XBeach structure array
   
   Example:
   xb = xb_read_dat('.')
   xb = xb_read_dat('H.dat')
   xb = xb_read_dat('path_to_model/')
   xb = xb_read_dat('path_to_model/H.dat')
   xb = xb_read_dat('.', 'vars', 'H')
   xb = xb_read_dat('.', 'vars', 'H*')
   xb = xb_read_dat('.', 'vars', '/_mean$')
   xb = xb_read_dat('path_to_model/', 'vars', {'H', 'u*', '/_min$'})
   
.. seealso:: xb_read_output_, xb_read_netcdf_


.. _xb_read_dims:

xb_read_dims
~~~~~~~~~~~~

.. code:: text

   XB_READ_DIMS  read dimensions from xbeach output
   
   Routine to read the dimension from either netcdf of .dat xbeach output.
   The input argument "filename" can be the directory of the xbeach
   Syntax:
   XBdims   = xb_read_dims(varargin)
   
   Input:
   filename = file name. This can either be a output folder, a dims.dat file
              or a xboutput.nc file.
   
   Output:
   "filename" is a directory, it is assumed that the dimensions should be
   read from the "dims.dat" file inside the given directory.
   
   Output:
   XBdims   = structure containing the dimensions of xbeach output
              variables
   
   Example:
   xb_read_dims
   
   See also

.. _xb_read_input:

xb_read_input
~~~~~~~~~~~~~

.. code:: text

   XB_READ_INPUT  Read XBeach parameter file and all files referred in it
   
   Reads the XBeach settings from the params.txt file and all files that
   are mentioned in the settings, like grid and wave definition files. The
   settings are stored in a XBeach structure. The referred files are
   stored in a similar sub-structure.
   
   Syntax:
   xb = xb_read_input(filename)
   
   Input:
   filename   = params.txt file name
   varargin   = read_paths:        flag to determine whether relative
                                   paths should be read and included in
                                   the result structure
   
   Output:
   xb         = XBeach structure array
   
   Example:
   xb = xb_read_input(filename)
   
.. seealso:: xb_read_params_, xb_read_waves_


.. _xb_read_mpi_dims:

xb_read_mpi_dims
~~~~~~~~~~~~~~~~

.. code:: text

   XB_READ_MPI_DIMS  Reads the mpi dimensions from an XBlog file.
   
   Scans the XBlog file and reads mpi domain dimensions when specified.
   This function throws an exception when it was not possible to read the
   dimensions (due to incorrect input or the fact that there is no
   XBlog.txt file available, or the calculation was not run in mpi mode).
   
           For example:
                    0    1  107    1   63
                    1  106  106    1   63
                    2  210  106    1   63
                    3    1  107   62   62
                    4  106  106   62   62
                    5  210  106   62   62
   
   Syntax:
   dims = xb_read_mpi_dims(dr)
   
   Input:
   dr  = Directory where the XBlog file resides
   
   Output:
   dims = the n * 5 matrix included in the XBlog file that describes the
           mpi domain dimensions in which:
               First column:  domain number
               Second column: position of the left boundary in m direction (cross-shore)
               Third column:  Length of the domain in m direction (cross-shore)
               Fourth column: position of upper boundary in n direction (alongshore)
               Third column:  Length of the domain in n direction (alongshore)
   
   Example:
   dims = xb_read_mpi_dims('D:\testrun\');
   
.. seealso:: xb_plot_mpi_


.. _xb_read_netcdf:

xb_read_netcdf
~~~~~~~~~~~~~~

.. code:: text

   XB_READ_NETCDF  Reads NetCDF formatted output files from XBeach
   
   Reads NetCDF formatted output file from XBeach in the form of an
   XBeach structure. Specific variables can be requested in the varargin
   by means of an exact match, dos-like filtering or regular expressions
   (see strfilter)
   
   Syntax:
   variables = xb_read_netcdf(fname, varargin)
   
   Input:
   fname       = filename of the netcdf file
   varargin    = vars:     variable filters
                 start:    Start positions for reading in each dimension,
                           first item is zero
                 length:   Number of data items to be read in each
                           dimension, negative is unlimited
                 stride:   Stride to be used in each dimension
                 index:    Cell array with indices to read in each
                           dimension (overwrites start/length/stride)
   
   Output:
   variables   = XBeach structure array
   
   Example:
   xb = xb_read_netcdf('xboutput.nc')
   xb = xb_read_netcdf('xboutput.nc', 'vars', 'H')
   xb = xb_read_netcdf('xboutput.nc', 'vars', 'H*')
   xb = xb_read_netcdf('xboutput.nc', 'vars', '/_mean$')
   xb = xb_read_netcdf('path_to_model/xboutput.nc', 'vars', {'H', 'u*', '/_min$'})
   
.. seealso:: xb_read_output_, xb_read_dat_


.. _xb_read_output:

xb_read_output
~~~~~~~~~~~~~~

.. code:: text

   XB_READ_OUTPUT  Reads output files from XBeach
   
   Reads output files from XBeach. The actual work is done by either the
   xb_read_dat or xb_read_netcdf function. This function only determines
   which one to use. Specific variables can be requested in the varargin
   by means of an exact match, dos-like filtering or regular expressions
   (see strfilter)
   
   Syntax:
   varargout = xb_read_output(fname, varargin)
   
   Input:
   fname       = Path to the directory containing the dat files, a dat
                 file or the netcdf file to be read. This can also be a
                 XBeach run structure, which is translated to a path.
   varargin    = vars:         variable filters
   
   Output:
   varargout = XBeach structure array
   
   Example:
   xb = xb_read_output('path_to_model/')
   xb = xb_read_output('path_to_model/H.dat')
   xb = xb_read_output('path_to_model/', 'vars', {'H', 'u*', '/_mean$'})
   xb = xb_read_output('xboutput.nc')
   
.. seealso:: xb_read_input_, xb_write_input_, xb_read_dat_, xb_read_netcdf_


.. _xb_read_params:

xb_read_params
~~~~~~~~~~~~~~

.. code:: text

   XB_READ_PARAMS  read XBeach params.txt file
   
   Routine to read the xbeach settings from the params.txt file. The
   settings are stored in a XBeach structure.
   
   Syntax:
   xb = xb_read_params(filename)
   
   Input:
   filename   = params.txt file name
   varargin   = none
   
   Output:
   xb         = XBeach structure array
   
   Example:
   xb_read_params
   
.. seealso:: xb_read_input_, xb_read_waves_


.. _xb_read_parlist:

xb_read_parlist
~~~~~~~~~~~~~~~

.. code:: text

   XB_READ_PARLIST  Read list of available parameters from output file or directory
   
   Read list of available parameters from output file or directory.
   Returns a cell array with names.
   
   Syntax:
   parlist = xb_read_parlist(fname, varargin)
   
   Input:
   fname     = output file or directory
   varargin  = None
   
   Output:
   parlist   = cell array with parameter names
   
   Example:
   parlist = xb_read_parlist(pwd);
   
.. seealso:: xb_read_output_


.. _xb_read_ship:

xb_read_ship
~~~~~~~~~~~~

.. code:: text

   XB_READ_SHIP  One line description goes here.
   
   More detailed description goes here.
   
   Syntax:
   varargout = xb_read_ship(varargin)
   
   Input:
   varargin  =
   
   Output:
   varargout =
   
   Example:
   xb_read_ship
   
   See also

.. _xb_read_tide:

xb_read_tide
~~~~~~~~~~~~

.. code:: text

   XB_READ_TIDE  Reads tide definition file for XBeach input
   
   Reads a tide definition file containing a nx3 matrix of which the first
   column is the time definition and the second and third column the water
   level definition at respectively the seaward and landward boundary of
   the model.
   
   Syntax:
   xb  = xb_read_tide(filename)
   
   Input:
   filename    = filename of tide definition file
   varargin    = none
   
   Output:
   xb          = XBeach structure array
   
   Example:
   xb  = xb_read_tide(filename)
   
.. seealso:: xb_read_params_, xb_write_tide_


.. _xb_read_waves:

xb_read_waves
~~~~~~~~~~~~~

.. code:: text

   XB_READ_WAVES  Reads wave definition files for XBeach input
   
   Determines the type of wave definition file and reads it into a XBeach
   structure. If a filelist is given, also the underlying files are read
   and stored. The resulting struct can be inserted into the generic
   XBeach structure.
   
   Syntax:
   xb  = xb_read_waves(filename, varargin)
   
   Input:
   filename    = filename of wave definition file
   varargin    = none
   
   Output:
   xb          = XBeach structure array
   
   Example:
   xb  = xb_read_waves(filename)
   
.. seealso:: xb_read_params_, xb_write_waves_


.. _xb_write_bathy:

xb_write_bathy
~~~~~~~~~~~~~~

.. code:: text

   XB_WRITE_BATHY  Writes XBeach bathymetry files from XBeach structure
   
   Writes XBeach bathymetry files x, y, depth and non-erodable layers
   based on a XBeach structure.
   
   Syntax:
   [xfile yfile depfile ne_layer] = xb_write_bathy(xb, varargin)
   
   Input:
   xb          = XBeach structure array
   varargin    = path:         path to output directory
                 xfile:        filename of x definition file
                 yfile:        filename of y definition file
                 depfile:      filename of depth definition file
                 ne_layerfile: filename of non-erodable layer definition
                               file
   
   Output:
   varargout   = filenames of created definition files, if used
   
   Example:
   [xfile yfile depfile ne_layer] = xb_write_bathy(xb)
   
.. seealso:: xb_read_bathy_, xb_write_input_


.. _xb_write_input:

xb_write_input
~~~~~~~~~~~~~~

.. code:: text

   XB_WRITE_INPUT  Write XBeach params.txt file and all files referred in it
   
   Writes the XBeach settings from a XBeach structure in a parameter file.
   Also the files that are referred to in the parameter file are written,
   like grid and wave definition files.
   
   Syntax:
   xb_write_input(filename, xb, varargin)
   
   Input:
   filename  = filename of parameter file
   xb        = XBeach structure array
   varargin  = write_paths:  flag to determine whether definition files
                             should be written or just referred
               xbdir:  option to parse xbeach code directory (to read
                       parameter info)
   
   Output:
   none
   
   Example:
   xb_write_input(filename, xb)
   
.. seealso:: xb_read_input_, xb_write_params_


.. _xb_write_params:

xb_write_params
~~~~~~~~~~~~~~~

.. code:: text

   XB_WRITE_PARAMS  Write XBeach settings to params.txt file
   
   Routine to create a XBeach settings file. The settings in the XBeach
   structure are written to "filename". Optionally an alternative header
   line or directory containing params.f90 can be defined.
   
   Syntax:
   varargout = xb_write_params(filename, xb, varargin)
   
   Input:
   filename   = file name of params file
   xb         = XBeach structure array
   varargin   = header:    option to parse an alternative header string
                xbdir :    option to parse an alternative xbeach code directory
   
   Output:
   none
   
   Example:
   xb_write_params(filename, xb)
   
.. seealso:: xb_write_input_, xb_read_params_


.. _xb_write_plot:

xb_write_plot
~~~~~~~~~~~~~

.. code:: text

   SAVEPLOT: provide figure handle, output directory and filename (without extension)

.. _xb_write_ship:

xb_write_ship
~~~~~~~~~~~~~

.. code:: text

   XB_WRITE_SHIP  One line description goes here.
   
   More detailed description goes here.
   
   Syntax:
   varargout = xb_write_ship(varargin)
   
   Input:
   varargin  =
   
   Output:
   varargout =
   
   Example:
   xb_write_ship
   
   See also

.. _xb_write_tide:

xb_write_tide
~~~~~~~~~~~~~

.. code:: text

   XB_WRITE_TIDE  Writes tide definition file for XBeach input
   
   Writes a tide definition file containing a nx3 matrix of which the
   first column is the time definition and the second and third column the
   waterlevel definition at respectively the seaward and landward boundary
   of the model. Returns the filename of the tide file.
   
   Syntax:
   filename = xb_write_tide(xb)
   
   Input:
   xb          = XBeach structure array
   varargin    = path:     path to output directory
                 filename: filename of tide definition file
   
   Output:
   filename    = filename to be referred in parameter file
   
   Example:
   filename = xb_read_tide(xb)
   
.. seealso:: xb_read_params_, xb_read_tide_


.. _xb_write_waves:

xb_write_waves
~~~~~~~~~~~~~~

.. code:: text

   XB_WRITE_WAVES  Writes wave definition files for XBeach input
   
   Writes JONSWAP or variance density spectrum files for XBeach input. In
   case of conditions changing in time, a file list file is created
   refering to multiple wave definition files. In case of a JONSWAP
   spectrum, the file list file can be omitted and a single matrix
   formatted file is created. Returns the filename of the file to be
   referred in the params.txt file.
   
   In order to generate time varying wave conditions, simply add an extra
   dimension to the input arguments specifying the spectrum. The
   single-valued parameters Hm0, Tp, dir, gammajsp, s fnyq, duration and
   timestep then become one-dimensional. The one- and two-dimensional
   parameters freqs, dirs and vardens then become two- and
   three-dimensional respectively. It is not necessary to provide
   time-varying values for all parameters. In case a specific parameter is
   constant, simply provide the constant value. The value is reused in
   each period of time. However, it is not possible to provide for one
   parameter more than one value and for another too, while the number of
   values is not the same.
   
   Syntax:
   filename = xb_write_waves(xb, varargin)
   
   Input:
   xb          = XBeach structure array that overwrites the
                 default varargin options (optional)
   varargin    = path:             path to output directory
                 filelist_file:    name of filelist file without extension
                 jonswap_file:     name of jonswap file without extension
                 vardens_file:     name of vardens file without extension
                 unknown_file:     name of unknown wave file without
                                   extension
                 omit_filelist:    flag to omit filelist generation in
                                   case of jonswap spectrum
   
   Output:
   filename = filename to be referred in parameter file
   
   Example:
   filename = xb_write_waves(xb)
   
.. seealso:: xb_write_input_, xb_read_waves_


.. _sec-xb_io/xb_dat:

xb_io/xb_dat
------------

.. _xb_dat_dims:

xb_dat_dims
~~~~~~~~~~~

.. code:: text

   XB_DAT_DIMS  Returns the lengths of all dimensions of a XBeach DAT file
   
   Returns an array with the lengths of all dimensions of a XBeach DAT
   file. The functionality works similar to the Matlab size() function on
   variables.
   
   Syntax:
   dims = xb_dat_dims(filename, varargin)
   
   Input:
   filename    = Filename of DAT file
   varargin    = ftype:    datatype of DAT file (double/single)
   
   Output:
   dims        = Array with lengths of dimensions
   names       = Cell array with names of dimensions (x/y/t/d/gd/theta)
   type        = String identifying the type of DAT file
                 (wave/sediment/graindist/bedlayers/point/drifter/2d)
   
   Example:
   dims = xb_dat_dims(filename)
   
.. seealso:: xb_dat_read_, xb_dat_type_, xb_read_dat_


.. _xb_dat_read:

xb_dat_read
~~~~~~~~~~~

.. code:: text

   XB_DAT_READ  Bytewise reading of XBeach DAT files using strides
   
   Reading of XBeach DAT files. Two read methods are available: minimal
   reads and minimal memory. The former minimizes the number of fread
   calls, while the latter minimizes the amount of data read into memory.
   In case the number of reads is for both methods equal, the memory
   method is used. This method is also used if the average number of reads
   per item is less than with the read method. The method used can also be
   forced. The requested data can be determined using start and end
   indices for each dimension and strides. This approach is similar to the
   netCDF toolbox. The dimensions of the DAT file provided are in general
   in the order x,y,t. The dimension order of the output is t,y,x to match
   the netCDF conventions. The start and end indices and strides should be
   provided in t,y,x order. The result is a matrix containing the
   requested data.
   
   Preferences:
   dat_method  = Force read method (read/memory)
   
               Preferences overwrite default options (not explicitly
               defined options) and can be set and retrieved using the
               xb_setpref and xb_getpref functions.
   
   Syntax:
   dat = xb_dat_read(fname, dims, varargin)
   
   Input:
   fname       = Filename of DAT file
   dims        = Array with lengths of all dimensions in DAT file
   varargin    = start:    Start positions for reading in each dimension,
                           first item is zero
                 length:   Number of data items to be read in each
                           dimension, negative is unlimited
                 stride:   Stride to be used in each dimension
                 index:    Cell array with indices to read in each
                           dimension (overwrites start/length/stride)
                 threshold:Fraction of items to read in order to switch to
                           read method
                 maxreads: Maximum reads in memory method
                 force:    Force read method (read/memory)
   
   Output:
   dat         = Matrix with dimensions defined in dims containing
                 requested data from DAT file
   
   Example:
   dat = xb_dat_read(fname, [100 3 20]);
   dat = xb_dat_read(fname, [100 3 20], 'start', 10, 'length', 90, 'stride', 2);
   dat = xb_dat_read(fname, [100 3 20], 'start', [10 1 1], 'length', [20 -1 -1], 'stride', [2 2 2]);
   
.. seealso:: xb_read_dat_, xb_read_output_, xb_dat_dims_, xb_dat_type_


.. _sec-xb_lib:

xb_lib
------

.. _xb_axes:

xb_axes
~~~~~~~

.. code:: text

   XB_AXES  Returns the data axes corresponding to a certain variable
   
   Returns all data axes corresponding to a certain variable, taking into
   account the start, stride, length and index used.
   
   Syntax:
   varargout = xb_axes(xb, var, varargin)
   
   Input:
   xb        = XBeach output structure
   var       = Variable name
   varargin  = none
   
   Output:
   varargout = Data axes
   
   Example:
   [t y x] = xb_axes(xb, 'H');
   
.. seealso:: xb_show_


.. _xb_bathy2input:

xb_bathy2input
~~~~~~~~~~~~~~

.. code:: text

   XB_BATHY2INPUT  Adds bathymetry to XBeach input structure
   
   Adds bathymetry to XBeach input structure. Also supports non-erodible
   layers.
   
   Syntax:
   xb = xb_bathy2input(xb, x, y, z, ne)
   
   Input:
   x   = x-coordinates of bathymetry
   y   = y-coordinates of bathymetry
   z   = z-coordinates of bathymetry
   ne  = non-erodible layers in bathymetry
   
   Output:
   xb  = XBeach input structure array
   
   Example:
   xb = xb_bathy2input(xb, x, y, z)
   
.. seealso:: xb_input2bathy_, xb_read_bathy_, xb_read_input_


.. _xb_check_stagger:

xb_check_stagger
~~~~~~~~~~~~~~~~

.. code:: text

   XB_CHECK_STAGGER  Compare xb_stagger output to xbeach spaceparams
   
   More detailed description goes here.
   
   Syntax:
   varargout = xb_check_stagger(xb, varargin)
   
   Input:
   varargin  =
   
   Output:
   varargout =
   
   Example:
   xb = xb_read_output('xboutput.nc', 'length', 1)
   xb_check_stagger(xb)
   
   See also

.. _xb_defpref:

xb_defpref
~~~~~~~~~~

.. code:: text

   XB_DEFPREF  Sets default preferences for XBeach Toolbox
   
   Sets default preferences for XBeach Toolbox
   
   Syntax:
   xb_defpref()
   
   Input:
   none
   
   Output:
   none
   
   Example:
   xb_defpref;
   
.. seealso:: xb_setpref_, xb_getpref_


.. _xb_dims2nc:

xb_dims2nc
~~~~~~~~~~

.. code:: text

   XB_DIMS2NC  Convert DAT formatted dimensions to NC formatted dimensions
   
   Does some dimension flipping.
   
   Syntax:
   [nc_dims dat_dims idx_dims] = xb_dims2nc(dat_dims)
   
   Input:
   dat_dims  = dimensions in DAT file
   
   Output:
   nc_dims   = dimensions in NC file
   dat_dims  = dimensions in DAT file, guaranteed to have 3 or more items
   idx_dims  = index vectors to convert DAT dims to NC dims
   
   Example:
   nc_dims = xb_dims2nc(dat_dims)
   
.. seealso:: xb_read_dims_


.. _xb_getpref:

xb_getpref
~~~~~~~~~~

.. code:: text

   XB_GETPREF  Gets values for XBeach Toolbox preferences
   
   Gets values for XBeach Toolbox preferences and initlises default
   preferences if not done yet.
   
   Syntax:
   varargout = xb_getpref(varargin)
   
   Input:
   varargin  = list of preference names
   
   Output:
   varargout = list of corresponding preference values
   
   Example:
   version = xb_getpref('version');
   [user pass] = xb_getpref('ssh_user', 'ssh_pass');
   
.. seealso:: xb_defpref_, xb_setpref_


.. _xb_getprefdef:

xb_getprefdef
~~~~~~~~~~~~~

.. code:: text

   XB_GETPREFDEF  Gets values for XBeach Toolbox preferences or supplied default
   
   Gets values for XBeach Toolbox preferences and initialises default
   preferences if not done yet. Returns supplied default if no preference
   is found.
   
   Syntax:
   varargout = xb_getprefdef(varargin)
   
   Input:
   varargin  = list of pairs with preference names and defaults
   
   Output:
   varargout = list of corresponding preference values
   
   Example:
   value = xb_getprefdef(name, default)
   
.. seealso:: xb_setpref_


.. _xb_index:

xb_index
~~~~~~~~

.. code:: text

   XB_INDEX  Makes sure that start/len/stride are of equal and right len and contain no invalid values
   
   Makes sure that start/len/stride are of equal and right len and
   contain no negative or otherwise invalid values.
   
   Syntax:
   [start len stride] = xb_index(dims, start, len, stride)
   
   Input:
   dims      = Array with dimension sizes of original data (result of
               size())
   start     = Starting indices per dimension
   len       = Length per dimension
   stride    = Strides per dimension
   
   Output:
   start     = Starting indices per dimension
   len    = len per dimension
   stride    = Strides per dimension
   
   Example:
   [start len stride] = xb_index(dims, start, len, stride)
   
.. seealso:: xb_dat_read_, xb_read_netcdf_


.. _xb_input2bathy:

xb_input2bathy
~~~~~~~~~~~~~~

.. code:: text

   XB_INPUT2BATHY  Reads bathymetry from XBeach input structure
   
   Converts XBeach input structure to a bathymetry with x, y and z values.
   Also supports reading of non-erodible layers.
   
   Syntax:
   [x y z ne] = xb_input2bathy(xb)
   
   Input:
   xb  = XBeach input structure array
   
   Output:
   x   = x-coordinates of bathymetry
   y   = y-coordinates of bathymetry
   z   = z-coordinates of bathymetry
   ne  = non-erodible layers in bathymetry
   
   Example:
   [x y z] = xb_input2bathy(xb)
   
.. seealso:: xb_bathy2input_, xb_read_bathy_, xb_read_input_


.. _xb_reference:

xb_reference
~~~~~~~~~~~~

.. code:: text

   XB_REFERENCE  Creates a WIKI page with a params.txt reference
   
   Creates a WIKI page with a params.txt reference
   
   Syntax:
   xb_reference(filename, varargin)
   
   Input:
   filename  = Filename of generated file
   varargin  = type:       Type of file (wiki)
   
   Output:
   none
   
   Example:
   xb_reference('reference.txt')
   
.. seealso:: xb_get_params_


.. _xb_scale:

xb_scale
~~~~~~~~

.. code:: text

   XB_SCALE  Scales XBeach model input according to Vellings (1986)
   
   Scales XBeach model input according to Vellings (1986). All scaling
   dependent parameters should be present in the model input structure.
   
   Syntax:
   xb = xb_scale(xb, varargin)
   
   Input:
   xb          = XBeach input structure
   varargin    = depthscale:   depthscale nd
                 contraction:  horizontal contraction S
                 zmin:         minimal z-value
   
   Output:
   xb          = Scaled XBeach input structure
   
   Example:
   xb = xb_scale(xb, 'depthscale', 40, 'contraction', 1.68)
   xb = xb_scale(xb, 'depthscale', 40, 'contraction', 1.68, 'zmin', 0)
   
.. seealso:: xb_generate_model_


.. _xb_setpref:

xb_setpref
~~~~~~~~~~

.. code:: text

   XB_SETPREF  Sets customised preferences for XBeach Toolbox
   
   Sets customised preferences for XBeach Toolbox and initialises default
   preferences if not done yet.
   
   Syntax:
   xb_setpref(varargin)
   
   Input:
   varargin  = name/value pairs of preferences
   
   Output:
   none
   
   Example:
   xb_setpref('interactive', false);
   xb_setpref('interactive', false, 'ssh_user', ' ... ', 'ssh_pass', ' ... ');
   
.. seealso:: xb_defpref_, xb_getpref_


.. _xb_stagger:

xb_stagger
~~~~~~~~~~

.. code:: text

   XB_STAGGER  Computes staggered grids and gridcell properties from XBeach grid
   
   Computes staggered grids for u-, v- and c-points and gridcell
   properties like boundary lengths, cell areas and orientations from
   XBeach grid. Works similar to the gridprops function from XBeach
   itself.
   
   The last character of each variable name indicates the location within
   a grid cell for which the value is computed (z, u, v or c, see
   illustration). The variable name further consists of x or y,
   indicating a location in world coordinates, or of ds or dn, indicating
   distances or surfaces in grid coordinates. The dsdn* variables are grid
   cell surfaces. The alfa* variables are orientations of the specified
   points.
   
        coast
   
      |   |   |
   ---+---c-u-+---  ^     ^
      |   v z |     | ds  |
   ---+---+---+---  v     | s
      |   |   |           |
   ---+---+---+---
      |   |   |
   
          <--->
            dn
   
      <-------
          n
   
         sea
   
   Syntax:
   g = xb_stagger(x, y)
   
   Input:
   x       = x-coordinates of z-points
   y       = y-coordinates of z-points
   
   Output:
   g       = structure with grid information
   
   Example:
   g = xb_stagger(x,y);
   
.. seealso:: xb_generate_bathy_


.. _xb_swap:

xb_swap
~~~~~~~

.. code:: text

   XB_SWAP  Swap dimensions of XBeach output matrices from new to old convention and back
   
   Swaps dimensions of matrices in XBeach output structure from the order
   t,y,x to x,y,t and back. Current order is determined based on dimension
   information included in the XBeach structure (DIMS). If the current
   dimension order cannot be determined, the given order is used. Usage of
   the given order can also be forced. Returns the modified XBeach
   structure.
   
   Syntax:
   xb = xb_swap(xb, varargin)
   
   Input:
   xb          = XBeach output structure
   varargin    = order:        Current dimension order (tyx/xyt)
                 force:        Boolean to determine whether given order
                               should be used in all cases or not
   
   Output:
   xb          = Modified XBeach output structure
   
   Example:
   xb = xb_swap(xb)
   xb = xb_swap(xb, 'order', 'xyt', 'force', true)
   
.. seealso:: xb_read_output_


.. _xb_transects2grid:

xb_transects2grid
~~~~~~~~~~~~~~~~~

.. code:: text

   UNTITLED  Interpolate Jarkus transects on grid to efficiently setup XBeach model.
   
   Jarkus resolution is much higher in cross shore direction than in longshore
   direction. This script spans a grid with y-direction parallel to the
   mean coastline and the x-axis shore normal to it. The resulting grid with 
   course lonsghore resolution and fine crosshore resolution is used in griddata 
   to efficiently interpolate transects measurments on a grid. 
   
   Syntax:
   [xJ,yJ,zJ] = xb_transects2grid(transects)
   
   Input:
   varargin  = transects structure 
   
   Output:
   varargout = gridded Jarkus data
   
   Example:
   transects = jarkus_transects('id', [6002000:6002900],'output',{'id','time','x','y','cross_shore','altitude','angle'}); %
   transects2 = jarkus_interpolatenans(transects);
   transects3 = jarkus_merge(transects2,'dim','time');
   [xJ,yJ,zJ] = xb_transects2grid(transects3)
   
   See also 

.. _xb_upgrade_1d:

xb_upgrade_1d
~~~~~~~~~~~~~

.. code:: text

   XB_UPGRADE_1D  Converts an old 1D model to a superfast 1D model
   
   Converts an old 1D model with 2 lateral grids to a superfast 1D model
   withou lateral grids. An XBeach input structure or path to an XBeach
   model can be supplied. In the latter case, the model is overwritten
   with the new setup.
   
   Syntax:
   xb = xb_upgrade_1d(xb, varargin)
   
   Input:
   xb        = XBeach input structure or path to XBeach model
   varargin  = none
   
   Output:
   xb        = modified XBeach input structure
   
   Example:
   xb = xb_upgrade_1d(xb)
   xb = xb_upgrade_1d(path)
   
.. seealso:: xb_write_input_


.. _xb_verbose:

xb_verbose
~~~~~~~~~~

.. code:: text

   XB_VERBOSE  Writes verbose messages fro XBeach toolbox
   
   Writes verbose messages fro XBeach toolbox. Excepts any fprintf like
   Syntax:
   xb_verbose(varargin)
   
   Input:
   
   Input:
   varargin  = fprintf like input
   
   Output:
   none
   
   Example:
   xb_verbose('File not found [%s]', filename)
   
   See also fprintf

.. _sec-xb_modelsetup:

xb_modelsetup
-------------

.. _xb_generate_model:

xb_generate_model
~~~~~~~~~~~~~~~~~

.. code:: text

   XB_GENERATE_MODEL  Generates a XBeach structure with a full model setup
   
   Generates a XBeach structure with a full model setup. By default this
   is a minimal setup with default bathymetry, boundary conditions and
   settings. The defaults can be overwritten by supplying cell arrays with
   settings for either the bathymetry, waves, tide or model settings. The
   result is a XBeach structure, which can be written to disk easily.
   
   Syntax:
   varargout = xb_generate_model(varargin)
   
   Input:
   varargin  = bathy:      cell array of name/value pairs of bathymetry
                           settings supplied to xb_generate_grid
               waves:      cell array of name/value pairs of waves
                           settings supplied to xb_generate_waves
               tide:       cell array of name/value pairs of tide
                           settings supplied to xb_generate_tide
               wavegrid:   cell array of name/value pairs of tide
                           settings supplied to xb_generate_wavedirgrid
               settings:   cell array of name/value pairs of model
                           settings supplied to xb_generate_settings
               write:      boolean that indicates whether model setup
                           whould be written to disk (default: false)
               path:       destination directory of model setup, if
                           written to disk
               createwavegrid: bool used to determine whether this
                           function calls the xb_generate_wavegrid
                           function. In case of long crested waves one can
                           think of turning off automatic generation of
                           the wave grid.
   
   Output:
   xb        = XBeach structure array
   
   Example:
   xb = xb_generate_model();
   xb = xb_generate_model('write', false);
   xb = xb_generate_model('bathy', {'x', [ ... ], 'z', [ ... ]}, 'waves', {'Hm0', 9, 'Tp', 18});
   
.. seealso:: xb_generate_settings_, xb_generate_bathy_, xb_generate_waves_, xb_generate_tide_, xb_write_input_


.. _xb_generate_settings:

xb_generate_settings
~~~~~~~~~~~~~~~~~~~~

.. code:: text

   XB_GENERATE_SETTINGS  Generates a XBeach structure with model settings
   
   Generates a XBeach structure with model settings. A minimal set of
   default settings is used, unless otherwise provided. Settings can be
   provided by a varargin list of name/value pairs.
   
   Syntax:
   xb = xb_generate_settings(varargin)
   
   Input:
   varargin  = Name/value pairs of model settings (e.g. 'nx',100,'ny',200)
   
   Output:
   xb        = XBeach structure array
   
   Example:
   xb = xb_generate_settings()
   xb = xb_generate_settings('nx', 100, 'ny', 200)
   
.. seealso:: xb_generate_model_


.. _sec-xb_modelsetup/xb_bc:

xb_modelsetup/xb_bc
-------------------

.. _xb_bc_extracttp:

xb_bc_extracttp
~~~~~~~~~~~~~~~

.. code:: text

   XB_BC_EXTRACTTP  Extracts wave period from XBeach input structure
   
   Extracts wave period from XBeach input structure
   
   Syntax:
   Tp = xb_bc_extracttp(xb)
   
   Input:
   xb          = XBeach input structure
   
   Output:
   Tp          = peak wave period
   
   Example:
   Tp = xb_bc_extracttp(xb)
   
.. seealso:: xb_bc_extractwl_, xb_generate_model_


.. _xb_bc_extractwl:

xb_bc_extractwl
~~~~~~~~~~~~~~~

.. code:: text

   XB_BC_EXTRACTWL  Extracts water level from XBeach input structure
   
   Extracts water level from XBeach input structure
   
   Syntax:
   wl = xb_bc_extractwl(xb)
   
   Input:
   xb          = XBeach input structure
   
   Output:
   wl          = water level
   
   Example:
   wl = xb_bc_extractwl(xb)
   
.. seealso:: xb_bc_extracttp_, xb_generate_model_


.. _xb_generate_tide:

xb_generate_tide
~~~~~~~~~~~~~~~~

.. code:: text

   XB_GENERATE_TIDE  Generates XBeach structure with tide data
   
   Generates a XBeach input structure with tide settings. A minimal set of
   default settings is used, unless otherwise provided. Settings can be
   provided by a varargin list of name/value pairs.
   
   WARNING: to define two different timeseries on the lateral boundaries
            without specifying waterlevels at the back, set the paulrevere
            option in your parameter setting to 1 and define the lateral
            boundary conditions as if front and back in this function.
   
   Syntax:
   xb = xb_generate_tide(varargin)
   
   Input:
   varargin  = time:   array of starttimes of tide period in seconds
               front:  array of waterlevels at seaward model border
               back:   array of waterlevels at landward model border
   
   Output:
   xb        = XBeach structure array
   
   Example:
   xb = xb_generate_tide()
   xb = xb_generate_tide('front', 10, 'back', 5)
   xb = xb_generate_tide('time', [0 1800 3600], 'front', [5 10 5], 'back', [5 5 5])
   xb = xb_generate_tide('time', [0 1800 3600], 'front', [5 10 5; 4 10 4]', 'back', [5 5 5; 5 5 5]')
   
.. seealso:: xb_generate_model_


.. _xb_generate_waves:

xb_generate_waves
~~~~~~~~~~~~~~~~~

.. code:: text

   XB_GENERATE_WAVES  Generates XBeach structure with waves data
   
   Generates a XBeach input structure with waves settings. A minimal set
   of default settings is used, unless otherwise provided. Settings can be
   provided by a varargin list of name/value pairs. The settings depend on
   the type of waves genarated (jonswap or vardens), which is indicated by
   the type parameter. The result is a XBeach structure, an instat number
   and, if necessary another XBeach structure containing the swtable.
   
               options for jonswap:
               Hm0:        significant wave height (default: 7.6)
               Tp:         peak wave period (default: 12)
               mainang     main wave direction (default: 270)
               gammajsp:   peak-enhancement factor (default: 3.3)
               s:          power in cosinus wave spreading (default: 20)
               fnyq:       Nyquist frequency (default: 1)
   
               options for vardens:
               freqs:      array of frequencies
               dirs:       array of directions
               vardens:    matrix of the size [length(dirs) length(freqs)]
                           containing variance densities
   
   Syntax:
   [xb instat swtable] = xb_generate_waves(varargin)
   
   Input:
   varargin  = type:       type of waves to be generated (jonswap/vardens)
               duration:   array with durations in seconds
               timestep:   array with timesteps in seconds
   
   Output:
   xb        = XBeach structure array
   
   Example:
   xb = xb_generate_waves()
   xb = xb_generate_waves('Hm0', 9, 'Tp', 18)
   xb = xb_generate_waves('Hm0', [7 9 7], 'Tp', [12 18 12], 'duration', [1800 3600 1800])
   xb = xb_generate_waves('type', 'vardens', 'freqs', [ ... ], 'dirs', [ ... ])
   
.. seealso:: xb_generate_model_


.. _xb_set_start_time:

xb_set_start_time
~~~~~~~~~~~~~~~~~

.. code:: text

   XB_SET_START_TIME  sets the required start and stop times
   
   Routine to estimate the time needed for the waves to reach te coast.
   Based on that the tstart, morstart and tstop can be adjusted. The
   required time is estimated assuming a wave celerity of sqrt(gh).
   
   Syntax:
   varargout = xb_set_start_time(varargin)
   
   Input:
   xb        = XBeach input structure
   varargin  = propertyname-propertyvaluepairs
               - waterlevel : maximum water level during the simulation
   
   Output:
   xb        = XBeach input structure
   
   Example:
   xb_set_start_time
   
   See also

.. _sec-xb_modelsetup/xb_grid:

xb_modelsetup/xb_grid
---------------------

.. _xb_generate_bathy:

xb_generate_bathy
~~~~~~~~~~~~~~~~~

.. code:: text

   XB_GENERATE_BATHY  Creates a model bathymetry
   
   Creates a model bathymetry in either one or two dimensions based on a
   given bathymetry. The result is a XBeach input structure containing
   three matrices of equal size containing a rectilinear grid in x, y and
   z coordinates.
   
   Syntax:
   xb = xb_generate_bathy(varargin)
   
   Input:
   varargin  = x:          x-coordinates of bathymetry
               y:          y-coordinates of bathymetry
               z:          z-coordinates of bathymetry
               ne:         vector or matrix of the size of z containing
                           either booleans indicating if a cell is
                           non-erodable or a numeric value indicating the
                           thickness of the erodable layer on top of a
                           non-erodable layer
               xgrid:      options for xb_grid_xgrid
               ygrid:      options for xb_grid_ygrid
               rotate:     boolean flag that determines whether the
                           coastline is located in line with y-axis
               crop:       either a boolean indicating if grid should be
                           cropped to obtain a rectangle or a [x y w h]
                           array indicating how the grid should be cropped
               finalise:   either a boolean indicating if grid should be
                           finalised using default settings or a cell
                           array indicating the finalisation actions to
                           perform
               posdwn:     boolean flag that determines whether positive
                           z-direction is down
               zdepth:     extent of model below mean sea level, which is
                           used if non-erodable layers are defined
               optimize:   boolean to enable optimization of the grid, if
                           switched on, the provided grid is interpreted as
                           bathymetry and an optimal grid is defined, if
                           switched off, the provided grid is used as is.
               world_coordinates:  
                           boolean to enable a grid defined in
                           world coordinates rather than XBeach model 
                           coordinates
               superfast:  boolean to enable superfast 1D mode
   
   Output:
   xb        = XBeach structure array
   
   Example:
   xb = xb_generate_bathy('x', x, 'y', y, 'z', z)
   
.. seealso:: xb_grid_xgrid_, xb_grid_ygrid_, xb_generate_model_


.. _xb_generate_wavedirgrid:

xb_generate_wavedirgrid
~~~~~~~~~~~~~~~~~~~~~~~

.. code:: text

   UNTITLED  Generates wave directional grid for wave action balance.
   
   More detailed description goes here.
   
   Syntax:
   xb          = xb_directional_wavegrid(xb,varargin)
   
   Input:
   varthr      = threshold (percentage of maximum variance) that is taken into account to set-up wavedir grid
   nbins       = number of directional bins
   normal      = direction normal to the shore line (optional)
   plot        = 0 = no plot, 1 = plot of directional wave grid
   
   Output:
   xb          = XBeach structure array (with adapted theta
   
   Example:
   %
   waves = xb_generate_waves
   xb = xs_join(xb, waves);
   exmaple 1
   xb = xs_empty(); xb = xs_set(xb,'alpha',0,'dir',[270],'s',[5]); xb = xb_generate_wavedirgrid(xb);
   example 2
   xb = xs_empty(); xb = xs_set(xb,'alpha',-33,'dir',[273],'s',[5]); xb = xb_generate_wavedirgrid(xb);
   example 3
   xb = xs_empty(); xb = xs_set(xb,'alpha',-33,'dir',[231 241 258 273],'s',[2 10 2 5]); xb = xb_generate_wavedirgrid(xb);
   
   See also

.. _xb_grid_add:

xb_grid_add
~~~~~~~~~~~

.. code:: text

   XB_GRID_ADD  Finalise grid and determine properties
   
   Finalizes a given grid and determines dimensions and other properties.
   The result is stored in an XBeach structure that can be used as model
   Syntax:
   xb = xb_grid_add(varargin)
   
   Input:
   
   Input:
   varargin  = x:          x-coordinates of bathymetry
               y:          y-coordinates of bathymetry
               z:          z-coordinates of bathymetry
               ne:         vector or matrix of the size of z containing
                           either booleans indicating if a cell is
                           non-erodable or a numeric value indicating the
                           thickness of the erodable layer on top of a
                           non-erodable layer
               posdwn:     boolean flag that determines whether positive
                           z-direction is down
               zdepth:     extent of model below mean sea level, which is
                           used if non-erodable layers are defined
               superfast:  boolean to enable superfast 1D mode
   
   Output:
   xb        = XBeach structure array
   
   Example:
   xb = xb_grid_add('x', x, 'y', y, 'z', z);
   
.. seealso:: xb_generate_bathy_, xb_grid_optimize_


.. _xb_grid_crop:

xb_grid_crop
~~~~~~~~~~~~

.. code:: text

   XB_GRID_CROP  Automatically crops grid minimzing the number of NaN's and specifies extent of cropped area
   
   Returns the extent of a cropped grid within a supplied grid. The
   cropped area can be supplied by a vector indicating the origin, width
   and height of the area. If no area is supplied, the largest area with a
   minimum of NaN's (approximately) is used.
   
   TODO: optimize the auto-crop algorithm
   
   Syntax:
   [xmin xmax ymin ymax] = xb_grid_crop(x, y, z, varargin)
   
   Input:
   x           = x-coordinates of grid to be cropped
   y           = y-coordinates of grid to be cropped
   z           = elevations of grid to be cropped
   varargin    = crop:     vector like [x y w h] containing the origin,
                           width and height of the cropped area
   
   Output:
   xmin        = minimum x-coordinate of cropped area
   xmax        = maximum x-coordinate of cropped area
   ymin        = minimum y-coordinate of cropped area
   xmax        = maximum y-coordinate of cropped area
   
   Example:
   [xmin xmax ymin ymax] = xb_grid_crop(x, y, z)
   [xmin xmax ymin ymax] = xb_grid_crop(x, y, z, 'crop', [x0 y0 w h])
   
.. seealso:: xb_generate_grid_, xb_grid_extent_


.. _xb_grid_delft3d:

xb_grid_delft3d
~~~~~~~~~~~~~~~

.. code:: text

   XB_GRID_DELFT3D  Convert XBeach grid to Delft3D and back
   
   Accepts a path to an XBeach model or an XBeach input structure. Either
   way it returns an XBeach structure with the grid definition swapped
   from XBeach format to Delft3D format or vice versa. In case a path is
   given, the written model is updated as well.
   
   Syntax:
   xb = xb_grid_delft3d(varargin)
   
   Input:
   varargin  = Either an XBeach input structure or path to XBeach model
   
   Output:
   xb        = Modified XBeach input structure
   
   Example:
   xb_grid_delft3d('path_to_model/')
   xb_grid_delft3d('path_to_model/')
   xb = xb_grid_delft3d('path_to_model/')
   xb = xb_grid_delft3d(xb)
   
.. seealso:: xb_generate_grid_


.. _xb_grid_extent:

xb_grid_extent
~~~~~~~~~~~~~~

.. code:: text

   XB_GRID_EXTENT  Determines the extent and minimum cellsize of a specified grid
   
   Determines the minimum and maximum values of the x and y coordinates of
   a specified grid and the minimum cell size as well.
   
   Syntax:
   [xmin xmax ymin ymax cellsize] = xb_grid_extent(x, y, varargin)
   
   Input:
   x           = x-coordinates of grid to be cropped
   y           = y-coordinates of grid to be cropped
   varargin    = none
   
   Output:
   xmin        = minimum x-coordinate of grid
   xmax        = maximum x-coordinate of grid
   ymin        = minimum y-coordinate of grid
   xmax        = maximum y-coordinate of grid
   cellsize    = minimum cellsize in grid
   
   Example:
   [xmin xmax ymin ymax cellsize] = xb_grid_extent(x, y)
   
.. seealso:: xb_grid_resolution_


.. _xb_grid_finalise:

xb_grid_finalise
~~~~~~~~~~~~~~~~

.. code:: text

   XB_GRID_FINALISE  Performs several finalisation actions on an XBeach grid
   
   Performs several finalisation actions on an XBeach grid, like extending
   and flattening boundaries to prevent numerical instabilities in the
   calculation.
   
                           currently available actions:
                               lateral_extend:     copy lateral boundaries
                               lateral_seawalls:   close dry lateral
                                                   boundaries with
                                                   sandwalls
                               seaward_flatten:    flatten offshore
                                                   boundary
                               landward_extend:    extend lanward border
                                                   with specified
                                                   elevation 
                               seaward_extend:     extend seaward border
                                                   to a certain depth
   
                 cells:    number of cells to use in each action
   
   Preferences:
   grid_finalise   = Cell array with finalisation options (see options)
   
               Preferences overwrite default options (not explicitly
               defined options) and can be set and retrieved using the
               xb_setpref and xb_getpref functions.
   
   Syntax:
   [x y z] = xb_grid_finalise(x, y, z, varargin)
   
   Input:
   x           = x-coordinates of grid to be finalised
   y           = y-coordinates of grid to be finalised
   z           = elevations of grid to be finalised
   varargin    = actions:  cell array containing strings indicating the
                           order and actions to be performed
   
   Output:
   x           = x-coordinates of finalised grid
   y           = y-coordinates of finalised grid
   z           = elevations of finalised grid
   
   Example:
   [x y z] = xb_grid_finalise(x, y, z)
   [x y z] = xb_grid_finalise(x, y, z, 'actions', {'landward_polder' 'lateral_sandwalls' 'lateral_extend' 'seaward_flatten'})
   
.. seealso:: xb_generate_grid_


.. _xb_grid_interpolate:

xb_grid_interpolate
~~~~~~~~~~~~~~~~~~~

.. code:: text

   XB_GRID_INTERPOLATE  Interpolates a 2D grid on another
   
   This function is equal to the INTERP2 function in case of an not
   rotated orthogonal input grid. However, if the input grid is rotated,
   the INTERP2 function fails. This function rotates both the input as the
   Syntax:
   zi = xb_grid_interpolate(x, y, z, zi, yi, varargin)
   
   Input:
   x           = x-coordinates of input grid
   y           = y-coordinates of input grid
   z           = elevations of input grid
   xi          = x-coordinates of output grid
   yi          = y-coordinates of output grid
   varargin    = precision:    Rotation precision in case of rotated grids
   
   Output:
   aligned with the coordinate axes. Subsequently, the INTERP2 function is
   used for interpolation.
   
   Output:
   zi          = elevations of output grid
   
   Example:
   zi = xb_grid_interpolate(x, y, z, xi, yi)
   
.. seealso:: xb_grid_merge_


.. _xb_grid_merge:

xb_grid_merge
~~~~~~~~~~~~~

.. code:: text

   XB_GRID_MERGE  Merges two or more 2D grids together
   
   Merges two or more 2D grids together by defining an output rectangular,
   orthogonal and equidistant output grid based on the smallest grid size
   in the input grids. The input grids are then interpolated on the output
   grid. The first grid will end up below the others, the last grid will
   end up on top of the others.
   
   Syntax:
   [x y z] = xb_grid_merge(varargin)
   
   Input:
   varargin  = x:          cell array with x-coordinate vectors or
                           matrices of input grids
               y:          cell array with y-coordinate vectors or
                           matrices of input grids
               z:          cell array with elevation matrices of input
                           grids
   
   Output:
   x           = x-coordinates of merged grid
   y           = y-coordinates of merged grid
   z           = elevations in merged grid
   
   Example:
   [x y z] = xb_grid_merge('x',{x1 x2 x3},'y',{y1 y2 y3},'z',{z1 z2 z3})
   
.. seealso:: xb_grid_rotate_


.. _xb_grid_optimize:

xb_grid_optimize
~~~~~~~~~~~~~~~~

.. code:: text

   XB_GRID_OPTIMIZE  Creates a model grid based on a given bathymetry
   
   Creates a model grid in either one or two dimensions based on a given
   bathymetry. The result is three matrices of equal size containing a
   rectilinear grid in x, y and z coordinates.
   
   Syntax:
   xb = xb_grid_optimize(varargin)
   
   Input:
   varargin  = x:          x-coordinates of bathymetry
               y:          y-coordinates of bathymetry
               z:          z-coordinates of bathymetry
               ne:         vector or matrix of the size of z containing
                           either booleans indicating if a cell is
                           non-erodable or a numeric value indicating the
                           thickness of the erodable layer on top of a
                           non-erodable layer
               xgrid:      options for xb_grid_xgrid
               ygrid:      options for xb_grid_ygrid
               rotate:     boolean flag that determines whether the
                           coastline is located in line with y-axis
               crop:       either a boolean indicating if grid should be
                           cropped to obtain a rectangle or a [x y w h]
                           array indicating how the grid should be cropped
               finalise:   either a boolean indicating if grid should be
                           finalised using default settings or a cell
                           array indicating the finalisation actions to
                           perform
               posdwn:     boolean flag that determines whether positive
                           z-direction is down
               world_coordinates:  
                           boolean to enable a grid defined in
                           world coordinates rather than XBeach model 
                           coordinates
               zdepth:     extent of model below mean sea level, which is
                           used if non-erodable layers are defined
   
   Output:
   xb        = XBeach structure array
   
   Example:
   xb = xb_grid_optimize('x', x, 'y', y, 'z', z)
   
.. seealso:: xb_grid_xgrid_, xb_grid_ygrid_, xb_generate_bathy_, xb_generate_model_


.. _xb_grid_orientation:

xb_grid_orientation
~~~~~~~~~~~~~~~~~~~

.. code:: text

   XB_GRID_ORIENTATION  Determines the orientation of a 2D bathymetric grid
   
   Determines the dimension and direction which runs from sea landward in
   a 2D bathymetric grid. Based on the maximum mean minimum slope, the
   dimension which runs cross-shore is chosen. Subsequently, the direction
   is determined. If the grid runs from sea landward in negative
   direction, the direction is negative, otherwise positive.
   
   Syntax:
   [dim dir] = xb_grid_orientation(x, y, z, varargin)
   
   Input:
   x           = x-coordinates of bathymetric grid
   y           = y-coordinates of bathymetric grid
   z           = elevations in bathymetric grid
   varargin    = none
   
   Output:
   dim         = corss-shore dimension in bathymetric grid (1/2)
   dir         = direction in bathymetric grid that runs from sea landward
                 (1/-1)
   
   Example:
   [dim dir] = xb_grid_orientation(x, y, z)
   
.. seealso:: xb_generate_grid_


.. _xb_grid_resolution:

xb_grid_resolution
~~~~~~~~~~~~~~~~~~

.. code:: text

   XB_GRID_RESOLUTION  Determines the maximum cellsize of a regular grid of a certain size fitting the extent of the specified grid
   
   Determines the extent of the specified grid and the maximum cellsize to
   generate a regular grid spanning this extent without exceeding a
   certain size (bytes).
   
   Syntax:
   [cellsize xmin xmax ymin ymax] = xb_grid_resolution(x, y, varargin)
   
   Input:
   x           = x-coordinates of grid to be cropped
   y           = y-coordinates of grid to be cropped
   varargin    = maxsize:  maximum grid size in bytes (default 10MB), in
                           case the value is 'max' the maximum size in the
                           currently available memory space is used.
   
   Output:
   cellsize    = maximum cellsize in regular grid
   xmin        = minimum x-coordinate of regular grid
   xmax        = maximum x-coordinate of regular grid
   ymin        = minimum y-coordinate of regular grid
   xmax        = maximum y-coordinate of regular grid
   
   Example:
   [cellsize xmin xmax ymin ymax] = xb_grid_resolution(x, y)
   [cellsize xmin xmax ymin ymax] = xb_grid_resolution(x, y, 'maxsize', 1024^3)
   
.. seealso:: xb_grid_extent_


.. _xb_grid_rotate:

xb_grid_rotate
~~~~~~~~~~~~~~

.. code:: text

   XB_GRID_ROTATE  Rotates a grid around an origin
   
   Rotates a grid around an origin. The origin can be specified.
   
   Syntax:
   [xr yr] = xb_grid_rotate(x, y, alpha, varargin)
   
   Input:
   x           = x-coordinates of bathymetric grid
   y           = y-coordinates of bathymetric grid
   alpha       = rotation angle
   varargin    = origin:   origin for rotation
                 units:    input units (degrees/radians)
   
   Output:
   xr          = x-coordinates of rotated grid
   yr          = y-coordinates of rotated grid
   
   Example:
   [xr yr] = xb_grid_rotate(x, y, alpha)
   
.. seealso:: xb_grid_rotation_


.. _xb_grid_rotation:

xb_grid_rotation
~~~~~~~~~~~~~~~~

.. code:: text

   XB_GRID_ROTATION  Determines rotation of a 2D grid based on the coastline
   
   Determines the location of a 2D grid based on the coastline by
   detecting the coastline and determining the angle of the coastline.
   
   Syntax:
   [alpha a b] = xb_grid_rotation(x, y, z, varargin)
   
   Input:
   x           = x-coordinates of bathymetric grid
   y           = y-coordinates of bathymetric grid
   z           = elevations in bathymetric grid
   varargin    = units:    output units (degrees/radians)
   
   Output:
   alpha       = rotation of grid
   a           = linear regression parameter of coastline (y=a+b*x)
   b           = linear regression parameter of coastline (y=a+b*x)
   
   Example:
   alpha = xb_grid_rotation(x, y, z)
   
.. seealso:: xb_grid_rotate_


.. _xb_grid_trim:

xb_grid_trim
~~~~~~~~~~~~

.. code:: text

   XB_GRID_TRIM  Removes all empty rows and columns from a 2D grid
   
   Removes all rows and columns in a 2D grid containing NaN's only.
   
   Syntax:
   [x y z] = xb_grid_trim(x, y, z, varargin)
   
   Input:
   x           = x-coordinates of grid to be trimmed
   y           = y-coordinates of grid to be trimmed
   z           = elevations of grid to be trimmed
   varargin    = none
   
   Output:
   x           = x-coordinates of trimmed grid
   y           = y-coordinates of trimmed grid
   z           = elevations of trimmed grid
   
   Example:
   [x y z] = xb_grid_trim(x, y, z)
   
.. seealso:: xb_grid_extent_


.. _xb_grid_world2xb:

xb_grid_world2xb
~~~~~~~~~~~~~~~~

.. code:: text

   XB_GRID_WORLD2XB  Rotates a grid in world coordinates to XBeach coordinates
   
   Rotates a grid in world coordinates to XBeach coordinates
   
   Syntax:
   [x y] = xb_grid_world2xb(x, y, xori, yori, alpha)
   
   Input:
   x           = x-coordinates
   y           = y-coordinates
   xori        = x-origin
   yori        = y-origin
   alpha       = grid rotation
   
   Output:
   x           = x-coordinates
   y           = y-coordinates
   
   Example:
   [x y] = xb_grid_world2xb(x, y, xori, yori, alpha)
   
.. seealso:: xb_grid_xb2world_


.. _xb_grid_xb2world:

xb_grid_xb2world
~~~~~~~~~~~~~~~~

.. code:: text

   XB_GRID_XB2WORLD  Rotates a grid in XBeach coordinates to world coordinates
   
   Rotates a grid in XBeach coordinates to world coordinates
   
   Syntax:
   [x y] = xb_grid_xb2world(x, y, xori, yori, alpha)
   
   Input:
   x           = x-coordinates
   y           = y-coordinates
   xori        = x-origin
   yori        = y-origin
   alpha       = grid rotation
   
   Output:
   x           = x-coordinates
   y           = y-coordinates
   
   Example:
   [x y] = xb_grid_xb2world(x, y, xori, yori, alpha)
   
.. seealso:: xb_grid_world2xb_


.. _xb_grid_xgrid:

xb_grid_xgrid
~~~~~~~~~~~~~

.. code:: text

   XB_GRID_XGRID  Creates a model grid in x-direction based on bathymetry
   
   Function to interpolate (no extrapolation) profile measurements to
   cross shore consant or varying grid for an XBeach profile model. Cross
   shore grid size is limited by user-defined minimum grid size in shallow
   water and land, long wave resolution on offshore boundary, depth to
   grid size ratio and grid size smoothness constraints. The function uses
   the Courant condition to find the optimal grid size given these
   constraints.
   
   Optional input in keyword,value pairs
     - xgrid    :: [m] vector which defines the variable xgrid
     - Tm       :: [s] incident short wave period (used for maximum grid size at offshore boundary) 
                       if you impose time series of wave conditions use the min(Tm) as input (default = 5)
     - dxmin    :: [m] minimum required cross shore grid size (usually over land) (default = 1)
     - dxmax    :: [m] user-specified maximum grid size, when usual wave
                       period / CFL condition does not suffice (default Inf)
     - vardx    :: [-] 0 = constant dx, 1 = varying dx (default = 1)
     - g        :: [ms^-2] gravity constant (default = 9.81)
     - CFL      :: [-] Courant number in grid generator (default = 0.9)
     - dtref    :: [-] Ref value for dt in computing dx from CFL (default = 4)
     - maxfac   :: [-] Maximum allowed grid size ratio between adjacent cells (default = 1.15)
     - dy, 5    :: [m] dy (default = 5)
     - wl,0     :: [m] water level elevation relative to bathymetry used to estimate water depth (default = 0)
     - depthfac :: [-] Maximum gridsize to water depth ratio (default = 2)
   
   Syntax:
   [xgr zgr] = xb_grid_xgrid(xin, zin, varargin)
   
   Input:
   xin   = vector with cross-shore coordinates; increasing from zero
   towards shore
   zin   = vector with bed levels; positive up
   
   Output:
   xgr   = x-grid coordinates
   zgr   = bed elevations
   
   Example:
   [xgr zgr] = xb_grid_xgrid([0:1:200], 0.1*[0:1:200]-15);
   
.. seealso:: xb_generate_grid_, xb_grid_ygrid_


.. _xb_grid_ygrid:

xb_grid_ygrid
~~~~~~~~~~~~~

.. code:: text

   XB_GRID_YGRID  Creates a model grid in y-direction based on minimum and maximum cell size and area of interest
   
   Generates a model grid in y-direction using two grid cellsizes. The
   minimum grid cellsize is used for the area of interest. The maximum is
   used near the lateral borders. A gradual transition between the grid
   cellsizes over a specified distance is automatically generated. The
   area of interest can be defined in several manners. By default, this is
   a distance of 100m in the center of the model.
   
   Syntax:
   ygr = xb_grid_ygrid(yin, varargin)
   
   Input:
   yin       = range of y-coordinates to be included in the grid
   varargin  = dymin:                  minimum grid cellsize
               dymax:                  maximum grid cellsize
               area_type:              type of definition of the area of
                                       interest (center/range)
               area_size:              size of the area of interest
                                       (length or fraction in case of
                                       area_type center, from/to range in
                                       case of area_type range)
               transition_distance:    distance over which the grid
                                       cellsize is gradually changed from
                                       mimumum to maximum, a negative
                                       value means the distance may be
                                       adapted to limit the error made in
                                       the fit
   
   Output:
   ygr       = generated grid in y_direction
   
   Example:
   ygr = xb_grid_ygrid(yin)
   
.. seealso:: xb_generate_grid_, xb_grid_xgrid_


.. _sec-xb_modelsetup/xb_grid/xb_delft3d:

xb_modelsetup/xb_grid/xb_delft3d
--------------------------------

.. _xb_delft3d_addpath:

xb_delft3d_addpath
~~~~~~~~~~~~~~~~~~

.. code:: text

   XB_DELFT3D_ADDPATH  One line description goes here.
   
   More detailed description goes here.
   
   Syntax:
   varargout = xb_delft3d_addpath(varargin)
   
   Input:
   varargin  =
   
   Output:
   varargout =
   
   Example:
   xb_delft3d_addpath
   
   See also

.. _xb_delft3d_wldep2xb:

xb_delft3d_wldep2xb
~~~~~~~~~~~~~~~~~~~

.. code:: text

   XB_DELFT3D_WLDEP2XB  One line description goes here.
   
   More detailed description goes here.
   
   Syntax:
   varargout = xb_delft3d_wldep2xb(varargin)
   
   Input:
   varargin  =
   
   Output:
   varargout =
   
   Example:
   xb_delft3d_wldep2xb
   
   See also

.. _xb_delft3d_wlgrid2xb:

xb_delft3d_wlgrid2xb
~~~~~~~~~~~~~~~~~~~~

.. code:: text

   XB_DELFT3D_WLGRID2XB  One line description goes here.
   
   More detailed description goes here.
   
   Syntax:
   varargout = xb_delft3d_wlgrid2xb(varargin)
   
   Input:
   varargin  =
   
   Output:
   varargout =
   
   Example:
   xb_delft3d_wlgrid2xb
   
   See also

.. _xb_delft3d_xb2wldep:

xb_delft3d_xb2wldep
~~~~~~~~~~~~~~~~~~~

.. code:: text

   XB_DELFT3D_XB2WLDEP  One line description goes here.
   
   More detailed description goes here.
   
   Syntax:
   varargout = xb_delft3d_xb2wldep(varargin)
   
   Input:
   varargin  =
   
   Output:
   varargout =
   
   Example:
   xb_delft3d_xb2wldep
   
   See also

.. _xb_delft3d_xb2wlgrid:

xb_delft3d_xb2wlgrid
~~~~~~~~~~~~~~~~~~~~

.. code:: text

   XB_DELFT3D_XB2WLGRID  One line description goes here.
   
   More detailed description goes here.
   
   Syntax:
   varargout = xb_delft3d_xb2wlgrid(varargin)
   
   Input:
   varargin  =
   
   Output:
   varargout =
   
   Example:
   xb_delft3d_xb2wlgrid
   
   See also

.. _sec-xb_nesting:

xb_nesting
----------

.. _xb_nest_delft3d:

xb_nest_delft3d
~~~~~~~~~~~~~~~

.. code:: text

   %XB_NEST_DELFT3D  Obtain Delft3D boundary conditions and bathymetry from another model
   
   Obtain Delft3D boundary conditions and/or bathymetry from a WW3,
   Delft3D. This information can be used for nesting a Delft3D model into
   one of these two models.
   
   Syntax:
   info = xb_nest_delft3d(varargin)
   
   Input:
   varargin  = type:   Model type to nest into (ww3/delft3d)
               flow:   Cell array with parameters for the flow nest
                       procedure
               wave:   Cell array with parameters for the wave nest
                       procedure
               bathy:  Cell array with parameters for the bathymetry nest
                       procedure
   
   Output:
   info      = Struct with filenames of boundary condition and bathymetry
               files
   
   Example:
   info = xb_nest_delft3d('delft3d', 'flow', { ... })
   
.. seealso:: xb_nest_xbeach_, xb_delft3d_flow_, xb_delft3d_wave_, xb_delft3d_bathy_


.. _xb_nest_xbeach:

xb_nest_xbeach
~~~~~~~~~~~~~~

.. code:: text

   XB_NEST_XBEACH  Obtain XBeach boundary conditions and bathymetry from another model
   
   Obtain XBeach boundary conditions and/or bathymetry from a WW3, Delft3D
   or another XBeach model. This information can be used for nesting an
   XBeach model into one of these three models.
   
   Syntax:
   info = xb_nest_xbeach(varargin)
   
   Input:
   varargin  = type:   Model type to nest into (ww3/delft3d/xbeach)
               flow:   Cell array with parameters for the flow nest
                       procedure
               wave:   Cell array with parameters for the wave nest
                       procedure
               bathy:  Cell array with parameters for the bathymetry nest
                       procedure
   
   Output:
   info      = Struct with filenames of boundary condition and bathymetry
               files
   
   Example:
   info = xb_nest_xbeach('delft3d', 'flow', { ... })
   
.. seealso:: xb_nest_delft3d_, xb_delft3d_flow_, xb_delft3d_wave_, xb_delft3d_bathy_


.. _sec-xb_nesting/xb_delft3d:

xb_nesting/xb_delft3d
---------------------

.. _xb_bct2xb:

xb_bct2xb
~~~~~~~~~

.. code:: text

   XB_BCT2XB  Converts Delft3D-FLOW BCT file to an XBeach tide file
   
   Converts Delft3D-FLOW BCT file to an XBeach tide file, possibly
   cropping the timeseries.
   
   Syntax:
   [tidefile tidelen tideloc] = xb_bct2xb(fname, varargin)
   
   Input:
   fname     = Path to BCT file
   varargin  = tidefile:       Path to output file
               tstart:         Datenum indicating simulation start time
               tlength:        Datenum indicating simulation length
   
   Output:
   tidefile  = Path to output file
   tidelen   = Length of timeseries in output file
   tideloc   = Number of locations in output file
   
   Example:
   [tidefile tidelen tideloc] = xb_bct2xb(fname)
   [tidefile tidelen tideloc] = xb_bct2xb(fname, 'tstart', datenum('2007-11-05'))
   
.. seealso:: xb_delft3d_flow_, xb_sp22xb_


.. _xb_delft3d_bathy:

xb_delft3d_bathy
~~~~~~~~~~~~~~~~

.. code:: text

   XB_DELFT3D_BATHY  Nest Delft3D or XBeach bathymetry in another Delft3D model
   
   Nests a Delft3D bathymetry.
   
   Syntax:
   varargout = xb_delft3d_bathy(depfile, varargin)
   
   Input:
   depfile   = Path to bathymetry file
   varargin  = type:           Type of output (delft3d/xbeach)
               file:           Path to output file
   
   Output:
   varargout   = Path to output file
   
   Example:
   file = xb_delft3d_bathy(depfile)
   
.. seealso:: xb_nest_xbeach_, xb_nest_delft3d_, xb_delft3d_flow_, xb_delft3d_wave_


.. _xb_delft3d_depth:

xb_delft3d_depth
~~~~~~~~~~~~~~~~

.. code:: text

   XB_DELFT3D_DEPTH  One line description goes here.
   
   More detailed description goes here.
   
   Syntax:
   varargout = xb_delft3d_depth(varargin)
   
   Input:
   varargin  =
   
   Output:
   varargout =
   
   Example:
   xb_delft3d_depth
   
   See also

.. _xb_delft3d_flow:

xb_delft3d_flow
~~~~~~~~~~~~~~~

.. code:: text

   XB_DELFT3D_FLOW  Nest Delft3D-FLOW or XBeach model in another Delft3D-FLOW
   
   Nests a Delft3D-FLOW model based on the trih-*, BND and NST files.
   
   Syntax:
   varargout = xb_delft3d_flow(trihfile, bndfile, nstfile, varargin)
   
   Input:
   trihfile  = Path to the trih-* file to use
   bndfile   = Path to the BND file to use
   nstfile   = Path to the NST file to use
   varargin  = type:       Type of output (delft3d/xbeach)
               file:       Path to output file
               zcorr:      Water level correction
               tstart:     Datenum indicating simulation start time
               tlength:    Datenum indicating simulation length
   
   Output:
   varargout = Path to output file and dimensions of file, if applicable
   
   Example:
   file = xb_delft3d_flow(trihfile, bndfile, nstfile)
   [file tidelen tideloc] = xb_delft3d_flow('trih-csm.dat', 'kuststrook.bnd', 'kuststrook.nst', 'type', 'xbeach')
   
.. seealso:: xb_nest_xbeach_, xb_nest_delft3d_, xb_delft3d_wave_, xb_delft3d_bathy_, xb_bct2xb_


.. _xb_delft3d_grid:

xb_delft3d_grid
~~~~~~~~~~~~~~~

.. code:: text

   XB_DELFT3D_GRID  One line description goes here.
   
   More detailed description goes here.
   
   Syntax:
   varargout = xb_delft3d_grid(varargin)
   
   Input:
   varargin  =
   
   Output:
   varargout =
   
   Example:
   xb_delft3d_grid
   
   See also

.. _xb_delft3d_wave:

xb_delft3d_wave
~~~~~~~~~~~~~~~

.. code:: text

   XB_DELFT3D_WAVE  Nest Delft3D-WAVE or XBeach model in another Delft3D-WAVE model
   
   Nests a Delft3D-WAVE model by simply joining the output SP2 files into
   a single input SP2 file.
   
   Syntax:
   varargout = xb_delft3d_wave(sp2files, varargin)
   
   Input:
   sp2files  = Path to SP2 files to be joined
   varargin  = type:           Type of output (delft3d/xbeach)
               file:           Path to output file
               coordsys_in     Coordinate system of input
               coordtype_in    Coordinate system type of input
               coordsys_out    Coordinate system of output
               coordtype_out   Coordinate system type of output
               tstart:         Datenum indicating simulation start time
               tlength:        Datenum indicating simulation length
   
   Output:
   varargout   = Path to output file
   
   Example:
   file = xb_delft3d_wave(sp2files)
   file = xb_delft3d_wave('pet.loct*.sp2')
   
.. seealso:: xb_nest_xbeach_, xb_nest_delft3d_, xb_delft3d_flow_, xb_delft3d_bathy_, xb_sp22xb_


.. _xb_sp22xb:

xb_sp22xb
~~~~~~~~~

.. code:: text

   XB_SP22XB  Converts a set of Delft3D-WAVE SP2 files into XBeach wave boundary conditions
   
   Converts a set of Delft3D-WAVE SP2 files into XBeach wave boundary
   conditions by cropping the timeseries, writing the corresponding SP2
   files and a filelist file.
   
   Syntax:
   wavefile = xb_sp22xb(fnames, varargin)
   
   Input:
   fnames    = Path to Delft3D-WAVE SP2 files
   varargin  = wavefile:       Path to output file
               tstart:         Datenum indicating simulation start time
               tlength:        Datenum indicating simulation length
   
   Output:
   wavefile  = Path to output file
   
   Example:
   wavefile = xb_sp22xb('*.sp2')
   
.. seealso:: xb_delft3d_wave_, xb_bct2xb_


.. _sec-xb_nesting/xb_swan:

xb_nesting/xb_swan
------------------

.. _xb_swan_coords:

xb_swan_coords
~~~~~~~~~~~~~~

.. code:: text

   XB_SWAN_COORDS  Convert coordinates in SWAN struct
   
   Convert coordinates in SWAN struct
   
   Syntax:
   sp2 = xb_swan_coords(sp2, cs_in, ct_in, cs_out, ct_out, varargin)
   
   Input:
   sp2       = SWAN struct to be converted
   cs_in     = Coordinate system of input
   ct_in     = Coordinate system type of input
   cs_out    = Coordinate system of output
   ct_out    = Coordinate system type of output
   varargin  = none
   
   Output:
   sp2       = Converted SWAN struct
   
   Example:
   sp2 = xb_swan_coords(sp2, ... )
   
.. seealso:: xb_swan_read_


.. _xb_swan_join:

xb_swan_join
~~~~~~~~~~~~

.. code:: text

   XB_SWAN_JOIN  Join multiple files in SWAN struct to a single file
   
   Join multiple files in SWAN struct to a single file
   
   Syntax:
   sp2 = xb_swan_join(sp2, varargin)
   
   Input:
   sp2       = SWAN struct with multiple files
   varargin  = none
   
   Output:
   sp2       = SWAN struct with a single file
   
   Example:
   sp2 = xb_swan_join(sp2)
   
.. seealso:: xb_swan_split_


.. _xb_swan_read:

xb_swan_read
~~~~~~~~~~~~

.. code:: text

   XB_SWAN_READ  Read a SWAN file into a struct
   
   Read one or more SWAN files into a struct
   
   Syntax:
   sp2 = xb_swan_read(fname, varargin)
   
   Input:
   fname     = Path to SWAN file
   varargin  = none
   
   Output:
   sp2       = SWAN struct
   
   Example:
   sp2 = xb_swan_read('waves.sp2')
   sp2 = xb_swan_read('*.sp2')
   
.. seealso:: xb_swan_write_, xb_swan_struct_


.. _xb_swan_split:

xb_swan_split
~~~~~~~~~~~~~

.. code:: text

   XB_SWAN_SPLIT  Split a SWAN struct with a single file into a struct with multiple files
   
   Split a SWAN struct with a single file into a struct with multiple
   files.
   
   Syntax:
   sp2 = xb_swan_split(sp2, dim, varargin)
   
   Input:
   sp2       = SWAN struct with single file
   dim       = Dimension to be splitted (time/location/etc)
   varargin  = none
   
   Output:
   sp2       = SWAN struct with multiple files
   
   Example:
   sp2 = xb_swan_split(sp2, 'time')
   
.. seealso:: xb_swan_join_


.. _xb_swan_struct:

xb_swan_struct
~~~~~~~~~~~~~~

.. code:: text

   XB_SWAN_STRUCT  Create empty SWAN structure
   
   Create empty SWAN structure
   
   Syntax:
   sp2 = xb_swan_struct()
   
   Input:
   none
   
   Output:
   sp2       = Empty SWAN structure
   
   Example:
   sp2 = xb_swan_struct
   
.. seealso:: xb_swan_read_, xb_swan_write_


.. _xb_swan_write:

xb_swan_write
~~~~~~~~~~~~~

.. code:: text

   XB_SWAN_WRITE  Write SWAN file from struct
   
   Write one or more SWAN files from struct. In case of multiple files,
   the filename is extended with an indexing number.
   
   Syntax:
   fnames = xb_swan_write(fname, sp2, varargin)
   
   Input:
   fname     = Path to output file
   sp2       = SWAN struct to write
   varargin  = none
   
   Output:
   fnames    = Path(s) to output file(s)
   
   Example:
   fnames = xb_swan_write('waves.sp2', sp2)
   
.. seealso:: xb_swan_read_, xb_swan_struct_


.. _sec-xb_run:

xb_run
------

.. _xb_check_list:

xb_check_list
~~~~~~~~~~~~~

.. code:: text

   XB_CHECK_LIST  One line description goes here.
   
   More detailed description goes here.
   
   Syntax:
   varargout = xb_check_list(varargin)
   
   Input:
   varargin  =
   
   Output:
   varargout =
   
   Example:
   xb_check_list
   
   See also 

.. _xb_check_run:

xb_check_run
~~~~~~~~~~~~

.. code:: text

   XB_CHECK_RUN  Checks whether a run from an XBeach input structure is still running
   
   Checks whether a run started using either xb_run or xb_run_remote is
   still running. The check can be repeated using a certain time interval
   to get a warning once the run is finished.
   
   Syntax:
   runs = xb_check_run(xb, varargin)
   
   Input:
   xb        = XBeach run structure
   varargin  = repeat:     boolean flag to determine if the check is
                           repeated
               interval:   seconds between checks, if repeated
               display:    boolean to determine whether a message is
                           displayed after process has finished
               sound:      boolean to determine whether a sound is made
                           after process has finished
               callback:   callback function fired after process has
                           finished
   
   Output:
   runs      = boolean that indicates whether the job is still running
   
   Example:
   xb_check_run(xb)
   
.. seealso:: xb_run_, xb_run_remote_


.. _xb_get_bin:

xb_get_bin
~~~~~~~~~~

.. code:: text

   XB_GET_BIN  Retrieves a XBeach binary from a remote source
   
   Retrieves a XBeach binary from a remote source. By default this is the
   latest binary from the TeamCity build server. Several flavours of
   binaries exist. By default the normal win32 binary is downloaded. A
   custom host can be provided as well. Returns the location where the
   downloaded binary can be found.
   
   WARNING: SOME BINARY TYPES ARE STILL MISSING, SINCE NOT AVAILABLE IN
   TEAMCITY YET
   
   Syntax:
   fpath = xb_get_bin(varargin)
   
   Input:
   varargin  = type:       Type of binary (win32/unix/mpi/netcdf).
                           Multiple qualifiers separated by a space can be
                           used. Specifying "custom" will use the host
                           provided in the equally named varargin
                           parameter.
               host:       Host to be used in case of custom type.
   
   Output:
   fpath     = Path to downloaded executable
   
   Example:
   fpath = xb_get_bin()
   fpath = xb_get_bin('type', 'win32 mpi')
   fpath = xb_get_bin('type', 'win32 netcdf')
   fpath = xb_get_bin('type', 'win32 netcdf mpi')
   fpath = xb_get_bin('type', 'unix netcdf mpi')
   fpath = xb_get_bin('type', 'custom', 'host', ' ... ')
   
.. seealso:: xb_run_


.. _xb_resume:

xb_resume
~~~~~~~~~

.. code:: text

   XB_RESUME  Resumes XBeach model run from certain moment
   
   Combines an existing model setup with a bathymetry from the
   corresponding model output into a new model setup. Given a certain
   simulation time, the bathymetry closest to this time is obtained from
   the model output. The input bathymetry from the existing model setup is
   replaced with this evolved bathymetry. The boundary conditions and
   simulation times are adjusted to fit the bathymetry used.
   If the given simulation time is not present in the model output, the
   available simulation time closest to the given value is used.
   
   Syntax:
   xb = xb_resume(xb, varargin)
   
   Input:
   xb        = XBeach run structure or path to model setup
   varargin  = t:      simulation time at which bathymetry is read
               spinup: spinup time to use in newly generated model setup
   
   Output:
   xb        = New model setup with evolved bathymetry
   
   Example:
   xbm = xb_resume(xbr, 't', 1000)
   xbm = xb_resume(xbr, 't', 1000, 'spinup', 100)
   xbm = xb_resume('path/to/model/setup, 't', Inf)
   
.. seealso:: xb_run_, xb_run_remote_, xb_read_input_, xb_generate_model_


.. _xb_run:

xb_run
~~~~~~

.. code:: text

   XB_RUN  Runs a XBeach model locally
   
   Writes a XBeach structure to disk, retrieves a XBeach binary file and
   runs it at a certain location. Supports the use of MPI using MPICH2.
   
   TODO: MPI support
   
   Syntax:
   xb_run(xb)
   
   Input:
   xb        = XBeach input structure
   varargin  = binary:     XBeach binary to use
               nodes:      Number of nodes to use in MPI mode (1 = no mpi)
               netcdf:     Flag to use netCDF output (default: false)
               path:       Path to the XBeach model
   
   Output:
   xb        = XBeach structure array
   
   Example:
   xb_run(xb)
   xb_run(xb, 'path', 'path_to_model/')
   
.. seealso:: xb_run_remote_, xb_get_bin_


.. _xb_run_list:

xb_run_list
~~~~~~~~~~~

.. code:: text

   XB_RUN_LIST  Lists registered runs
   
   Lists runs registered by the xb_run_register function.
   
   Syntax:
   xb_run_list(varargin)
   
   Input:
   varargin  = n:      Number of log lines to show
   
   Output:
   none
   
   Example:
   xb_run_list
   xb_run_list('n', 10)
   
.. seealso:: xb_run_, xb_run_remote_, xb_run_register_


.. _xb_run_parselog:

xb_run_parselog
~~~~~~~~~~~~~~~

.. code:: text

   XB_RUN_PARSELOG  One line description goes here.
   
   More detailed description goes here.
   
   Syntax:
   varargout = xb_run_parselog(varargin)
   
   Input:
   varargin  =
   
   Output:
   varargout =
   
   Example:
   xb_run_parselog
   
   See also 

.. _xb_run_queue:

xb_run_queue
~~~~~~~~~~~~

.. code:: text

   XB_RUN_QUEUE  Queues XBeach runs for local execution
   
   Queues XBeach runs for local execution. All models are executed one at
   a time, unless set otherwise. Runs are stored in an XBeach preference
   variable and monitored using the xb_check_run function. A callback to
   this function will fire the next model.
   
   Syntax:
   xb_run_queue(varargin)
   
   Input:
   varargin  = First argument may be an XBeach input structure to be
               queued. Other arguments are expected to be from the
               following list:
               action:     Queue action to perform (add/start/next/clear)
               options:    Options to be passed to the xb_run function
               nr:         Number of parallel runs
   
   Output:
   none
   
   Example:
   xb_run_queue(xb)
   xb_run_queue(xb, 'options', {'binary', 'C:\xbeach.exe'})
   xb_run_queue
   xb_run_queue('action', 'next')
   xb_run_queue('action', 'clear')
   
.. seealso:: xb_run_, xb_check_run_


.. _xb_run_register:

xb_run_register
~~~~~~~~~~~~~~~

.. code:: text

   XB_RUN_REGISTER  Registers XBeach run in preferences database
   
   Registers an XBeach run started from xb_run or xb_run_remote in a
   preference database for future use.
   
   Syntax:
   xb_run_register(xb, varargin)
   
   Input:
   xb        = XBeach struct resulting from xb_run or xb_run_remote
   varargin  = none
   
   Output:
   none
   
   Example:
   xb_run_register(xb)
   
.. seealso:: xb_run_, xb_run_remote_, xb_run_list_


.. _xb_run_remote:

xb_run_remote
~~~~~~~~~~~~~

.. code:: text

   XB_RUN_REMOTE  Runs a XBeach model remote on the H4 cluster
   
   Writes a XBeach structure to disk, retrieves a XBeach binary file and
   runs it at a remote location accessed by SSH (by default, H4 cluster).
   Supports the use of MPI.
   
   TODO: UNIX SUPPORT
   
   Preferences:
   ssh_user        = Username for remote computer
   ssh_pass        = Password for remote computer
   path_local      = Local path to the XBeach model
   path_remote     = Path to XBeach model seen from remote computer
   
               Preferences overwrite default options (not explicitly
               defined options) and can be set and retrieved using the
               xb_setpref and xb_getpref functions.
   
   Syntax:
   xb_run_remote(xb)
   
   Input:
   xb        = XBeach structure array
   varargin  = name:       Name of the model run
               binary:     XBeach binary to use
               nodes:      Number of nodes to use in MPI mode (1 = no mpi)
               netcdf:     Flag to use netCDF output (default: false)
               ssh_host:   Host name of remote computer
               ssh_user:   Username for remote computer
               ssh_pass:   Password for remote computer
               ssh_prompt: Boolean indicating if password prompt should be
                           used
               path_local: Local path to the XBeach model
               path_remote:Path to XBeach model seen from remote computer
   
   Output:
   xb        = XBeach structure array
   
   Example:
   xb_run_remote(xb)
   xb_run_remote(xb, 'path_local', 'u:\', 'path_remote', '~/')
   
.. seealso:: xb_run_, xb_get_bin_


.. _xb_run_sh_scripts:

xb_run_sh_scripts
~~~~~~~~~~~~~~~~~

.. code:: text

   XB_RUN_SH_SCRIPTS  Run SH scripts on H4 cluster
   
   Run SH scripts on H4 cluster
   
   Preferences:
   ssh_user        = Username for remote computer
   ssh_pass        = Password for remote computer
   
               Preferences overwrite default options (not explicitly
               defined options) and can be set and retrieved using the
               xb_setpref and xb_getpref functions.
   
   Syntax:
   [job_id job_name messages] = xb_run_sh_scripts(rpath, script, varargin)
   
   Input:
   rpath     = Path where script is located seen from remote source
   script    = Name of the script to run
   varargin  = ssh_host:   Host name of remote computer
               ssh_user:   Username for remote computer
               ssh_pass:   Password for remote computer
               ssh_prompt: Boolean indicating if password prompt should be
                           used
               cd:         Boolean flag to determine if directory should
                           be changed to rpath
   
   Output:
   job_id    = Job number of process started
   job_name  = Name of process started
   messages  = Messages returned by remote source
   
   Example:
   job_id = xb_run_sh_scripts('~/', 'run.sh', 'ssh_prompt', true)
   
.. seealso:: xb_run_remote_, xb_write_sh_scripts_


.. _xb_run_unregister:

xb_run_unregister
~~~~~~~~~~~~~~~~~

.. code:: text

   XB_RUN_UNREGISTER  Unregister registered run
   
   Unregister run registered by xb_run_register
   
   Syntax:
   xb_run_unregister(id)
   
   Input:
   id        = Run index to be unregistered
   
   Output:
   none
   
   Example:
   xb_run_unregister(1234)
   
.. seealso:: xb_run_register_


.. _xb_write_sh_scripts:

xb_write_sh_scripts
~~~~~~~~~~~~~~~~~~~

.. code:: text

   XB_WRITE_SH_SCRIPTS  Writes SH scripts to run applications on H4 cluster using MPI
   
   Writes SH scripts to run applications on H4 cluster. Optionally
   includes statements to run applications using MPI.
   
   Preferences:
   mpitype   = Type of MPI application (MPICH2/OpenMPI)
   
               Preferences overwrite default options (not explicitly
               defined options) and can be set and retrieved using the
               xb_setpref and xb_getpref functions.
   
   Syntax:
   fname = xb_write_sh_scripts(lpath, rpath, varargin)
   
   Input:
   lpath     = Local path to store scripts
   rpath     = Path to store scripts seen from h$ cluster
   varargin  = name:       Name of the run
               binary:     Binary to use
               nodes:      Number of nodes to use (1 = no MPI)
               mpitype:    Type of MPI application (MPICH2/OpenMPI)
   
   Output:
   fname     = Name of start script
   
   Example:
   fname = xb_write_sh_scripts(lpath, rpath, 'binary', 'bin/xbeach')
   fname = xb_write_sh_scripts(lpath, rpath, 'binary', 'bin/xbeach', 'nodes', 4)
   
.. seealso:: xb_run_remote_


.. _sec-xb_testbed:

xb_testbed
----------

.. _xb_testbed_check:

xb_testbed_check
~~~~~~~~~~~~~~~~

.. code:: text

   XB_TESTBED_CHECK  Checks if running as part of the XBeach testbed
   
   Checks if running as part of the XBeach testbed by checking if
   preferences are set.
   
   Syntax:
   istestbed = xb_testbed_check()
   
   Input:
   none
   
   Output:
   istestbed = Boolean indicating if testbed is running
   
   Example:
   if xb_testbed_check; disp('YES!'); end;
   
.. seealso:: xb_check_


.. _xb_testbed_getpref:

xb_testbed_getpref
~~~~~~~~~~~~~~~~~~

.. code:: text

   XB_TESTBED_GETPREF  Read testbed preferences
   
   Read testbed preferences
   
   Syntax:
   pref = xb_testbed_getpref()
   
   Input:
   none
   
   Output:
   pref    = Struct with testbed preferences
   
   Example:
   pref = xb_testbed_getpref
   
.. seealso:: xb_testbed_check_


.. _xb_testbed_includes:

xb_testbed_includes
~~~~~~~~~~~~~~~~~~~

.. code:: text

   XB_TESTBED_INCLUDES  Dummy function to automatically include OET functions in XBeach toolbox release
   
   Dummy function to automatically include OET functions in XBeach toolbox
   release. These functions can be used in the XBeach testbed, for
   Syntax:
   varargout = xb_testbed_includes(varargin)
   
   Input:
   varargin  = none
   
   Output:
   varargout = none
   Example:
   

.. _xb_testbed_loadskill:

xb_testbed_loadskill
~~~~~~~~~~~~~~~~~~~~

.. code:: text

   XB_TESTBED_LOADSKILL  Loads testbed skill history of specific variable
   
   Loads testbed skill history of specific variable
   
   Syntax:
   s = xb_testbed_loadskill(var)
   
   Input:
   var       = Variable name
   
   Output:
   s         = Struct with skill history
   
   Example:
   s = xb_testbed_loadskill('zb')
   
.. seealso:: xb_testbed_storeskill_


.. _xb_testbed_plotskill:

xb_testbed_plotskill
~~~~~~~~~~~~~~~~~~~~

.. code:: text

   XB_TESTBED_PLOTSKILL  Plot testbed skill history
   
   Plot testbed skill history
   
   Syntax:
   xb_testbed_plotskill(var)
   
   Input:
   var       = Variable to plot
   varargin  = r2:         Boolean indicating whether r^2 score should be
                           plotted
               sci:        Boolean indicating whether Sci score should be
                           plotted
               relbias:    Boolean indicating whether realtive bias should
                           be plotted
               bss:        Boolean indicating whether Brier Skill Score
                           should be plotted
   
   Output:
   none
   
   Example:
   xb_testbed_plotskill(var)
   
.. seealso:: xb_testbed_loadskill_, xb_testbed_storeskill_


.. _xb_testbed_storeskill:

xb_testbed_storeskill
~~~~~~~~~~~~~~~~~~~~~

.. code:: text

   XB_TESTBED_STORESKILL  Updates testbed skill history for specific variable
   
   Updates testbed skill history for specific variable
   
   Syntax:
   xb_testbed_storeskill(var, r2, sci, relbias, bss)
   
   Input:
   var       = Variable name
   r2        = R^2 score
   sci       = Sci score
   relbias   = Relative bias
   bss       = Brier Skill score
   
   Output:
   none
   
   Example:
   xb_testbed_storeskill('zb', r2, sci, relbias, bss)
   
.. seealso:: xb_testbed_loadskill_


.. _sec-xb_visualise:

xb_visualise
------------

.. _xb_get_handles:

xb_get_handles
~~~~~~~~~~~~~~

.. code:: text

   XB_GET_HANDLES  Return valid axes handles for creating plots
   
   Return valid axes handles for creating plots. By default, a single
   figure window is created with the requested number of axes as subplots.
   If a specific figure handle is provided, this figure window is used. If
   a list of specific axes handles are provided, these are checked for
   validity and appended with another figure window with subplots to
   arrive at the requested number of axes.
   
   Syntax:
   ax = xb_get_handles(n, varargin)
   
   Input:
   n         = Number of axes needed
   varargin  = handles:    Figure handle or list of axes handles
   
   Output:
   ax        = List with axes handles
   
   Example:
   ax = xb_get_handles(4)
   ax = xb_get_handles(4, 'handles', gcf)
   ax = xb_get_handles(2, 'handles', [ax1 ax2])
   
.. seealso:: xb_plot_hydro_, xb_plot_morpho_, xb_plot_sedtrans_


.. _xb_plot:

xb_plot
~~~~~~~

.. code:: text

   XB_PLOT  Basic visualisation tool for XBeach output structure
   
   Creates a basic interface to visualise contents of a XBeach output
   structure. The function is meant to be basic and provide first
   visualisations of model results.
   
   Syntax:
   varargout = xb_plot(xb, varargin)
   
   Input:
   xb        = XBeach structure array
   varargin  = none
   
   Output:
   none
   
   Example:
   xb_plot(xb)
   
.. seealso:: xb_read_dat_


.. _xb_plot_hydro:

xb_plot_hydro
~~~~~~~~~~~~~

.. code:: text

   XB_PLOT_HYDRO  Create uniform wave transformation plots
   
   Create uniform wave transformation plots from an XBeach hydrodynamics
   structure. Depending on the amount of information provided, different
   plots over the x-axis are plotted. Measurements are provided in a nx2
   matrix in which the first column is the x-axis and the second the data
   axis.
   
   Syntax:
   fh = xb_plot_hydro(xb, varargin)
   
   Input:
   xb        = XBeach output structure
   varargin  = handles:    Figure handle or list of axes handles
               zb:         Measured final profile
               Hrms_hf:    Measured high frequency wave height
               Hrms_lf:    Measured low frequency wave height
               Hrms_t:     Measured total wave height
               s:          Measured water level setup
               urms_hf     Measured high frequency oribtal velocity
               urms_lf     Measured low frequency oribtal velocity
               urms_t      Measured total oribtal velocity
               umean       Measured mean cross-shore flow velocity
               vmean       Measured mean longshore flow velocity
               units_dist: Units used for x- and z-axis
               units_vel:  Units used for secondary z-axis
               showall:    Show all data available instead of only show
                           data matched by measurements
   
   Output:
   none
   
   Example:
   xb_plot_hydro(xb)
   xb_plot_hydro(xb, 'zb', zb, 'Hrms_t', H)
   
.. seealso:: xb_plot_profile_, xb_plot_morpho_, xb_get_hydro_


.. _xb_plot_hydro2:

xb_plot_hydro2
~~~~~~~~~~~~~~

.. code:: text

   XB_PLOT_HYDRO2  Create uniform wave transformation plots (advanced)
   
   Create uniform wave transformation plots from an XBeach hydrodynamics
   structure. Depending on the amount of information provided, different
   plots over the x-axis are plotted. Measurements are provided in a nx2
   matrix in which the first column is the x-axis and the second the data
   axis. This function plots information additional to the xb_plot_hydro
   function and can be seen as the advanced part of the hydrodynamic
   analysis.
   
   Syntax:
   fh = xb_plot_hydro2(xb, varargin)
   
   Input:
   xb        = XBeach output structure
   varargin  = handles:    Figure handle or list of axes handles
               urms_hf     Measured high frequency oribtal velocity
               urms_lf     Measured low frequency oribtal velocity
               urms_t      Measured total oribtal velocity
               uavg        Measured mean undertow velocity
               rho         Measured correlation between short wave
                           variance and long wave surface elevation
               SK          Measured wave skewness
               AS          Measured wave asymmetry
               B           Measured wave nonlinearity
               beta        Measured wave phase
               units_dist: Units used for x- and z-axis
               units_corr: Units used for correlation z-axis
               units_skas: Units used for wave shape z-axis
               units_b:    Units used for wave nonlinearity z-axis
               units_vel:  Units used for flow velocity z-axis
               showall:    Show all data available instead of only show
                           data matched by measurements
   
   Output:
   none
   
   Example:
   xb_plot_hydro2(xb)
   xb_plot_hydro2(xb, 'zb', zb, 'Hrms_t', H)
   
.. seealso:: xb_plot_hydro_, xb_plot_profile_, xb_plot_morpho_, xb_get_hydro_


.. _xb_plot_klopman:

xb_plot_klopman
~~~~~~~~~~~~~~~

.. code:: text

   XB_PLOT_KLOPMAN  Plots LF spectrum based on HF spectrum
   
   Implementation of the Klopman and Dingemans (2000) paper based on the
   original secspecorg function made by Ap van Dongeren in 2001. The
   function computes a low-frequency spectrum from a JONSWAP spectrum.
   
   Syntax:
   varargout = xb_plot_klopman(varargin)
   
   Input:
   varargin  = Hm0:        significant wave height
               fp:         peak frequency
               df:         frequency resolution
               gamma:      peak-enhancement factor
               h:          water depth
               s:          directional spreading
               K:          number of wave components
               ntheta:     number of directions
               g:          gravitational acceleration
               plot:       boolean flag to enable plotting
   
   Output:
   f1        = LF frequencies
   S1        = LF energy density
   f2        = HF frequencies
   S2        = HF energy density
   
   Example:
   xb_plot_klopman;
   xb_plot_klopman('Hm0', 9, 'fp', 14);
   
.. seealso:: xb_plot_spectrum_


.. _xb_plot_morpho:

xb_plot_morpho
~~~~~~~~~~~~~~

.. code:: text

   XB_PLOT_MORPHO  Create uniform morphology plots
   
   Create uniform morphology plots from an XBeach morphology
   structure. Depending on the amount of information provided, different
   plots over the x-axis and time are plotted. Measurements are provided
   in a nx2 matrix in which the first column is the x-axis and the second
   the data axis.
   
   Syntax:
   fh = xb_plot_morpho(xb, varargin)
   
   Input:
   xb        = XBeach output structure
   varargin  = handles:    Figure handle or list of axes handles
               dz:         Measured bed level change
               sed:        Measured sedimentation volume in time
               ero:        Measured erosion volume in time
               R:          Measured retreat distance
               Q:          Measured separation point between erosion and
                           accretion
               P:          Measured delimitation of active zone
               units_dist: Units used for distance axes
               units_vol:  Units used for volume axis
               units_time: Units used for time axis
               showall:    Show all data available instead of only show
                           data matched by measurements
   
   Output:
   none
   
   Example:
   xb_plot_morpho(xb)
   xb_plot_morpho(xb, 'dz', dz, 'sed', sed)
   
.. seealso:: xb_plot_profile_, xb_plot_hydro_, xb_get_morpho_


.. _xb_plot_mpi:

xb_plot_mpi
~~~~~~~~~~~

.. code:: text

   XB_PLOT_MPI  One line description goes here.
   
   More detailed description goes here.
   
   Syntax:
   varargout = xb_plot_mpi(varargin)
   
   Input:
   varargin  =
   
   Output:
   varargout =
   
   Example:
   xb_plot_mpi
   
   See also

.. _xb_plot_profile:

xb_plot_profile
~~~~~~~~~~~~~~~

.. code:: text

   XB_PLOT_PROFILE  Create uniform profile plots
   
   Create uniform profile plots from XBeach output structure and other
   sources. Standardized coloring is used per source. Automatically
   computes Brier Skill Scores and makes sure the dune is located right on
   the figure. Additional profiles are supplied in the form of matrices
   from which the first column contains x-coordinates and successive
   columns contain z-values.
   
   Syntax:
   xb_plot_profile(xb, varargin)
   
   Input:
   xb       =  XBeach output structure
   varargin =  handle:         Figure or axes handle
               measured:       Measured post storm profile
               xbeach:         Profile computed by another version of
                               XBeach
               durosta:        Profile computed by DurosTA
               duros:          Profile computed by DUROS
               duros_p:        Profile computed by DUROS+
               duros_pp:       Profile computed by D++
               nonerodible:    Non-erodible layer
               units:          Units used for x- and z-axis
               BSS:            Boolean indicating whether BSS should be
                               included
               flip:           Boolean indicated whether figure should be
                               flipped in case dune is located left
   
   Output:
   none
   
   Example:
   xb_plot_profile(xb, 'initial', profile0, 'measured', profile1)
   
.. seealso:: xb_plot_hydro_, xb_plot_morpho_


.. _xb_plot_sedtrans:

xb_plot_sedtrans
~~~~~~~~~~~~~~~~

.. code:: text

   XB_PLOT_SEDTRANS  Create uniform wave sediment transport plots
   
   Create uniform sediment transport transformation plots from an XBeach
   sediment transport structure. Depending on the amount of information
   provided, different plots over the x-axis are plotted. Measurements are
   provided in a nx2 matrix in which the first column is the x-axis and
   the second the data axis.
   
   Syntax:
   fh = xb_plot_sedtrans(xb, varargin)
   
   Input:
   xb        = XBeach output structure
   varargin  = handles:    Figure handle or list of axes handles
               c:          Measured sediment concentration
               k_wavg:     Measured wave averaged turbulence
               k_bavg:     Measured bore averaged turbulence
               k_bavg_nb:  Measured bore averaged near-bed tubulence
               S_dz:       Measured sediment transport from bed level changes
               S_av:       Measured sediment transport from avalancing
               S_s:        Measured suspended sediment transports
               S_b:        Measured bedload sediment transports
               S_lf:       Measured low frequency related sediment
                           transports
               S_ut:       Measured short wave undertow related sediment
                           transports
               S_as:       Measured sediment transports related to wave
                           asymmetry
               units_dist: Units used for x- and z-axis
               units_vol:  Units used for transport z-axis
               showall:    Show all data available instead of only show
                           data matched by measurements
   
   Output:
   none
   
   Example:
   xb_plot_sedtrans(xb)
   xb_plot_sedtrans(xb, 'zb', zb, 'Hrms_t', H)
   
.. seealso:: xb_plot_profile_, xb_plot_morpho_, xb_plot_hydro_, xb_get_sedtrans_


.. _xb_plot_spectrum:

xb_plot_spectrum
~~~~~~~~~~~~~~~~

.. code:: text

   XB_PLOT_SPECTRUM  Create uniform spectrum plots
   
   Create uniform spectrum plots from XBeach spectrum structure. Creates a
   subplot per location.
   
   Syntax:
   fh = xb_plot_spectrum(xb, varargin)
   
   Input:
   xb        = XBeach spectrum structure
   varargin  = measured:   Measured spectra
               units:      Units used for x- and z-axis
               units2:     Units used for secondary z-axis
   
   Output:
   fh        = Figure handle
   
   Example:
   fh = xb_plot_spectrum(xb)
   
.. seealso:: xb_plot_hydro_, xb_plot_morpho_, xb_get_spectrum_


.. _xb_select_crop:

xb_select_crop
~~~~~~~~~~~~~~

.. code:: text

   function to draw a rubberband line and return the start and end points
   Usage: [xi,yi] = select_oblique_rectangle;     uses current axes
   or    [xi,yi] = select_oblique_rectangle(h);  uses axes refered to by handle, h

.. _xb_view:

xb_view
~~~~~~~

.. code:: text

   XB_VIEW  Generic GUI to visualize XBeach input and output
   
   Graphical User Interface to visualize input and output. The viewer is
   *complementary* to the Delft3D QuickPlot (QP) viewer. Most options you
   find in this viewer are lacking in QP and vise versa.
   
   This viewer accepts XBeach input, output and run structures or a path
   string to an XBeach output directory. In case no input is given, the
   current directory is assumed to contain XBeach output. It is also
   possible to provide a cell array with different datasources. These
   sources can be a mix of types (structures and paths), but need to have
   the same spatial grid and equal timestep.
   
   The most important options are:
       * 1D, 2D and 3D plots
       * Multiple variable selection
       * Time sliders and animation
       * Time difference plots
       * Time comparison plots
       * Multiple datasource comparison plots
       * Transect view for 2D models with transect slider
       * Fixation and alignment of caxis
       * Measurement tool
   
   xb_view;
   xb_view(xbi);
   xb_view(xbr);
   xb_view(xbo);
   xb_view('path/to/output/dir');
   xb_view(pwd, 'width', 1024, 'height', 768, 'modal', true);
   
   Syntax:
   xb_view(data, varargin)
   
   Input:
   data      = XBeach structure (input, output or run) or path to XBeach
   Output:
   varargin  = width:  Width of window at startup (inf = screensize)
               height: Height of window at startup (inf = screensize)
               model:  Boolean indicating modal state
   
   Output:
   none
   
   Example:
   xbi = xb_generate_model;
   xbr = xb_run(xbi, 'path', 'path/to/output/dir');
   xbo = xb_read_output('path/to/output/dir');
   
.. seealso:: xb_plot_profile_


.. _xb_visualize_modelsetup:

xb_visualize_modelsetup
~~~~~~~~~~~~~~~~~~~~~~~

.. code:: text

   XB_VISUALIZE_MODELSETUP  One line description goes here.
   
   More detailed description goes here.
   
   Syntax:
   varargout = xb_visualize_modelsetup(varargin)
   
   Input:
   varargin  =
   
   Output:
   varargout =
   
   Example:
   xb_visualize_modelsetup
   
   See also

.. _xb_write_skilltable:

xb_write_skilltable
~~~~~~~~~~~~~~~~~~~

.. code:: text

   XB_WRITE_SKILLTABLE  Writes Latex table with skills to file
   
   Writes Latex table with skills to file
   
   Syntax:
   varargout = xb_write_skilltable(xb, measured, varargin)
   
   Input:
   xb        = XBeach output structure
   measured  = Cell array containing measurement data with in the first
               column the x-axis values and the second column z-axis
               values
   initial   = Cell array containing initial data with in the first
               column the x-axis values and the second column z-axis
               values
   varargin  = file:   Filename
               title:  Table title
               vars:   Cell array with variable names to plot
   
   Output:
   none
   
   Example:
   xb_write_skilltable(xb, measured, initial)
   
.. seealso:: xb_skill_


