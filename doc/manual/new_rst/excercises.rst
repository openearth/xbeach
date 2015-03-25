Excercises
==========

The hands-on exercises can be downloaded via subversion. Subversion is a
well-known version management system that allows you to always have the
most recent source code at hand. It also allows developers to commit
changes to the source code, without interfering with other developers.
In order to use Subversion, you will need a Subversion client. A
well-known client for Windows is Tortoise. If you have registered, you
can download the source code via the following URL:
*https://svn.oss.deltares.nl/repos/xbeach/Courses/DSD\_2014/Examples –
Basic*. For the tools like Quickplot and Quickin the Delft3D environment
is needed.

Dune erosion at Delfland, Netherlands (1D)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first case we will run is a relative simple 1D case. It concerns a
profile along the Dutch coast and the hydraulic boundary conditions are
based on the 1953 storm surge that caused substantial flooding in the
Netherlands.

You can work on the following assignments:

#. Go to the folder “Examples\\DelflandStorm” and double click the file
   “run\_model.bat”. The simulation will start. The model will run for a
   few minutes, but in the meantime you can already work on question 2
   to 5.

#. Open params.txt in which you specify the model input files and
   settings. Check the number of grid-points in x-direction (keyword:
   *nx*) and y-direction (keyword: *ny*). Check the filenames in which
   you specify the wave conditions (keyword: *bcfile*) and the storm
   surge level (SSL) (keyword: *zs0file*).

#. Do the wave conditions change during the simulation? What is/are the
   wave height(s) and wave period(s) applied in the simulation?

#. Does the storm surge level change during the simulation? What is the
   maximum surge height in the simulation. Surge height is defined with
   respect to the mean sea level (MSL)?

#. What is the simulation time (keyword: *tstop*)? Do we apply a
   morphological acceleration factor (keyword: *morfac*)? What variables
   are stored as output and with what time interval? How much
   hydrodynamic time is simulated?

#. Probably the simulation has finished. When you start the model, it
   generates a file named XBlog.txt. Open this file and check what is
   stored in the file. What was the total simulation time?

#. To check out the simulation results we make use of the Quickplot tool
   (A brief tutorial is attached to this document). You can start
   Quickplot via the Delft3D environment we installed (Start Programs
   Deltares Delft3D Delft3D). In the Delft 3D menu choose Utilities
   Quickplot. Choose Files of type “NetCDF files and GRIB files” and
   open “xboutput.nc” in the simulation folder.

#. Use the Quickplot tutorial and try to make an animation in which you
   plot short wave height (H), water level (including long wave
   variations, zs) and bed level (zb) as function of time.

#. Plot the offshore water level as function of time. Also open the file
   “tide.tek” (Tekal data files format), which contains the imposed
   surge level. Did the model correctly simulate the imposed surge
   level?

#. Copy all model files to a new folder named “superfast”. Edit
   params.txt and set ny=0 (instead of ny=2), and run the model. What is
   the simulation time compare to the original simulation?

#. Compare simulation results for the “superfast” and “default”
   simulation. Are these the same? What option will you use in the
   future?

Nourishment scenarios near Kijkduin, Holland (1D) 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This case concerns the exploration of a nourishment strategy near
Kijkduin along the Holland coast in the Netherlands. At this location a
mega nourishment of 21 Mm\ :math:`{}^{3}` named the Sand Engine was
constructed. In this case we will explore to what extent nourishments
can reduce the (dune and beach) erosion during a storm event.

You can work on the following assignments:

#. Go to the folder “Examples\\Nourishment case” and double click on the
   file “runall.bat”. This batch file will run three simulations
   sequentially in which the profile configuration varies and
   corresponds with the undisturbed profile (folder reference), a
   shoreface nourishment (folder shoreface) and a beach nourishment
   (folder beach) respectively. Each model will run for a few minutes.
   While running you can already answer question 2 to 6.

#. For the reference case open the params.txt in which you specify model
   input files and settings. Check the number of grid-points in
   x-direction (keyword: *nx*) and y-direction (keyword: *ny*). How many
   directional wave bins are defined and what is their width (keywords:
   *thetamin*, *thetamax*, *dtheta*).

#. Do the wave conditions change during the simulation? What is/are the
   wave height(s) and wave period(s) applied in the simulation?

#. Does the storm surge level change during the simulation? What is the
   maximum surge height in the simulation. Surge height is defined with
   respect to the mean sea level (MSL)?

#. What is the simulation time (keyword: *tstop*)? Do we apply a
   morphological acceleration factor (keyword: *morfac*)? What variables
   are stored as output and with what time interval? How much
   hydrodynamic time is simulated?

#. Probably the simulation has finished. When you start the model, it
   generates a file named XBlog.txt. Open this file and check what is
   stored in the file. What was the total simulation time?

#. Inspect the initial bathymetries of each simulation with QUICKPLOT.
   Choose Files of type “NetCDF files and GRIB files” and open
   “xboutput.nc” in the simulation folder).

   #. At what cross-shore position were the shoreface nourishment and
      beach nourishment placed?

   #. What is the (average) thickness of the nourishments?

   #. Is the volume of the nourishments comparable?

   #. Plot the reference profile with markers; does the grid resolution
      vary in cross-shore direction?

#. Use the Quickplot tutorial and try to make an animation in which you
   plot short wave height (H), water level (including long wave
   variations, zs) and bed level (zb) as function of time.

#. Plot the offshore water level as function of time. Also open the file
   “tide.tek” (Tekal data files format), which contains the imposed
   surge level. Did the model correctly simulate the imposed surge
   level?

#. Inspect the final bathymetries of each simulation.

   #. What is the dune face retreat in the three simulations you have
      carried out?

   #. Where does the eroded sediment form the dunes deposit?

   #. What nourishment type is most effective in reducing the impact of
      a storm and do you have an explanation for this?

#. In the folder “banquette” you find a final simulation in which a
   special beach nourishment type is evaluated named a banquette. This
   beach nourishment has a highly elevated flat area that connects to
   the dune foot on which beach restaurants can be build.

   #. Run the model and compare in Quickplot the banquette design with
      the beach nourishment design we have evaluated before. Do you
      expect more or less erosion?

   #. Check your hypothesis by comparing the final profile of the
      banquette simulation to the other simulations.

   #. What would be your approach to further reduce beach and dune
      erosion?

Overwash at Santa Rosa Island , USA (2DH) 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This case concerns overwash at Santa Rosa island in the Gulf of Mexico
during hurricane Ivan in 2004.

You can work on the following assignments.

#. For the reference case open the params.txt in which you specify
   model input files and settings. Check the number of grid-points in
   x-direction (keyword: *nx*) and y-direction (keyword: *ny*). How many
   directional wave bins are defined and what is their width (keywords:
   *thetamin*, *thetamax*, *dtheta*).

#. In this simulation the grid is specified in Delft3D format. Open
   Quickin in the Delft 3D menu (Grid Quickin) and use the brief
   tutorial to read in the grid and bathymetry. Does the grid resolution
   vary in cross-shore direction? And in longshore direction? What are
   the minimum dx and dy? Why can the grid be coarse offshore?

#. How many wave conditions do we apply in this simulation? What is the
   offshore mean wave direction? Does the surge level change in the
   simulation?

#. What is the simulation time (hydrodynamic and morphologic)?

#. Inspect the model results and make an animation of the short wave
   height (H) and the water levels (including long wave, zs). Describe
   what is happening.

#. For the water levels set the color limits manual between -0.5 and
   3.5.

#. Make an animation of cumulative sedimentation/erosion. Describe what
   is happening.

#. For the sedimentation/erosion set the color limits manual between -3
   and 3

#. Look at the mean flow field. Plot the flow field in colored vectors.
   Where are the flow velocities highest and what is the direction of
   the flow (cross-shore or longshore)? Is there (also) a longshore
   current present and what is its intensity?

If you have time left feel free to:

#. Narrow or broaden the imposed spectrum by changing the parameter
   directional spreading (*s*) in ‘jonswap.inp’ (you could for example
   set s = 100 and s = 2 respectively). Make animations of the
   instantaneous short wave height to see what is happening to the size
   of the wave groups.

#. Design a nourishment in Quickin to reduce the impact of the storm on
   Santa Rosa Island. Change the depth file in params.txt to make a
   simulation with the updated bathymetry.

Yanchep perched beach and natural breakwater (2DH) 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This case is an example of a beach 60km north of Perth most commonly
known as Yanchep lagoon. Many beaches in WA like Yanchep are fronted by
shallow reef and here we are investigating the effects of the reef on
the morphodynamics.

You can work on the following assignments:

#. Go to the folder “Examples\\YanchepBeach” and double click the
   file “run\_model.bat”. The simulation will start (and will run about
   15 minutes).

#. Meanwhile, inspect the bathymetry file and the structure file (using
   Quickin). What is the depth in the lagoon? Is the reef enclosing the
   lagoon below or above the model initial water level? What is the wave
   height at the boundary condition?

#. Use Quickplot and try to make an animation in which you plot short
   wave height (H), water level (including long wave variations) (zs)
   and Eulerian velocities (ue and ve) as function of time. What happens
   in the lagoon?

#. Use Quickplot and try to make an animation of cumulative
   sedimentation/erosion. What happens in the lagoon?

#. How is the lagoon affected by the mean water level? Increase or
   decrease the mean water level condition (‘tide.tx’), run the model
   again (maybe for a shorter time by reducing keyword: *tstop*). How
   are the circulation and sediment transport affected?

#. What would happen if the lagoon was open at the southern end? Open
   the structure file (keyword: ne\_layer=’reef.dep’) with the Quickin
   tool and modify it to allow the southern end of the lagoon to be
   eroded. Modify the param.txt file to use this new structure file and
   run the model. Alternatively, remove the reef from the bathymetry and
   rerun the model without the structure file, by setting the keyword
   *struct=0*.

If you still have time;

#. Reefs are very rough what happens in the model when the friction is
   increased? Reduce the Chezy roughness and increase the value of
   *f\ :math:`{}_{w}`*. Rerun the model what do you observe?

#. Is wave/current interaction (keyword: *wci=1*) switched on? Rerun the
   model with the wave/current switch on/off. Compare the output with
   model you ran previously. How much effect do you see on the
   morphology?
