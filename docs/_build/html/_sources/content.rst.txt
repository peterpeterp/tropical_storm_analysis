=======
Modules
=======

---------------------------------------
Tropical cyclone detection and tracking
---------------------------------------

.. currentmodule:: TC_detection

This part of the documentation shows the full API reference of all public
functions. 

.. autoclass:: tc_tracks

~~~~~~~~~~~~~~~~~~~~~~~~~
Detecting storm snapshots
~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: tc_tracks.detect_contours
.. automethod:: tc_tracks.detect_knutson2007
.. automethod:: tc_tracks.find_closed_contours

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Combine storm snapshots to tracks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: tc_tracks.combine_tracks

~~~~~~~~~~~~~~~~~~
Plotting functions
~~~~~~~~~~~~~~~~~~
.. automethod:: tc_tracks.init_map
.. automethod:: tc_tracks.plot_on_map
.. automethod:: tc_tracks.plot_track_path
.. automethod:: tc_tracks.plot_all_tracks
.. automethod:: tc_tracks.plot_detect_summary
.. automethod:: tc_tracks.plot_surrounding

~~~~~~~~~~~~~~~~~~~~
Supporting functions
~~~~~~~~~~~~~~~~~~~~
.. automethod:: tc_tracks.get_box
.. automethod:: tc_tracks.area_around
.. automethod:: tc_tracks.circle_around
.. automethod:: tc_tracks.degree_to_step
.. automethod:: tc_tracks.tc_cat

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Analysis using observed track data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: tc_tracks.init_obs_tcs
.. automethod:: tc_tracks.obs_track_info
