Teledyne RDI Instruments
************************

Currents data
=============

WorkHorse instruments
---------------------

#. Open Velocity software.

#. Click Options -> ensure metric units are selected.

#. Click Tools -> Check PD0 Data File. Ensure *Output Date/Time* is checked. This will generate a .log file and a .sensors.txt file.

#. Click Open -> select your raw .000 file.

#. Click the gears icon in upper-right (Processing settings).

   * Ensure coordinate system is Earth. Keep magnetic variation at 0. Under Range to Boundary, select None. Under Averaging, uncheck Average data.

#. Click Export -> Export to Matlab to save a .mat file.

#. Click Export -> Export to ASCII.

   * Ensure all cells and times are selected.

   * Under Output: click Descriptions, Sample number, Date and time, and Sensors data. Save as .txt file.

Now you will have a directory of files that looks like:

.. list-table::

  * - PCT20000.000
    - original raw data file
  * - PCT20000.log
    - file exported by Check PD0 Data File
  * - PCT20000.sensors.txt
    - file exported by Check PD0 Data File
  * - PCT20000.000.mat
    - file exported by Export to Matlab
  * - PCT20000.000.txt
    - file exported by Export to ASCII containing pressure data

Use :doc:`runots.py </runots>` to process using the two :doc:`configuration files </config>`.

For basefile, in the above example, you would use PCT20000.

Velocity will split files at 25,000 ensembles, so you may end up with multiple files. As long as they all have the same base filename (e.g., PCT20000.000.mat, PCT20000.001.mat, PCT20000.000.txt, PCT20000.001.txt), stglib can handle these files. 



Sentinel V instruments
----------------------

Coming soon.

Waves processing
================

Coming soon.
