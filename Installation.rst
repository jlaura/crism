.. _installation:

********************************
Installing CSAS
********************************

.. toctree::
    :maxdepth: 2
    
CSAS Dependencies
====================

CSAS requires Python and depends on several other freely available Python
modules. Prior to installing CSAS, you should make sure its dependencies are met.

.. list-table:: CSAS Dependencies
   :header-rows: 1
   :widths: 20, 25
   :class: center

   * - Dependency
     - Requirement
   * - `Python 2.5+ <http://www.python.org>`_
     - Required
   * - `NumPy <http://numpy.scipy.org/>`_
     - Required
   * - `GDAL <http://gdal.org//>`_
     - Required  
     
Note that CSAS does not work with Python 3.x.

Installing 
==========
CSAS is distributed as a stand alone script and therefore does not require installation.  This documentation is distributed locally, with CSAS.  Simply place the script in a convenient directory.

Should you wish to access CSAS from any directory on your machine, it is necessary to append your PYTHONPATH with the path to CSAS.  To do this, first not the absolute path to the directory where CSAS is installed.  For example: '/Users/Jay/CSAS'.  Then run the following from a command line (run cmd on a Windows machine / terminal on an OS X machine).

OS X
----
To allow CSAS to be called from any directory via the command line utilize the following.  To have these changes be persistent append them to ~/.profile or ~/.bashrc.::

	$ PYTHONPATH="/Users/Jay/CSAS:$PYTHONPATH"
	$ export PYTHONPATH
	
	where, /Users/Jay/CSAS is replaced with the path to your CSAS installation
	
Windows
------- 
To allow CSAS to be called from any directory via the command line utilize the following.  To have these changes be persistent append them to autoexec.bat.::

	$ set PYTHONPATH=%PYTHONPATH%;C:\CSASInstallDirectory
	
You can also append autoexec.bat graphically.  Run *msconfig* from the start menu.
