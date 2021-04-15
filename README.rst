.. image:: https://raw.githubusercontent.com/sammmlow/ORBITM/master/gui/orbm_logo_large.png

.. |docs| image:: https://img.shields.io/badge/docs-latest-brightgreen.svg?style=flat-square
   :target: https://orbitm.readthedocs.io/en/latest/

.. |license| image:: https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square
   :target: https://github.com/sammmlow/ORBITM/blob/master/LICENSE
   
.. |orcid| image:: https://img.shields.io/badge/ID-0000--0002--1911--701X-a6ce39.svg
   :target: https://orcid.org/0000-0002-1911-701X/
   
.. |linkedin| image:: https://img.shields.io/badge/LinkedIn-sammmlow-blue.svg
   :target: https://www.linkedin.com/in/sammmlow

:Project: Orbit.M
:Github: https://github.com/sammmlow/ORBITM
:Documents: https://orbitm.readthedocs.io/en/latest/
:Version: 1.0 (Stable)

|docs| |license|

:Author: Samuel Y. W. Low

|linkedin| |orcid|



Welcome to Orbit.M!
-------------------

Orbit.M is an open-source and free orbit maintenance simulator and propulsion sizing tool, recommended for use with near-circular low Earth orbits up to 1,000 km altitude, on Python.

It stands for the **Orbit Maintenance and Propulsion Sizing Tool**, and it comes with its own built-in orbit decay model and maintenance simulator. Alternatively, Orbit.M can also be used to interface with and automate orbit maintenance simulations using AGI's Systems Tool Kit (STK) as an alternative simulator, through the STK Integration Object Model libraries.

A valid license for STK 10 and 11, with Astrogator and Integration modules, is needed for interfacing with STK.

The objective of Orbit.M is to allow for a quick sizing of low Earth orbit (LEO) mission lifetimes, sized against propulsion units of the user's choosing. The user can enter the parameters of their intended mission, the orbital elements, the spacecraft characteristics, and Orbit.M would compute the ΔV necessary to counteract drag forces throughout the mission, while sizing it against your choices of thrusters (depending on your inputs in the **"thruster_shortlist.txt"** file).



First Steps in Using Orbit.M
----------------------------

First things first, check that you have these Python libraries: **TKinter, NumPy, Matplotlib**, and among other standard libraries such as **os, datetime, comtypes etc**. If you do not have Python, I recommend using the `Anaconda installer <https://www.anaconda.com/>`_ with Spyder as your IDE, as it comes with all the installed packages by default.

Now, the software can be started by running the Python file **orbitm.py** in the main directory (equivalent to the directory you see on the master branch on Orbit.M's github page). You should see a GUI, like below, pop up:

.. image:: https://raw.githubusercontent.com/sammmlow/ORBITM/master/gui/orbm_screenshot.png

You can fill in your spacecraft and orbit parameters through the GUI, check that you have your desired thruster specified in the "thruster_shortlist.txt" file on the main Orbit.M directory.

Now, you can hit **Run ORBITM**, and it will output a ΔV report as a text file, and plot the altitude profile and the locus of thruster requirements that satisfy your mission (in terms of the required fuel mass against Isp). Example plots are shown below, with plots in blue ran under the fast algorithm (Sam's), and plots in orange ran under the more precise orbit model in STK.

.. image:: https://raw.githubusercontent.com/sammmlow/ORBITM/master/gui/orbm_outp_550km.png

The above results were run for a 200kg satellite at 550km circular orbit. Note that the blue and orange plots in the thruster requirements profile almost overlap.

For full documentation, please refer to the `Orbit.M Read-The-Docs. <https://orbitm.readthedocs.io/en/latest/>`_



Final Notes
-----------

This project is free, and open-source. If you would like to contribute to this project, add in new features, or enhance existing atmospheric models etc, please free free to branch the master repository and collaborate with me. I am also hoping to interface Orbit.M's use with other astrodynamics softwares like GMAT, and improve on existing features such as improving the accuracy of the atmospheric density models etc.

If you had felt that OrbitM was useful in your research, please do give the due credit and cite my paper or this project!

Low, S. Y. W., &; Chia, Y. X. (2018). “Assessment of Orbit Maintenance Strategies for Small Satellites”, 32nd Annual AIAA/USU Conference on Small Satellites, Logan, Utah, Utah State University, USA.

Contact
-------

If you have any other queries feel free to reach out to me at:

sammmlow@gmail.com

|linkedin| |orcid|

*Last Modified on 03-04-2021*

