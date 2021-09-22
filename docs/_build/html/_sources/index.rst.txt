.. 
   ###########################################################################
   ###########################################################################
   ##                                                                       ##
   ##     _____ ___  ___  ___  _____      __  __                            ##
   ##    |  _  | _ \| _ \|_ _||_   _|    |  \/  |                           ##
   ##    | |_| |   <| _ < | |   | |   _  | \  / |                           ##
   ##    |_____|_|\_|___/|___|  |_|  |_| |_|\/|_|                           ##
   ##                                                     v 1.1             ##
   ##                                                                       ##
   ###########################################################################
   ###########################################################################
   
.. |docs| image:: https://img.shields.io/badge/docs-latest-brightgreen.svg?style=flat-square
   :target: https://orbitm.readthedocs.io/en/latest/

.. |license| image:: https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square
   :target: https://github.com/sammmlow/ORBITM/blob/master/LICENSE
   
.. |orcid| image:: https://img.shields.io/badge/ID-0000--0002--1911--701X-a6ce39.svg
   :target: https://orcid.org/0000-0002-1911-701X/
   
.. |linkedin| image:: https://img.shields.io/badge/LinkedIn-sammmlow-blue.svg
   :target: https://www.linkedin.com/in/sammmlow

.. image:: /_static/orbitm_logo.png

|

:Github: https://github.com/sammmlow/ORBITM
:Documents: https://orbitm.readthedocs.io/en/latest/
:Version: 1.1 (Latest)
:Author: Samuel Y. W. Low

|docs| |license| |linkedin| |orcid|

ORBITM
======

ORBITM stands for the **Orbit Maintenance and Propulsion Sizing Tool**, and it is an open-source, easy-to-use, free-ware orbit maintenance simulator and propulsion sizing tool, for anyone with a Python 3 installation. The software comes with its own built-in orbit decay model and maintenance simulator, a life-time estimator, and it is also interface-able with AGI's Systems Tool Kit (STK) as an alternative mode of simulation.

The objective of ORBITM is to allow for a quick sizing of low Earth orbit (LEO) mission lifetime, while sizing the mission against propulsion units of the user's choosing. The user may input into the UI their orbital parameters, spacecraft properties (mass, area, drag parameters etc), and the software will compute your desired ΔV necessary for the mission, track the altitude profile over time, and plot the thruster-fuel profile that satisfies the orbit maintenance needs of your mission.

Document tree for ORBITM is listed out below.

.. toctree::
   :maxdepth: 1
   :caption: Getting Started
   
   docs_overview.rst
   docs_install.rst
   docs_firststep.rst
   docs_interpreting.rst

.. toctree::
   :maxdepth: 1
   :caption: Orbit Maintenance
   
   docs_orbmaint.rst
   docs_functions.rst

|

If ORBITM has made your life somewhat easier, please do cite the work behind it!

`Low, S. Y. W., &; Chia, Y. X. (2018). “Assessment of Orbit Maintenance Strategies for Small Satellites”, 32nd Annual AIAA/USU Conference on Small Satellites, Logan, Utah, Utah State University, USA. <https://digitalcommons.usu.edu/smallsat/2018/all2018/364/>`_


For bugs, raise the issues in the `GitHub repository <https://github.com/sammmlow/ORBITM/issues>`_.

For collaborations, reach out to me: sammmlow@gmail.com (Samuel Y. W. Low)

|linkedin| |orcid|

The project is licensed under the MIT license.


.. toctree::
   :maxdepth: 2
   :caption: Contents:
