![Orbit.M - Orbit Maintenance Analysis for LEO in Python](https://raw.githubusercontent.com/sammmlow/ORBITM/master/gui/orbm_logo_large.png)

### What is Orbit.M?

Orbit.M is an open-source, easy-to-use, free-ware orbit maintenance simulator and propulsion sizing tool, meant for circular orbits, for anyone and anywhere.

It is a software with a graphical user interface where you can fill in your orbital parameters, satellite masses, satellite areas, drag parameters, and the software will compute your desired Delta-V necessary for the mission, sizing it against your choices of thrusters (based on the Isp). For beginners who are not familiar with the terms "Delta-V" and "Isp", the Wikipedia articles about them, below, do a fantastic job of explaining what they are.

[WikiPedia article on Delta-V.](https://en.wikipedia.org/wiki/Delta-v)

[WikiPedia article on Specific Impulse, or Isp.](https://en.wikipedia.org/wiki/Specific_impulse)

Once you input your orbit and satellite parameters in the GUI, and you have configured your desired thrusters that you want to use in your mission Delta-V sizing profile, the software will output a mission Delta-V profile against the Isp of your thruster choices, as well as an altitude and mean semi-major axis chart. There are examples showing these charts down below.

### How do I use Orbit.M specifically?

First things first, check that you have these Python libraries: **TKinter, NumPy, Matplotlib**, and among other standard libraries such as **os, datetime, comtypes etc**. If you do not have Python, I recommend using the [Anaconda installer](https://www.anaconda.com/), and running it on Spyder 4 simply because that's what I'm using right now. At the time of the latest commit and push, NumPy version is 1.18, and Matplotlib is 3.1.3. Older versions should work fine as ORBIT.M doesn't use any of the newer features (let me know if it does not).

If you would like to use the STK AstroGator module in STK10 or STK11 as a high precision orbit propagator and simulator instead of my orbit maintenance simulator, select the program depending on which STK license you have. **You will need a valid STK license for STK Integration and STK Astrogator in order to use this module in ORBITM.**

Now, the software can be started by running the Python file **orbitm.py** in the main directory (which is equivalent to the directory you see on the master branch on ORBITM's github page). You should see a GUI that looks like the one below, pop up:

![Orbit.M - GUI running on TKinter](https://raw.githubusercontent.com/sammmlow/ORBITM/master/gui/orbm_screencap.png)

You can fill in your spacecraft and orbit parameters through the GUI (quite self-explanatory), which actually updates it in a "config.txt" file. You can alternatively update it in the "config.txt" file manually, but it is not recommended since then the software can't check for errors (e.g. a negative drag coefficient was typed into the config file which would crash the program).

Next, before you run ORBITM, you should open up the "thruster_shortlist.txt" file on the main ORBITM directory, and fill in any thrusters you wish to size your missions against. Some thrusters I had previously been aware of are written inside this shortlist as an example. The purpose of this shortlist is to graphically compare its Isp and fuel capacity to the suitability of your mission later on.

Now, you can run ORBITM.

### How do I use Orbit.M specifically?

Below, we have ran three scenarios for ORBITM - each is a circular orbit at 63.4 degrees inclination, at 450km, 500km, and 550km mean altitude respectively.

The plots in blue correspond to "Sam's Decay Model", a quick-computed decay model for circular orbits which I derived and wrote as part of a paper, using some back-of-the-envelope physics [(Assessment of Orbit Maintenance Strategies for Small Satellites)](https://digitalcommons.usu.edu/smallsat/2018/all2018/364/). This computation is done in closed-form, and thus the solution to your orbit maintenance Delta-V needs are computed at a fraction of a second, as compared to running a full STK simulation. It also gives a decently accurate ball-park figure.

The plots in orange correspond to using STK 10's Astrogator, with a high precision orbit propagator, and it uses Astrogator's Automatic Sequences feature to boost the orbit up whenever it hits the minimum of the maintenance tolerance band. Astrogator's Target sequence was not used to compute the Delta-V. The Delta-V is computed offline using a first order Taylor expansion to the Vis-Visa equation. **Again, you will need a valid STK license for STK Integration and STK Astrogator in order to use this module in ORBITM.**

For the circular orbit @ 450km, with a tolerance band of 5km:

![Orbit.M - Results for a circular orbit @ 450km](https://raw.githubusercontent.com/sammmlow/ORBITM/master/gui/orbm_outp_450km.png)

For the circular orbit @ 500km, with a tolerance band of 5km:

![Orbit.M - Results for a circular orbit @ 500km](https://raw.githubusercontent.com/sammmlow/ORBITM/master/gui/orbm_outp_500km.png)

For the circular orbit @ 550km, with a tolerance band of 5km:

![Orbit.M - Results for a circular orbit @ 550km](https://raw.githubusercontent.com/sammmlow/ORBITM/master/gui/orbm_outp_550km.png)

Probably the most useful graph you'd need is the third one from the left - you can input your thruster specs, and compare the height of the bar chart against the desired fuel mass needed for the locus of all specific impulse values (Isp, units in seconds). For each thruster's Isp value, if the height of the bar chart (which represents the max fuel mass) exceeds the fuel mass requirements of your mission at that particular Isp, then that thruster works for you.

The program will also output a time-schedule for orbit maintenance in the main directory of ORBITM, as a text file "deltaV.txt". This text file gives you the time stamps of when each thrust should occur, which is primarily dependent on the decay computation (drag effects) and your chosen tolerance band (i.e. how far do you descend before you fire your thrusters).

**Do note that you can scale your fuel mass requirements in the line plots by increasing the Maintenance Mission Margin parameter in the GUI.** A value of 1.0 implies that the program will compute the fuel mass and Delta-V equivalently needed to counter 100% of the drag effects you encounter. If you would like to have a margin of 3x amount of fuel, then input 3.0 as your mission margin.

### What is the difference between using Sam's and STK's simulator?:

My own model, as aforementioned, is a simplified decay model that solves for the decay at each time step in closed-form. Thus, no true orbit propagation is computed (i.e. no computation of state vectors, no Runge-Kutta integration and all that). In fact, only the orbit radial position scalar is computed by solving for Kepler's equation for the two-body case in each time step. This was done because we actually only need the altitude (and not the full set of state vectors) in order to estimate the atmospheric density. This saves us a lot of computational time using my model, albeit with some accuracy loss. The atmospheric density model used is the US Standard Atmosphere 1976 Model [(see PDF here)](https://github.com/sammmlow/ORBITM/blob/master/docs/USSA1976.pdf). A rough flow chart of using Sam's Model works is as follows:

![Orbit.M - Flow Chart for Sam's Decay Computer](https://raw.githubusercontent.com/sammmlow/ORBITM/master/gui/orbm_flow_sams.png)

On the other hand, STK computes a full orbit propagation, inclusive of all the high precision perturbing forces as per the diagram below. The default Jacchia-Roberts atmospheric density from STK's HPOP was used. ORBITM takes only the computed Delta-V values from STK, as well as the thrust time-schedule "deltaV.txt" file. A rough flow chart, using STK's Astrogator module, with a high precision orbit propagator (HPOP), works as follows:

![Orbit.M - Flow Chart for Orbit Maintenance with STK Astrogator](https://raw.githubusercontent.com/sammmlow/ORBITM/master/gui/orbm_flow_stk.png)

### Some Final Notes:

I am also looking for users experienced in NASA's GMAT program to automate orbit maintenance routines uing Python. If you like open source, and you enjoy orbital mechanics, do consider adding your contributions here for other GMAT users too!

This project is open source. If you would like to contribute to this project, add in new features, or enhance existing atmospheric models etc, please free free to fork my repo. Or, if you had felt that OrbitM was useful in your research, please do give the due credit and cite my paper or this project (the APA format is below)!

Low, S. Y. W., &; Chia, Y. X. (2018). “Assessment of Orbit Maintenance Strategies for Small Satellites”, 32nd Annual AIAA/USU Conference on Small Satellites, Logan, Utah, Utah State University, USA.

### Contact:

If you have any other queries feel free to reach out to me at:
sammmlow@gmail.com
