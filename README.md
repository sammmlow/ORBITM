![Orbit.M - Orbit Maintenance Analysis for LEO in Python](https://raw.githubusercontent.com/sammmlow/ORBITM/master/gui/orbm_logo_large.png)

### What is Orbit.M?

Orbit.M is an open-source easy-to-use orbit maintenance simulator and propulsion sizing tool for anyone and anywhere! It is a software with a graphical user interface where you can fill in your orbital parameters, satellite masses, satellite areas, drag parameters, and the software will compute your desired Delta-V necessary for the mission, sizing it against your choices of thrusters (based on the Isp). For beginners who are not familiar with the terms "Delta-V" and "Isp", the Wikipedia articles about them, below, do a fantastic job of explaining what they are.

[WikiPedia article on Delta-V.](https://en.wikipedia.org/wiki/Delta-v)

[WikiPedia article on Specific Impulse, or Isp.](https://en.wikipedia.org/wiki/Specific_impulse)

Once you input your orbit and satellite parameters in the GUI, and you have configured your desired thrusters that you want to use in your mission Delta-V sizing profile, the software will output a mission Delta-V profile against the Isp of your thruster choices, as well as an altitude and mean semi-major axis chart. There are examples showing these charts down below.

### How do I use Orbit.M specifically?

First things first, check that you have these Python libraries: **TKinter, NumPy, Matplotlib**, and among other standard libraries such as **os, datetime, comtypes etc**. If you do not have Python, I recommend using the [Anaconda installer](https://www.anaconda.com/), and running it on Spyder 4 simply because that's what I'm using right now. At the time of the latest commit and push, NumPy version is 1.18, and Matplotlib is 3.1.3. Older versions should work fine as ORBIT.M doesn't use any of the newer features (let me know if it does not).

Now, the software can be started by running the Python file **orbitm.py** in the main directory (which is equivalent to the directory you see on the master branch on ORBITM's github page). You should see a GUI that looks like the one below, pop up:

![Orbit.M - GUI running on TKinter](https://raw.githubusercontent.com/sammmlow/ORBITM/master/gui/orbm_screenshot.png)



### Contact:

If you have any queries feel free to reach out to me at:
sammmlow@gmail.com
