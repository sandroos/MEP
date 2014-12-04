MEP
===

A code for modeling Magnetospheric Energetic Particles.

GENERAL
=======

The simulation is compiled by typing
> make mercury

After a successful compilation, simulation can be run by typing
> ./mercury <parameter file> <output file>

where <parameter file> is the name of the parameter file (see below)
      <output file>    is the name of the output file

The output file contains the injection, future, and past states of each particle. Here 
"future" means that the trajectory has been calculated forward in time, and "past" that 
the trajectory has been calculated backwards in time. The output file is in ascii.

Per column, the contents are: 

1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20
gi  vxi vyi vzi t1  x1  y1  z1  g1 vx1 vy1 vz1 t2  x2  y2  z2  g2  vx2 vy2 vz2

Here i = injection state
     1 = future state
     2 = past state
	 
    t           = exit time
    g           = gamma-factor
    x,y,z       = coordinates of the exit point
    vx,vy,vz    = components of the velocity

There are a few example parameter files:

./mercury params.txt o.txt

If o.txt is used as the name of the output file, there are gnuplot macros that 
can be used to visualize the results:

> gnuplot
gnuplot> load "gp_obs"
gnuplot> load "gp_fobs"

The 3D plot can be rotated using mouse. There is a very crude wireframe model of 
the surface of the Mercury to aid the eye. BepiColombo's nadir- and velocity vectors 
are drawn using long vectors, and the instrument directions are drawn using shorter, 
numbered, vectors. The future exit points are drawn using red symbols, and past exit 
points using green symbols. "gp_fobs" is the same as "gp_obs", except that it also 
draws field lines that have been calculated to lines.txt (see below).

There are at the moment two magnetic field configurations that can be used. The 
selection between these is made at compile time, i.e., you need to edit the FIELD-parameter 
in the Makefile (either dipole.o or mirrordipole.o). Alternatively, you can set the 
magnetic field at the command line 

> make mercury FIELD=dipole.o
> make mercury FIELD=mirrordipole.o

You can also calculate and print the trajectory of a single particle. In this case one writes to the terminal
> make clean
> make mercury CXXFLAGS=-DTRAJ

The simulation is ran exactly in the same manner as above:
> ./mercury params_chaotic.txt o.txt
> ./mercury params_chaotic2.txt o.txt

If the output file is named as o.txt, there are again gnuplot macro for fast visualization.
> gnuplot
gnuplot> load "gp_traj"
gnuplot> load "gp_ftraj"

The BepiColombo vectors and Hermean surfaces are drawn as explained above. The particle future and 
past trajectories are drawn using a curve. NOTE: both the future and past trajectories are drawn using 
the same color! "gp_ftraj" is the same as "gp_traj", except that it also draws fields lines from 
"lines.txt" (see below).

The output file has a different format in this case. It is as follows:
`1  2  3  4  5  6  7  8  9  10`
`t  x  y  z  vx vy vz Bx By Bz`

The two example parameter files named as "_chaotic.txt" should show chaotic orbits of ~300 keV protons. 
The file "params_mirror.txt" does the same for the mirror dipole field configuration, but 
using 1 MeV protons.
> ./mercury params_mirror.txt o.txt

Note that the code needs to be compiled using 
> make mercury "FIELD=mirrordipole.o" "FLAGS=-DTRAJ"

INPUT PARAMETERS FOR THE SIMULATION
===================================

RandomNumbers
-------------

seed                      A negative integer, seed value for the random number generator. 
                          If not given, seed is calculated from system clock.
						  
Dipole
------

B_pole_reference(T)       Magnitude of the magnetic field, in Teslas, at the poles of the dipole.

R_pole(km)                Distance from dipole origin to the pole where B_pole_reference(T) is 
                          specified, in km. If MERCURY is given, the Hermean equatorial radius is used.

x_ori,y_ori,z_ori         These defines the orientation of the dipole moment, positive pole is 
                          to the direction of this vector. This does not need to be a unit vector, 
						  it will be normalized to unity. 

offset_angle(deg)         The Hermean dipole moment does not reside at the center of the planet. 
                          This defines the angle between the dipole origin and x-axis, in the xy-plane.
						  
offset(km)                The distance between the center of Mercury and the dipole origin, in km. 
                          The dipole origin will be (r cos(a),r sin(a), 0), where r is the value 
						  given here, and a is the offset_angle(deg). 

minradius(hermean)        Minimum allowed radius for the particles. If r <= radius, the particle is 
                          removed from the simulation.

maxradius(hermean)        Maximum allowed radius for the particles. If r >= radius, the particle is 
                          removed from the simulation.

Injector
--------

Instrument                The number of the instrument (1-5) that should be used in the simulation. 
                          Only the particles that have entered this instrument will be calculated.

conehalfwidth(deg)        Half width of the observation cone of the chosen instrument, in degrees.

particles(int)            Number of simulated particles.

q(elementary)             Charge of the simulated particles, in elementary charges (proton charge = 1.0 etc).
1
m(protons)                Mass of the simulated particles, in proton masses. 

energy_min(keV)           Minimum energy of the simulated particles, in keV.

energy_max(keV)           Maximum energy of the simulated particles, in keV. Must be > energy_min(keV).

pwrlawindex               Particles are injected using a power law distribution. This is the power law index.

x_spacecraft(km)          Coordinates of the injection position, i.e., coordinates of the spacecraft, in km.
y_spacecraft(km)
z_spacecraft(km)

x_velocity                A vector pointing to the direction of the spacecraft velocity. Currently this vector must 
y_velocity                lie in the xz-plane. The nadir direction is defined to be orthogonal to this vector in xz-plane
z_velocity                also, and point towards the Mercury. 

                          Velocity and nadir vectors are used to construct an orthogonal basis, where the z'-axis is 
						  the nadir direction, x'-axis is to the direction of the velocity, and y'-axis completes the set. 
						  
                          Instrument #1 will point to the opposite direction from the nadir vector.
						  Instruments #2-#5 will lie in the x'y' -plane.

MirrorDipole
------------

B0_pole_reference(T)              A reference value for the planetary dipole field at the poles of 
                                  the planetary dipole field.

R0_pole(hermean)                  Distance from the planetary dipole moment to the poles in hermean radii.
                                  This value, together with B0_pole_reference(T), are used to calculate the 
								  dipole moment of the planetary dipole.
								  
tiltangle(deg)                    The tilt angle of the planetary dipole in the yz-plane, in degrees. The 
                                  mirror dipole will have the same tilt angle.

B1_pole_reference(T)              A reference value for the mirror dipole field at the poles of the mirror dipole field.

R1_pole(hermean)                  Distance from the mirror dipole moment to the poles in hermean radii. See 
                                  R0_pole(hermean).

mirrordipole_distance(hermean)    The distance between the center of Mercury and the mirror dipole moment along 
                                  the x-axis, in hermean radii. In other words, the origin of the mirror dipole is 
								  at coordinates (x,0,0), where x is the value of this parameter. 

minradius(hermean)                Minimum allowed radius for the particles. If r <= radius, the particle is
                                  removed from the simulation.

maxradius(hermean)                Maximum allowed radius for the particles. However, the origin of this "exit sphere"
                                  is not at the center of the planet, but at the coordinates specified below.
								  
x_maxradius(hermean)              x-coordinate of the center of the outer "exit sphere".

y_maxradius(hermean)              y-coordinate of the center of the outer "exit sphere".

z_maxradius(hermean)              z-coordinate of the center of the outer "exit sphere".


Simulation
----------

maxtime_per_dt            Maximum simulation time for each particle, in time steps. The time step 
                          is calculated in the simulation. If this simulation time is reached, the 
						  particle will be discarded.

saveinterval_dt           Defines the interval the trajectory of the particle is written into the output file.
                          This parameter is used only when the simulation has been compiled with the 
						  -DTRAJ option (see above).

CALCULATION OF FIELD LINES FOR THE GNUPLOT VISUALIZATION
========================================================

There is a very crude and barbaric program that will calculate some field lines to an ascii 
file, which can then be included in gnuplot visualization. First, it needs to be compiled 
with

> make flines

Make sure that in the Makefile the value of the FIELD parameter is the same that has been used 
when compiling the particle mover simulation, or otherwise the field lines are calculated from 
a different magnetic field :)

NOTE: If you change the parameters of the magnetic field in the mercury-simulation, you also 
      need to rerun flines!

If you do not want to edit the Makefile, you can also set it as 

> make flines FIELD=dipole.o
> make flines FIELD=mirrordipole.o

The program itself is run by typing

> ./flines <tilt angle> <N_theta> <N_phi> <parameter file> > lines.txt

where
   <tilt angle>          is a tilt angle that is used when calculating the starting points 
                         for each field line. Suggestion is to use the same value that has 
						 been used for the magnetic field in the parameter file.
						 
   <N_theta>             Number of theta values used. 
   <N_phi>               Number of phi values used. The total number of calculated field lines 
                         is N_theta * N_phi. 
						 
   <parameter file>      Name of the parameter file for the mercury-simulation, this is needed 
                         to set the same parameters for the magnetic field that have been used 
						 in the mercury simulation.

Example:

> ./flines 14.0 3 4 params.txt > lines.txt
