// This file is part of MEP.
//
//    Copyright 2014 Arto Sandroos.
//
//    MEP is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    MEP is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with MEP.  If not, see <http://www.gnu.org/licenses/>.

#ifndef CONSTANTS_H
#define CONSTANTS_H

/* This file contains some physical constants and
 * some stellar body parameters.
*/

// Magnetic constant, permeability of vacuum.
// Unit N/A^2 (Newton per Ampere squared).
const double MU_0 = 12.566371e-07;
// Electric constant, permittivity of vacuum.
// Unit F/m (Faraday per meter).
const double EPSILON_0 = 8.854188e-12;
// Newtonian constant of gravitation. 
// Unit m^3/kg/s^2 (cubic meter per kilogram per second squared).
const double GRAV_CONSTANT = 6.6742e-11;
// Planck constant.
// Unit Js (Joule times second).
const double PLANCK = 6.626069e-34;
// Planck constant in eVs.
// Unit eVs (electron Volt times second).
const double PLANCK_eVs = 4.135667e-15;
// Speed of light in vacuum.
// Unit m/s (meters per second).
const double SPEED_OF_LIGHT = 299792458.0;
// Elementary charge.
// Unit C (Coulomb).
const double ELEMENTARY_CHARGE = 1.602177e-19;

// Astronomical unit.
// Unit m (meters).
const double AU = 149598000000.0;
const double DAY_IN_SECONDS = 86400.0;
const double YEAR_IN_SECONDS = 3.1536e+7;

// Mass of electron.
// Unit kg (kilogram).
const double ELECTRON_MASS = 9.109383e-31;
// Electron mass energy equivalent.
// Unit MeV (10^6 electron volts).
const double ELECTRON_MASS_MeV = 0.510998;
// Mass of proton.
// Unit kg (kilogram).
const double PROTON_MASS = 1.672622e-27;
// Proton mass energy equivalent.
// Unit MeV (10^6 electron volts).
const double PROTON_MASS_MeV = 938.272029;
// Mass of alpha particle.
// Unit kg (kilogram).  
const double ALPHA_MASS = 6.644657e-27;
// Alpha particle mass energy equivalent.
// Unit MeV (10^6 electron volts).
const double ALPHA_MASS_MeV = 3727.37917;


// Mass of Earth.
// Unit kg (kilograms).
const double EARTH_MASS = 5.9736e+24;
// Equatorial radius of Earth.
// Unit m (meters).
const double EARTH_EQUAT_RADIUS = 6378136.0;
// Mass of the sun.
// Unit kg (kilograms).
const double SUN_MASS = 1.989e+30;
// Equatorial radius of the sun.
// Unit m (meters).
const double SUN_EQUAT_RADIUS = 695000000.0;
const double MERCURY_EQUAT_RADIUS = 2439700.0; /**< Equatorial radius of the mercury (meters). */

#endif
