
#ifndef __CONSTANTS_HPP__
#define __CONSTANTS_HPP__

#include <cmath>

const double mass_sun   = 1.989e30; //kg
const double UA_m       =1.495978707e11; //m

const double G_SI       = 6.6743e-11; // m3 kg-1 s-2
const double G_UA_MSun  = G_SI / pow(1.495978707, 3) * 1.989e-3; // UA3 msun-1 s-2

const double mass_earth = 5.972e24; //kg
const double earth_obital_speed = 29.78e3; // m/s
const double earth_obital_speed_UA = 29.78e3/UA_m; // UA/s


const double moon_distance_to_earth = 384400e3; //m
const double moon_mass =  7.342e22; //kg
const double moon_obital_speed = 1.0e3; // m/s
const double moon_obital_speed_UA = moon_obital_speed/UA_m; // UA/s



const double saturn_distance = 9.5; //UA
const double saturn_mass= 5.683e26; //kg
const double saturn_orbital_speed =9.69e3; // m/s
const double saturn_orbital_speed_UA = saturn_orbital_speed/UA_m; // UA/s


const double distance_jupiter = 5.2; //UA
const double mass_jupiter = 1.898e27; //kg
const double jupiter_orbital_speed = 13.07e3; // m/s
const double jupiter_orbital_speed_UA = jupiter_orbital_speed/UA_m; // UA/s


const double seconds_per_year = 365.26*24*3600;
const double G_UA_MSun_yr = G_SI / pow(1.495978707, 3) * 1.989e-3 * pow(seconds_per_year,2); // UA3 msun-1 s-2
const double earth_obital_speed_UA_yr = earth_obital_speed/UA_m*seconds_per_year; // m/s

const double K_SI = 9e9; //N m2 / C2
const double e_ch = 1.602176634e-19; //C

const double K_SI_e = K_SI * e_ch * e_ch ; //N m2
const double K_SI_e_nm = K_SI_e *1e18 ; //N nm2

const double UAM_kg = 1.66e-27; //kg


#endif
