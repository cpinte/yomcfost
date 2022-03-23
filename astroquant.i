
/* Author :
 * $Log: astroquant.i,v $
 * Revision 1.5  2004/07/13 14:46:20  jblebou
 * test a new version
 *
 * Revision 1.4  2004/07/13 14:40:15  guieu
 *  - carrected a bug in the first WRITE function
 *
 * Revision 1.3  2004/07/09 14:26:06  jblebou
 * global revision of all the files :
 *  - replace the old file name (toto_v1.i) by new one (tot.i)
 *  - add some author name
 *  - clean the CVS history
 *
 * Revision 1.2  2004/07/09 10:27:45  fmillour
 *
 */

local astroquant
/* DOCUMENT astroquant.i -- main physical and astrophysical units
   units are accessed via the preferred unit system
   examples:
       3*CGS.N:      3 newtons (in SI)
       SI.mug:       1 microgram (in MKSA)
       SI.band.R:    wavelength of the R-band
       SI.band.Kwd:  width of the K-band
       SI.band.Uf0:  flux at U = 0 mag
       10*SI.mJy:    10 millijanskys
       SI.Rsun:      solar radius
       SI.grav:      gravitation constant
       SI.c:         speed of light

   and for conversion :
       RadtoDeg = 180./ pi;
       DegtoRad = pi / 180.;
       MastoRad = pi/648000000.;
       RadtoMas = 648000000./pi;
       AngtoMicron = 1.e-4;
       MicrontoAng = 1.e+4;
       AngtoMeter = 1.e-10;
       MetertoAng = 1.e+10;
       eVtoJoule = 1.602177e-19;
       eVtocm_1 = 8065.48;

   *Functions: - SI
               - CGS
               - SI_extinction
               - SI_Blambda
               ...

   *Require: nothing

   *History:
      - *** Malbet: creation.
      - 03/11/28 LeBouquin: put in this dir.

   SEE ALSO:
 */
/* ----------------------------------------------------------------------- */
if(_PRINT) write, "#include \"astroquant.i\"";
/* ----------------------------------------------------------------------- */

struct spectral_band
/* DOCUMENT spectral band

   Structure containing information on spectral bands
   - mean wavelength:      U, B, V, etc.
   - width:                Uwd, Bwd, Vwd, etc.
   - flux for magnitude 0: Uf0, Bf0, Vf0
 */
{
   double U, B, V, R, I, J, H, K, L, M, N, Q, K12, K25, K60, K100;
   double Uwd, Bwd, Vwd, Rwd, Iwd, Jwd, Hwd, Kwd, Lwd, Mwd, Nwd, Qwd, K12wd, K25wd, K60wd, K100wd;
   double Uf0, Bf0, Vf0, Rf0, If0, Jf0, Hf0, Kf0, Lf0, Mf0, Nf0, Qf0, K12f0, K25f0, K60f0, K100f0;
};

/* ----------------------------------------------------------------------- */
struct unit_system
/* DOCUMENT unit_system

   Structure containing values of main units (mostly CGS and SI units with
   their most used multiples) and some (astro-)physical quantities. This
   structure has two instances called SI and CGS, which express these units
   and quantities in SI (MKSA) and CGS.

   To access a unit, for instance kg, one uses SI.kg or CGS.kg according
   to your favorite system.

   Most units are called by their standard name.  Some have non-trivial
   names:
      hr (hour), ps suffix (per second), mu prefix (micro), angstroem,
      Ohm (ohm), tr (2pi rad), deg (degree) as (second of arc),
      mmHg (pressure of 1 mm of mercury), Jy (1e-26 W/m²), amu (atomic
      mass unit), Rsun, Msun, Lsun (radius, mass and luminosity of the
      sun), Mearth, Rearth (mass and radius of the earth), hbar (reduced
      Planck constant), kB (Boltzman constant), Grav (gravity constant)

   Note: candela, lux and lumen are not supported.

   Examples:

      Specific flux per unit of bandwidth from a body at temperature 170 mK
      near 220 µm.

         sflux_si  = ;SI_Blambda(220* SI.mum, 170 * SI.mK)
         spflux_cgs = SI_Blambda(220*CGS.mum, 170 * SI.mK)

      Express 100 light years in solar radii
         lyr_rsun   = (100 * SI.c  * SI.yr ) / SI.Rsun,  or

      Mass of 0.13 mmol of carbon 14.

         mass_mg    = (0.13 * SI.mmol)  * SI.Na * ( 14 * SI.u) / SI.mg
*/
{
   double tr, rad, deg, as, mas, muas;
   double mus, ms, s, mn, hr, day, yr, Myr, Gyr;
   double Hz, kHz, MHz, GHz, Ci;
   double fm, pm, nm, mum, mm, cm, m, km, Rearth, Rsun, AU, pc, kpc, Mpc, Gpc;
   double mL, cL, L, hL;
   double cmps, mps, kps, kph, c;
   double mug, mg, g, kg, t, me, mp, amu, Msun, Mearth;
   double N, daN, kgf, dyn;
   double eV, keV, MeV, GeV, erg, mJ, J, kJ, MJ, GJ;
   double muW, mW, W, kW, MW, GW, TW, ergps, Lsun;
   double muT, mT, T, muG, mG, G;
   double mumol, mmol, mol;
   double muK, mK, K;
   double mu0;
   double mJy, Jy;
   double C, epsilon0, e;
   double muA, mA, A;
   double muV, mV, V, kV;
   double pF, nF, muF, mF, F;
   double muH, mH, H;
   double Ohm, kOhm, MOhm;
   double muS, mS, S;
   double Pa, hPa, mbar, bar, atm, mmHg;
   double Grav, kB, Na, h, hbar, sigma, alpha;
   spectral_band band;
   double _c1, _c2;
};

/* ----------------------------------------------------------------------- */
SI = unit_system(
   tr        = 3.141592653590e+00,
   rad       = 1e+00,
   deg       = 1.745329251994e-02,
   as        = 4.848136811095e-06,
   mas       = 4.848136811095e-09,
   muas      = 4.848136811095e-12,
   ms        = 1e-03,
   s         = 1e+00,
   mn        = 6e+01,
   hr        = 3.6e+03,
   yr        = 3.15569259747e+07,
   Myr       = 3.15569259747e+13,
   Gyr       = 3.15569259747e+16,
   day       = 8.64e+04,
   Hz        = 1e+00,
   kHz       = 1e+03,
   MHz       = 1e+06,
   GHz       = 1e+09,
   Ci        = 3.7e+10,
   fm        = 1e-15,
   pm        = 1e-12,
   nm        = 1e-09,
   mum       = 1e-06,
   mm        = 1e-03,
   cm        = 1e-02,
   m         = 1e+00,
   km        = 1e+03,
   Rsun      = 6.9598e+08,
   Rearth    = 6.378164e+06,
   AU        = 1.495985e+11,
   pc        = 3.085e+16,
   kpc       = 3.085e+19,
   Mpc       = 3.085e+22,
   Gpc       = 3.085e+25,
   mL        = 1e-06,
   cL        = 1e-04,
   L         = 1e-03,
   hL        = 1e-01,
   cmps      = 1e-02,
   mps       = 1e+00,
   kps       = 1e+03,
   kph       = 2.777777777778e-01,
   c         = 2.99792453e+08,
   mug       = 1e-09,
   mg        = 1e-06,
   g         = 1e-03,
   kg        = 1e+00,
   t         = 1e+03,
   me        = 9.10938188e-33,
   mp        = 1.67262158e-27,
   amu       = 1.66053873e-27,
   Msun      = 1.9892e+30,
   Mearth    = 5.976e+24,
   N         = 1e+00,
   daN       = 1e+01,
   kgf       = 9.80665,
   dyn       = 1e-05,
   eV        = 1.602176462e-19,
   keV       = 1.602176462e-16,
   MeV       = 1.602176462e-13,
   GeV       = 1.602176462e-10,
   erg       = 1e-07,
   mJ        = 1e-03,
   J         = 1e+00,
   kJ        = 1e+03,
   MJ        = 1e+06,
   GJ        = 1e+09,
   muW       = 1e-06,
   mW        = 1e-03,
   W         = 1e+00,
   mJy       = 1e-29,
   Jy        = 1e-26,
   kW        = 1e+03,
   MW        = 1e+06,
   GW        = 1e+09,
   TW        = 1e+12,
   ergps     = 1e-07,
   Lsun      = 3.826e+26,
   muT       = 1e-06,
   mT        = 1e-03,
   T         = 1e+00,
   muG       = 1e-10,
   mG        = 1e-07,
   G         = 1e-04,
   mu0       = 1.2566370615e-6,
   C         = 1e+00,
   epsilon0  = 8.854187817e-12,
   e         = 1.602176462e-19,
   muA       = 1e-06,
   mA        = 1e-03,
   A         = 1e+00,
   muV       = 1e-06,
   mV        = 1e-03,
   V         = 1e+00,
   kV        = 1e+03,
   muH       = 1e-06,
   mH        = 1e-03,
   H         = 1e+00,
   pF        = 1e-12,
   nF        = 1e-09,
   muF       = 1e-06,
   mF        = 1e-03,
   F         = 1e+00,
   Ohm       = 1e+00,
   kOhm      = 1e+03,
   MOhm      = 1e+06,
   muS       = 1e-06,
   mS        = 1e-03,
   S         = 1e+00,
   Pa        = 1e+00,
   hPa       = 1e+02,
   mbar      = 1e+02,
   bar       = 1e+05,
   atm       = 1.01325e+05,
   mmHg      = 1.3332e+02,
   Na        = 6.02214199e23,
   Grav      = 6.673e-11,
   kB        = 1.3806503e-23,
   sigma     = 5.670400e-8,
   muK       = 1e-06,
   mK        = 1e-03,
   K         = 1e+00,
   mumol     = 1e-06,
   mmol      = 1e-03,
   mol       = 1e+00,
   h         = 6.62606876e-34,
   hbar      = 1.05457160e-34,
   alpha     = 1/137.04,
   band      = spectral_band(
      U    = 3.67e-7, Uwd = 6.60e-8, Uf0    = 3.9811e-02,
      B    = 4.36e-7, Bwd = 9.40e-8, Bf0    = 6.3096e-02,
      V    = 5.45e-7, Vwd = 8.80e-8, Vf0    = 3.6308e-02,
      R    = 6.38e-7, Rwd = 1.38e-7, Rf0    = 2.2387e-02,
      I    = 7.97e-7, Iwd = 1.49e-7, If0    = 1.1482e-02,
      J    = 1.22e-6, Jwd = 2.13e-7, Jf0    = 3.1623e-03,
      H    = 1.63e-6, Hwd = 3.07e-7, Hf0    = 1.1482e-03,
      K    = 2.19e-6, Kwd = 3.90e-7, Kf0    = 3.9811e-04,
      L    = 3.45e-6, Lwd = 4.72e-7, Lf0    = 7.0795e-05,
      M    = 4.75e-6, Mwd = 4.60e-7, Mf0    = 2.0417e-05,
      N    = 1.02e-5, Nwd = 4.00e-6, Nf0    = 1.2303e-06,
      Q    = 2.10e-5, Qwd = 5.00e-6, Qf0    = 6.7608e-08,
      K12  = 1.20e-5,                K12f0  = 5.89175e-07,
      K25  = 2.50e-5,                K25f0  = 3.22817e-08,
      K60  = 5.00e-5,                K60f0  = 9.90981e-10,
      K100 = 1.00e-4,                K100f0 = 1.28911e-10
   ),
   _c1       = 1.19107e-16,
   _c2       = 0.0143883

);

func SI_Blambda(lambda, T)
/* DOCUMENT SI_Blambda(lambda, T)

   Black body function in SI
   - lambda: wavelength (m)
   - T:      temperature (K)
   - return: W/m²/m
 */
{
   local B, nan;

   // nan = exp(1e+50);
   B   = SI._c1 / lambda^5 / (exp( min(SI._c2 / (T*lambda), 700)) - 1);
   //if (numberof((idx = where(B == nan))) != 0)
   //   B(idx) = 0;
   return B;
}

func SI_dBlambda_dT(lambda, T)
/* DOCUMENT SI_dBlambda_dT(lambda, T)

   Derivative of black body function in respect to T
   - lambda: wavelength (m)
   - T:      temperature (K)
   - return: W/m²/m
 */
{
   return SI._c1 / lambda^5 / (exp(SI._c2 / (T*lambda)) - 1)^2 \
          * exp(SI._c2/(T*lambda)) * SI._c2/(T^2*lambda);
}

func SI_B(T)
/* DOCUMENT SI_B(T)

   Black body function in SI
   - T:      temperature (K)
   - return: W/m²/m
 */
{
   return SI.sigma * T^4;
}

func SI_extinction(wavelength, Av)
/* DOCUMENT SI_extinction(wavelength, Av)

   Return the extinction at `wavelength' knowing the one in V `Av. SI system.
 */
{
  bigr = 3.1;
  wavlaw = [1.e-5,1.05e-5,1.11e-5,1.18e-5,1.25e-5,
	    1.39e-5,1.49e-5,1.6e-5,1.7e-5,1.8e-5,1.9e-5,2.e-5,2.1e-5,
	   2.19e-5,2.3e-5,2.4e-5,2.5e-5,2.74e-5,3.44e-5,4.e-5,4.4e-5,
	   5.5e-5,7.e-5,9.e-5,1.25e-4,2.2e-4,3.4e-4,1.e0];
  extlaw = [11.3e0,9.8e0,8.45e0,7.45e0,6.55e0,5.39e0,
	   5.05e0,5.02e0,4.77e0,4.65e0,4.9e0,5.52e0,6.23e0,6.57e0,5.77e0,
	   4.9e0,4.19e0,3.1e0,1.8e0,1.3e0,1.e0,0.e0,-0.78e0,-1.6e0,
	   -2.23e0,-2.72e0,-2.94e0,-3.1e0];

  return (interp(extlaw, wavlaw, wavelength/SI.cm)+bigr)*Av/bigr;

}

func SI_Jy2flux (jy, lambda)
/* DOCUMENT I_jy2flux (jy, lambda)

   Return the flux per unit of band width. (W/m³)

 */
{
   return jy / lambda^2 * SI.c * SI.Jy;
}

func SI_flux2Jy (flux, lambda)
/* DOCUMENT  SI_flux2Jy (flux, lambda)

   Return the flux in Jy.

 */
{
   return flux * lambda^2 / SI.c / SI.Jy;
}

func SI_mag2flux (mag, band) {

   local i, flux;

   flux = array(double, dimsof(mag));

   for (i = 1; i <= numberof(band); ++i) {
      if (band(i) == "U")    { flux(i) =  SI.band.Uf0    * 10 ^ (-0.4*mag(i)); continue; }
      if (band(i) == "B")    { flux(i) =  SI.band.Bf0    * 10 ^ (-0.4*mag(i)); continue; }
      if (band(i) == "V")    { flux(i) =  SI.band.Vf0    * 10 ^ (-0.4*mag(i)); continue; }
      if (band(i) == "R")    { flux(i) =  SI.band.Rf0    * 10 ^ (-0.4*mag(i)); continue; }
      if (band(i) == "I")    { flux(i) =  SI.band.If0    * 10 ^ (-0.4*mag(i)); continue; }
      if (band(i) == "J")    { flux(i) =  SI.band.Jf0    * 10 ^ (-0.4*mag(i)); continue; }
      if (band(i) == "H")    { flux(i) =  SI.band.Hf0    * 10 ^ (-0.4*mag(i)); continue; }
      if (band(i) == "K")    { flux(i) =  SI.band.Kf0    * 10 ^ (-0.4*mag(i)); continue; }
      if (band(i) == "L")    { flux(i) =  SI.band.Lf0    * 10 ^ (-0.4*mag(i)); continue; }
      if (band(i) == "M")    { flux(i) =  SI.band.Mf0    * 10 ^ (-0.4*mag(i)); continue; }
      if (band(i) == "N")    { flux(i) =  SI.band.Nf0    * 10 ^ (-0.4*mag(i)); continue; }
      if (band(i) == "Q")    { flux(i) =  SI.band.Qf0    * 10 ^ (-0.4*mag(i)); continue; }
      if (band(i) == "K12")  { flux(i) =  SI.band.K12f0  * 10 ^ (-0.4*mag(i)); continue; }
      if (band(i) == "K25")  { flux(i) =  SI.band.K25f0  * 10 ^ (-0.4*mag(i)); continue; }
      if (band(i) == "K60")  { flux(i) =  SI.band.K60f0  * 10 ^ (-0.4*mag(i)); continue; }
      if (band(i) == "K100") { flux(i) =  SI.band.K100f0 * 10 ^ (-0.4*mag(i)); continue; }
   }

   return flux;

}

func SI_flux2mag (flux, band) {

   local i, jy;

   jy = array(double, dimsof(flux));

   for (i = 1; i <= numberof(band); ++i) {
      if (band(i) == "U")    { jy(i) =  -2.5*log10(flux(i)/SI.band.Uf0   ); continue; }
      if (band(i) == "B")    { jy(i) =  -2.5*log10(flux(i)/SI.band.Bf0   ); continue; }
      if (band(i) == "V")    { jy(i) =  -2.5*log10(flux(i)/SI.band.Vf0   ); continue; }
      if (band(i) == "R")    { jy(i) =  -2.5*log10(flux(i)/SI.band.Rf0   ); continue; }
      if (band(i) == "I")    { jy(i) =  -2.5*log10(flux(i)/SI.band.If0   ); continue; }
      if (band(i) == "J")    { jy(i) =  -2.5*log10(flux(i)/SI.band.Jf0   ); continue; }
      if (band(i) == "H")    { jy(i) =  -2.5*log10(flux(i)/SI.band.Hf0   ); continue; }
      if (band(i) == "K")    { jy(i) =  -2.5*log10(flux(i)/SI.band.Kf0   ); continue; }
      if (band(i) == "L")    { jy(i) =  -2.5*log10(flux(i)/SI.band.Lf0   ); continue; }
      if (band(i) == "M")    { jy(i) =  -2.5*log10(flux(i)/SI.band.Mf0   ); continue; }
      if (band(i) == "N")    { jy(i) =  -2.5*log10(flux(i)/SI.band.Nf0   ); continue; }
      if (band(i) == "Q")    { jy(i) =  -2.5*log10(flux(i)/SI.band.Qf0   ); continue; }
      if (band(i) == "K12")  { jy(i) =  -2.5*log10(flux(i)/SI.band.K12f0 ); continue; }
      if (band(i) == "K25")  { jy(i) =  -2.5*log10(flux(i)/SI.band.K25f0 ); continue; }
      if (band(i) == "K60")  { jy(i) =  -2.5*log10(flux(i)/SI.band.K60f0 ); continue; }
      if (band(i) == "K100") { jy(i) =  -2.5*log10(flux(i)/SI.band.K100f0); continue; }
   }

   return jy;

}



/* ----------------------------------------------------------------------- */
CGS = unit_system(
   tr        = 3.141592653590e+00,
   rad       = 1e+00,
   deg       = 1.745329251994e-02,
   as        = 4.848136811095e-06,
   mas       = 4.848136811095e-09,
   muas      = 4.848136811095e-12,
   ms        = 1e-03,
   s         = 1e+00,
   mn        = 6e+01,
   hr        = 3.6e+03,
   day       = 8.64e+04,
   yr        = 3.15569259747e+07,
   Myr       = 3.15569259747e+13,
   Gyr       = 3.15569259747e+16,
   Hz        = 1e+00,
   kHz       = 1e+03,
   MHz       = 1e+06,
   GHz       = 1e+09,
   Ci        = 3.7e+10,
   fm        = 1e-13,
   pm        = 1e-10,
   nm        = 1e-07,
   mum       = 1e-04,
   mm        = 1e-01,
   cm        = 1e+00,
   m         = 1e+02,
   km        = 1e+05,
   Rearth    = 6.378164e+08,
   Rsun      = 6.9598e+10,
   AU        = 1.495985e+13,
   pc        = 3.085e+18,
   kpc       = 3.085e+21,
   Mpc       = 3.085e+24,
   Gpc       = 3.085e+27,
   mL        = 1e+00,
   cL        = 1e+01,
   L         = 1e+03,
   hL        = 1e+05,
   cmps      = 1e+00,
   mps       = 1e+02,
   kps       = 1e+05,
   kph       = 2.777777777778e+01,
   c         = 2.99792453e+10,
   mug       = 1e-06,
   mg        = 1e-03,
   g         = 1e+00,
   kg        = 1e+03,
   t         = 1e+06,
   me        = 9.10938188e-30,
   mp        = 1.67262158e-24,
   amu       = 1.66053873e-24,
   Mearth    = 5.976e+27,
   Msun      = 1.9892e+33,
   N         = 1e+05,
   daN       = 1e+06,
   kgf       = 9.80665e+05,
   dyn       = 1e+00,
   eV        = 1.602176462e-12,
   keV       = 1.602176462e-09,
   MeV       = 1.602176462e-06,
   GeV       = 1.602176462e-03,
   erg       = 1e+00,
   mJ        = 1e+04,
   J         = 1e+07,
   kJ        = 1e+10,
   MJ        = 1e+13,
   GJ        = 1e+16,
   muW       = 1e+01,
   mW        = 1e+04,
   W         = 1e+07,
   kW        = 1e+10,
   MW        = 1e+13,
   GW        = 1e+16,
   TW        = 1e+19,
   ergps     = 1e+00,
   Lsun      = 3.826e+33,
   mJy       = 1e-26,
   Jy        = 1e-23,
   muT       = 1e-02,
   mT        = 1e+01,
   T         = 1e+04,
   muG       = 1e-06,
   mG        = 1e-03,
   G         = 1e+00,
   mu0       = 1.,
   C         = 1e+00,
   epsilon0  = 8.854187817e-12,
   e         = 1.602176462e-19,
   muA       = 1e-06,
   mA        = 1e-03,
   A         = 1e+00,
   muV       = 1e-06,
   mV        = 1e-03,
   V         = 1e+00,
   kV        = 1e+03,
   muH       = 1e-06,
   mH        = 1e-03,
   H         = 1e+00,
   pF        = 1e-12,
   nF        = 1e-09,
   muF       = 1e-06,
   mF        = 1e-03,
   F         = 1e+00,
   Ohm       = 1e+00,
   kOhm      = 1e+03,
   MOhm      = 1e+06,
   muS       = 1e-06,
   mS        = 1e-03,
   S         = 1e+00,
   Pa        = 1e+01,
   hPa       = 1e+03,
   mbar      = 1e+03,
   bar       = 1e+06,
   atm       = 1.01325e+06,
   mmHg      = 1.3332e+03,
   Na        = 6.02214199e23,
   Grav      = 6.673e-08,
   kB        = 1.3806503e-16,
   sigma     = 5.670400e-5,
   muK       = 1e-06,
   mK        = 1e-03,
   K         = 1e+00,
   h         = 6.62606876e-27,
   hbar      = 1.05457160e-27,
   alpha     = 1/137.04,
   mumol     = 1e-06,
   mmol      = 1e-03,
   mol       = 1e+00,
   band      = spectral_band(
      U    = 3.67e-5, Uwd = 6.60e-8, Uf0    = 3.9811e-01,
      B    = 4.36e-5, Bwd = 9.40e-8, Bf0    = 6.3096e-01,
      V    = 5.45e-5, Vwd = 8.80e-8, Vf0    = 3.6308e-01,
      R    = 6.38e-5, Rwd = 1.38e-7, Rf0    = 2.2387e-01,
      I    = 7.97e-5, Iwd = 1.49e-7, If0    = 1.1482e-01,
      J    = 1.22e-4, Jwd = 2.13e-7, Jf0    = 3.1623e-02,
      H    = 1.63e-4, Hwd = 3.07e-7, Hf0    = 1.1482e-02,
      K    = 2.19e-4, Kwd = 3.90e-7, Kf0    = 3.9811e-03,
      L    = 3.45e-4, Lwd = 4.72e-7, Lf0    = 7.0795e-04,
      M    = 4.75e-4, Mwd = 4.60e-7, Mf0    = 2.0417e-04,
      N    = 1.02e-3, Nwd = 4.00e-6, Nf0    = 1.2303e-05,
      Q    = 2.10e-3, Qwd = 5.00e-6, Qf0    = 6.7608e-06,
      K12  = 1.20e-3,                K12f0  = 5.8918e-06,
      K25  = 2.50e-3,                K25f0  = 3.2282e-07,
      K60  = 5.00e-3,                K60f0  = 9.9098e-09,
      K100 = 1.00e-3,                K100f0 = 1.2891e-09
   ),
   _c1       = 1.19107e-5,
   _c2       = 1.43883
);

func CGS_Blambda(lambda, T)
/* DOCUMENT CGS_Blambda(lambda, T)

   Black body function in CGS
   - lambda: wavelength (cm)
   - T:      temperature (K)
   - return: W/cm²/cm
 */
{
   return CGS._c1 / lambda^5 / (exp(CGS._c2 / (T*lambda)) - 1);
}

func CGS_B(T)
/* DOCUMENT CGS_B(T)

   Integrated black body function in CGS
   - T:      temperature (K)
   - return: erg/s/m²/m
 */
{
   return CGS.sigma * T^4;
}

func CGS_extinction(wavelength, Av)
/* DOCUMENT SI_extinction(wavelength, Av)

   Return the extinction at `wavelength' knowing the one in V `Av. CGS system.
 */
{
  bigr = 3.1;
  wavlaw = [1.e-5,1.05e-5,1.11e-5,1.18e-5,1.25e-5,
	    1.39e-5,1.49e-5,1.6e-5,1.7e-5,1.8e-5,1.9e-5,2.e-5,2.1e-5,
	   2.19e-5,2.3e-5,2.4e-5,2.5e-5,2.74e-5,3.44e-5,4.e-5,4.4e-5,
	   5.5e-5,7.e-5,9.e-5,1.25e-4,2.2e-4,3.4e-4,1.e0];
  extlaw = [11.3e0,9.8e0,8.45e0,7.45e0,6.55e0,5.39e0,
	   5.05e0,5.02e0,4.77e0,4.65e0,4.9e0,5.52e0,6.23e0,6.57e0,5.77e0,
	   4.9e0,4.19e0,3.1e0,1.8e0,1.3e0,1.e0,0.e0,-0.78e0,-1.6e0,
	   -2.23e0,-2.72e0,-2.94e0,-3.1e0];

  return (interp(extlaw, wavlaw, wav/CGS.cm)+bigr)*av/bigr;

}
/* ----------------------------------------------------------------------- */
extern RadtoDeg;
extern DegtoRad;
RadtoDeg = 180./ pi;
DegtoRad = pi / 180.;

extern RadtoMas;
extern MastoRad;
MastoRad = pi/648000000.;
RadtoMas = 648000000./pi;

extern AngtoMicron;
extern MicrontoAng;
AngtoMicron = 1.e-4;
MicrontoAng = 1.e+4;
extern AngtoMeter;
extern MetertoAng;
AngtoMeter = 1.e-10;
MetertoAng = 1.e+10;

extern eVtoJoule;
eVtoJoule = 1.602177e-19;

extern eVtocm_1;
eVtocm_1 = 8065.48;
