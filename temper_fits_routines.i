# include "ieee.i"
# include "lmfit.i"


NaN = array(float,1) ; ieee_set, NaN,2 ; NaN = NaN(1) ;
_wcs_bin_dir = "~/Software/wcstools-3.8.7/bin/" ;


//*************************************************************************

func define_opacity(wl,P) {
  // Define opacity law in m^2.kg^-1 of gas
  extern _kappa, _kappa0, _l0, _fact ;

  // Draine at 0.03 produit kappa = 0.2 cm^2/g de poussiere = 2e-4 em m^2.kg^-1 of gas
  _kappa0 = 0.1 * 1e-3 ; _l0 = 1.2e-3 ;

  // 1m^2/kg = 10 cm^2/g

  //for consistency with temper_fits.i and Philippe's code.
  _kappa0 = 0.1 * 1e-1;  // 1 cm^2/g = 0.1 m^2/kg
  write, "Setting nu0 to 1000Ghz, kappa0 to ",_kappa0 * 10," (cm^2.g^-1)" ;
  nu0 = 1000e9 ; // Hz ;
  _l0 = SI.c / nu0 ;

  P.kappa0 = _kappa0 ; P.nu0 = nu0 ;

  _kappa = _kappa0 * (wl/_l0)^-2 ;
  _fact = _kappa * SI.Msun / ((P.distance*SI.pc)^2) ;
}


//*************************************************************************

func zero_to_NaN(img) {

  // On mets les NaN a 0 d'abord
  ou = where(ieee_test(img)) ;
  if (numberof(ou)) {
    img(ou) = 0. ;
  }

   // On fait le test sur 0
  ou = where(abs(img)<1e-300) ;
  NaN = array(float,1) ; ieee_set, NaN,2 ; NaN = NaN(1) ;

  if (numberof(ou)) {
    img(ou) = NaN ;
  }
}

//*************************************************************************

func NaN_to_zero(img) {
  // Warning : not a function

  ou = where(ieee_test(img)) ;

  if (numberof(ou)) {
    img(ou) = 0. ;
  }
}

//*************************************************************************

func get_beam(conf,scan_speed) {
/* DOCUMENT  beam(conf,scan_speed)
   Returns instrument beam size in arcsec
   SEE ALSO:
 */

  beam_size = [] ;
  if (conf.instrument=="PACS" || conf.instrument=="SPIRE") {
    if (abs(scan_speed - 10) < 1e-3) {
      if (abs(conf.lamb - 70) < 1e-3) beam_size = 5.4 ;
      if (abs(conf.lamb - 100) < 1e-3) beam_size = 6.7 ;
      if (abs(conf.lamb - 160) < 1e-3) beam_size = 11.2 ;
    }

    if (abs(scan_speed - 20) < 1e-3) {
      if (abs(conf.lamb - 70) < 1e-3) beam_size = 5.6 ;
      if (abs(conf.lamb - 100) < 1e-3) beam_size = 6.8 ;
      if (abs(conf.lamb - 160) < 1e-3) beam_size = 11.4 ;
    }

    if (abs(scan_speed - 60) < 1e-3) {
      if (abs(conf.lamb - 70) < 1e-3) beam_size = 8.4 ;
      if (abs(conf.lamb - 100) < 1e-3) beam_size = 9.4 ;
      if (abs(conf.lamb - 160) < 1e-3) beam_size = 13.5 ;
    }

    if (abs(conf.lamb - 250) < 1e-3) beam_size = 18.1 ;
    if (abs(conf.lamb - 350) < 1e-3) beam_size = 25.2 ;
    if (abs(conf.lamb - 500) < 1e-3) beam_size = 36.9 ;

    if (is_void(beam_size)) {
      "ERROR: cannot obtain beam_size. Exiting." ;
      error ;
    }
  } else if (conf.instrument=="SABOCA" ||  conf.instrument=="LABOCA") {

    // Reads the fits file and pixel size
    image = cfitsio_open(conf.file, fh) ;

    // Beam size (circular for APEX)
    beam_size = cfitsio_get(fh,"BMAJ") * 3600.; // en arcsec
    if (is_void(beam_size)) {
      "ERROR: cannot obtain beam_size. Exiting." ;
      error ;
    }
  }

  return beam_size ;
}

//*************************************************************************

func read_Herschel(filename,instrument,wavelength,scan_speed,&pix_size,&beam_size,&conversion_factor,&med) {
/* DOCUMENT  read_Herschel(filename,instrument,wavelength,scan_speed,&pix_size,&beam_size,&conversion_factor,&med)
   Reads a fits file filename containing Herschel data
   Convert the flux to MJy/sr assuming PACS data is in Jy/pixel and SPIRE data in Jy/beam if the info is not

   instrument = "PACS" or "SPIRE"

   Returns:
     -- Herschel beam size in arcsec
     -- pixel size in arcsec
   SEE ALSO:
 */

  // Obtain beam_size as a function of wavelength and scan_speed
  beam_size = [] ;
  if (abs(scan_speed - 10) < 1e-6) {
    if (abs(wavelength - 70e-6) < 1e-6) beam_size = 5.4 ;
    if (abs(wavelength - 100e-6) < 1e-6) beam_size = 6.7 ;
    if (abs(wavelength - 160e-6) < 1e-6) beam_size = 11.2 ;
  }

  if (abs(scan_speed - 20) < 1e-6) {
    if (abs(wavelength - 70e-6) < 1e-6) beam_size = 5.6 ;
    if (abs(wavelength - 100e-6) < 1e-6) beam_size = 6.8 ;
    if (abs(wavelength - 160e-6) < 1e-6) beam_size = 11.4 ;
  }

  if (abs(scan_speed - 60) < 1e-6) {
    if (abs(wavelength - 70e-6) < 1e-6) beam_size = 8.4 ;
    if (abs(wavelength - 100e-6) < 1e-6) beam_size = 9.4 ;
    if (abs(wavelength - 160e-6) < 1e-6) beam_size = 13.5 ;
  }

  if (abs(wavelength - 250e-6) < 1e-6) beam_size = 18.1 ;
  if (abs(wavelength - 350e-6) < 1e-6) beam_size = 25.2 ;
  if (abs(wavelength - 500e-6) < 1e-6) beam_size = 36.9 ;

  if (is_void(beam_size)) {
    "ERROR: cannot obtain beam_size. Exiting." ;
    error ;
  }

  // Reads the fits file and pixel size
  fh = cfitsio_open(filename) ;
  image = cfitsio_read_image(fh) ;

  if (is_void(image)) {
      "ERROR : can't read fits file. Exiting" ;
      "Is it a single HDU file ?" ;
      error ;
  }

  dims = dimsof(image) ;
  if ((dimsof(image))(1) > 2) image = image(,,1) ; // only considering the 1st plane

  //if(is_void(image)) {
  //  fits_goto_hdu, fh, 2;
  //  fits_list, fh ;
  //  image = fits_read_array(fh) ;
  //  if (is_void(image)) {
  //    "ERROR : can't read fits file. Exiting" ;
  //    error ;
  //  }
  //}

  if (typeof(image(1)) == "double") image = float(image) ;

  CDELT1 = cfitsio_get(fh,"CDELT1") ;
  if (is_void(CDELT1)) CDELT1 = cfitsio_get(fh,"CD1_1") ;
  pix_size = abs(CDELT1) * 3600.; // en arcsec

  write, "Number of pixels", dimsof(image)(2), "x", dimsof(image)(3) ;
  write, "pixel_scale=", pix_size,"arcsec" ;


  // Conversion en MJy/sr
  BUNIT = cfitsio_get(fh,"BUNIT") ;
  if (is_void(BUNIT)) BUNIT = cfitsio_get(fh,"QTTY____") ;
  if (BUNIT == "MJy/sr") {
    write, "Data in MJy/sr, no conversion needed" ;
    conversion_factor = 1.0 ;
  } else {
    if (instrument == "PACS") {
      // Jy/pixel --> MJy/sr
      write, "Assuming PACS data is in Jy/pixel" ;
      write, "Converting Jy/pixel to MJy/sr";
      conversion_factor = 1.0 / (1e6 * (pix_size * SI.as)^2) ;
    } else if (instrument == "SPIRE") {
      // Jy/beam --> MJy/sr
      // 1.13309 * thetaA^2 is the integral of a Gaussian beam
      // where thetaA is the HPBW (Half-Power Beamwidth)
      write, "Assuming SPIRE data is in Jy/beam" ;
      write, "Converting Jy/beam to MJy/sr";
      conversion_factor = 1.0 / (1e6 * (beam_size * SI.as)^2) * pi/4./log(2.) ;
    }
  }
  write,  "conversion factor = ", conversion_factor ;
  cfitsio_close, fh ;

  image *= float(conversion_factor) ;

  ou = where(!ieee_test(image)) ;

  // remove median for speed
  write, "Computing image median ..." ;
  //med = median(image(ou)) ;
  write, "Skipping median of image ..." ;
  med = 0
  write, "measured median [MJy/sr]", med ;

  return image ;
 }


//*************************************************************************

func read_APEX(filename,instrument,&pix_size,&beam_size,&conversion_factor) {
/* DOCUMENT  read_APEX(filename,instrument,&pix_size,&beam_size,&conversion_factor)
   Reads a fits file filename containing APEX (SABOCA or LABOCA) data
   Convert the flux to Mjy/sr assuming PACS data is in Jy/beam

   instrument = "SABOCA" or "LABOCA"

   Returns:
     -- beam size in arcsec
     -- pixel size in arcsec
   SEE ALSO:
 */


  // Reads the fits file and pixel size
  image = cfitsioRead(filename, fh) ;

  if (is_void(image)) {
      "ERROR : can't read fits file. Exiting" ;
      "Is it a single HDU file ?" ;
      error ;
  }

  dims = dimsof(image) ;
  if ((dimsof(image))(1) > 2) image = image(,,1) ; // only considering the 1st plane

  if (typeof(image(1)) == "double") image = float(image) ;

  CDELT1 = cfitsio_get(fh,"CDELT1") ;
  if (is_void(CDELT1)) CDELT1 = cfitsio_get(fh,"CD1_1") ;
  pix_size = abs(CDELT1) * 3600.; // en arcsec

  write, "Number of pixels", dimsof(image)(2), "x", dimsof(image)(3) ;
  write, "pixel_scale=", pix_size,"arcsec" ;

  // Beam size (circular for APEX)
  beam_size = cfitsio_get(fh,"BMAJ") * 3600.; // en arcsec
  if (is_void(beam_size)) {
    "ERROR: cannot obtain beam_size. Exiting." ;
    error ;
  }
  write, "beam_size=", beam_size,"arcsec" ;

  // Conversion en MJy/sr
  // Jy/beam --> MJy/sr
  // 1.13309 * thetaA^2 is the integral of a Gaussian beam
  // where thetaA is the HPBW (Half-Power Beamwidth)
  BUNIT = cfitsio_get(fh,"BUNIT") ;
  if (BUNIT == "MJy/sr") {
    write, "Data in MJy/sr, no conversion needed" ;
    factor = 1.0 ;
  } else if (BUNIT == "Jy/beam") {
    write, "data is in Jy/beam" ;
  write, "Converting Jy/beam to MJy/sr";
  factor = 1.0 / (1e6 * (beam_size * SI.as)^2) * pi/4./log(2.) ;
  conversion_factor = factor ;
  } else {
    "ERROR: unknown BUNIT. Exiting." ;
    write, "BUNIT=", BUNIT ;
    error ;
  }
  write,  "conversion factor = ", conversion_factor ;

  return image * float(factor);
 }


//*************************************************************************

func fit_sed(flux,sigma,upper,wl,&T,&dens,&chi2,&T0,&dens0, &T_error,&dens_error_factor, logfit=,plot=,stdev=,correl=) {
/* DOCUMENT
     - Fit a map by a black body * power-law opacity law
     - Perform the Bayesian analysis of the different parameters

     Units:
     - wavelengths in m

     Si le flux en W.m^-2.Hz^-1, ie 1e26 Jy et _fact en m^2.kg^-1 * Msun / pc^2,
     dens0 est renvoye en Msun


   SEE ALSO:
 */

  extern _kappa, _kappa0, _l0, _fact ;
  extern tab_flux_ratio, tab_T ;

  //flux = fact * 10.^1 * bb_nu(30.,wl) * (1+0.1*random(5))  ; // TEST
  //flux = [-0.00148352,-0.00120954,0.34486,0.502697,0.55466]

  if (is_void(logfit)) logfit=0 ;
  if (is_void(plot)) plot=0 ;

  if (is_void(stdev)) stdev=0 ;
  if (is_void(correl)) correl=0 ;


  flux0 = flux ;
  sigma0 = sigma ;
  ou_upper = where(upper) ;
  ou_data = where(!upper) ; //points qui ne sont des upper limits


  if (plot) {
    fma ; lxy, 1, 1 ;

    if (numberof(ou_data)) plp, flux0(ou_data), wl(ou_data) * 1e6, symbol = 6, dy=sigma0(ou_data)  ;
    if (numberof(ou_upper)) plp, 3*sigma0(ou_upper), wl(ou_upper) * 1e6, symbol = 7 ;
    limits ; relimits_log ;
  }

  // Using color if we only have 2 data points
  if ( (numberof(flux) == 2) && (numberof(ou_data) == 2)) {
    flux_ratio = flux(2)/flux(1) ;
    T = interp(tab_T, tab_flux_ratio, flux_ratio) ;

    dens = flux(2) / (_fact(2) * bb_nu(T,wl(2)))  ;
    // Using proper fit if more than 3 data points
  } else if (numberof(ou_data) >= 3) {  // do the fitting if we have at least 1 point which is not an upper limit


    // initial guess:
    // - estimate the temperature from the wavelength where the SED peaks
    // - compute the density from the last point given the estimated temperature
    T0 = 25. * 100e-6/wl(flux(mxx)) ; // Loi de Wien
    l =  ou_data(0) ;
    dens0 = log10( flux(l) / (_fact(l) * bb_nu(T0,wl(l))) ) ;
    a = [T0,dens0] ;


    if (logfit==1) {
      if (anyof(flux(ou_data)) <= 0.) {
        T = -99 ;
        dens = -99 ;
        chi2 = -99 ;
        write, "ERROR: negative flux in SED" ;
        return [] ;
      }

      // Flux en log
      sigma(ou_data) = sigma(ou_data)/flux(ou_data) ; // erreur en relatif
      flux(ou_data) = log(flux(ou_data)) ;

      // Poids
      W = sigma * 0. ;
      W(ou_data) = 1./sigma(ou_data)^2 ;

      // Upper limits
      if (numberof(ou_upper)) W(ou_upper) = 0. ;

      // fit
      result= lmfit(log_black_body, wl, a, flux, W, deriv=1, stdev=stdev, monte_carlo=500, correl=correl) ;

    } else if (logfit==2) {
      if (anyof(flux(ou_data)) <= 0.) {
        T = -99 ;
        dens = -99 ;
        chi2 = -99 ;
        write, "ERROR: negative flux in SED" ;
        return [] ;
      }

      _is_upper_limit = array(0,numberof(flux)) ;

      // Flux en log
      sigma(ou_data) = sigma(ou_data)/flux(ou_data) ; // erreur en relatif
      flux(ou_data) = log(flux(ou_data)) ;

      // Poids
      W = sigma * 0. ;
      W(ou_data) = 1./sigma(ou_data)^2 ;

      // Upper limits
      if (numberof(ou_upper)) {
        W(ou_upper) = 1./sigma(ou_upper)^2 ;
        flux(ou_upper) = 0.0 ;
        _is_upper_limit(ou_upper) = 1 ;
      }

      X = [wl,_is_upper_limit] ;

      // fit
      result= lmfit(log_black_body_upper, X, a, flux, W, deriv=1, stdev=stdev, monte_carlo=0, correl=correl) ;

    } else {

      // Poids
      W = 1./sigma^2 ;
      // Upper limits
      if (numberof(ou_upper)) {
        flux(ou_upper) = 0.0 ;
      }
      // fit
      result = lmfit(black_body2, wl, a, flux, W, deriv=1, stdev=stdev, monte_carlo=0, correl=correl) ;
    }

    if (plot) {
      //*** pltr, swrite(i)+swrite(j), 0, 90 ;

      wl2 = spanl(50.,1000.,100) * 1e-6 ;
      plg, _kappa0 * (wl2/_l0)^-2 *  SI.Msun / ((distance*SI.pc)^2)  * bb_nu(a(1),wl2) * 10.^a(2), wl2* 1e6, color="red" ;
      //plg, _kappa0 * (wl2/_l0)^-2 *  SI.Msun / ((distance*SI.pc)^2)  * bb_nu(T0,wl2) * 10.^dens0, wl2* 1e6, color="blue", type=2 ;

      //plp, _kappa *  SI.Msun / ((distance*SI.pc)^2)  * bb_nu(a(1),wl) * 10.^a(2), wl* 1e6, color="red" ;

      //*** pltr, swrite(a(1))+swrite(10^a(2)*_density_factor), 0, 70 ;
    }

    // Save results
    T = a(1) ;
    dens = 10.^a(2) ;
    chi2 = result.chi2_last ;

    dens0 = 10^dens0 ;

    T_error = (*result.stdev)(1) ;
    dens_error_factor = 10.^((*result.correl)(2)) ;

    //write, "T", T, T_error ;


  } else {
    T = NaN ;
    dens = NaN ;
    chi2 = NaN ;

    T_error = NaN ;
    dens_error_factor = NaN ;
  }


}
//*************************************************************************


func black_body(x,a,&grad,deriv=) {
  // a(1) is temperature
  // 10^a(2) is density or mass
  extern _fact, cst_th ;

  wl= x ;
  T = a(1) ;

  if (a(1) < 0) {
    zero = x * 0. ;
    grad = [zero, zero] +1e-3;
    return zero ;
  }

  hc_lkT = cst_th/(T*wl) ;
  fact_exp = exp(min(hc_lkT,700.) ) ;

  if ((min(fact_exp) -1.0) > 1e-20) {
    BB = 2.*SI.h*SI.c *  max(1./ ((fact_exp -1.) * wl^3), 1e-200) ; // same as bb_nu
  } else {
    BB = 1e50 ;
  }
  BB = 10.^min(a(2),100) * _fact * BB ;

  // derivee par rapport a T, teste Ok numeriquement + maxima
  if ((min(fact_exp) -1.0) > 1e-20) {
    dBB_dT = BB * hc_lkT / T * (fact_exp / (fact_exp -1.));
  } else {
    dBB_dT = 1e50 ;
  }

  if (deriv) {
    grad = [dBB_dT, BB * log(10.)] + 1e-300;
  }

  return BB ;
}

//*************************************************************************
func black_body2(x,a,&grad,deriv=) {
  // a(1) is temperature
  // 10^a(2) is density or mass
  extern _fact, cst_th ;

  wl= x ;
  T = a(1) ;

  if (a(1) < 0) {
    zero = x * 0. ;
    grad = [zero, zero] +1e-3;
    return zero ;
  }

  hc_lkT = cst_th/(T*wl) ;
  fact_exp = exp(min(hc_lkT,700.) ) ;

  if ((min(fact_exp) -1.0) > 1e-20) {
    BB = 2.*SI.h*SI.c *  max(1./ ((fact_exp -1.) * wl^3), 1e-200) ; // same as bb_nu
  } else {
    BB = 1e50 ;
  }
  BB = 1.0 - exp(-10.^min(a(2),100) * _fact) * BB ;

  // derivee par rapport a T, teste Ok numeriquement + maxima
  if ((min(fact_exp) -1.0) > 1e-20) {
    dBB_dT = BB * hc_lkT / T * (fact_exp / (fact_exp -1.)) ;
  } else {
    dBB_dT = 1e50 ;
  }

  if (deriv) {
    grad = [dBB_dT, BB * log(10.)] + 1e-300;
  }

  return BB ;
}

//*************************************************************************

func log_black_body(x,a,&grad,deriv=) {
  // a(1) is temperature
  // a(2) is density
  extern _fact, cst_th ;

  wl= x ;
  T = a(1) ;

  if (a(1) < 0) {
    zero = x * 0. ;
    grad = [zero, zero] +1e-3;
    return zero ;
  }

  hc_lkT = cst_th/(T*wl) ;
  fact_exp = min(exp(min(hc_lkT,700.) ), 1e100) ;

  if ((min(fact_exp) -1.0) > 1e-20) {
    BB = 2.*SI.h*SI.c *  max(1./ ((fact_exp -1.) * wl^3), 1e-200) ; // same as bb_nu
  } else {
    BB = 1e50 ;
  }
  log_BB = log(10.^min(a(2),100) * _fact * BB) ;

  // derivee par rapport a T, teste Ok numeriquement + maxima
  if ((min(fact_exp) -1.0) > 1e-20) {
    dBB_dT = hc_lkT / T * fact_exp / (fact_exp -1.) ;
  } else {
    dBB_dT = 1e50 ;
  }

  if (deriv) {
    grad = [dBB_dT, log(10.)] + 1e-100;
  }

  return log_BB ;
}

//*************************************************************************

func log_black_body_upper(x,a,&grad,deriv=) {
  // a(1) is temperature
  // a(2) is density
  extern _fact, cst_th;

  wl = x(,1) ;
  _is_upper_limit = x(,2) ;

  T = a(1) ;

  ou_upper = where(_is_upper_limit) ;
  ou_data = where(!_is_upper_limit) ;

  log_BB = dBB_dT = array(0.,numberof(wl)) ;

  if (a(1) < 0) {
    zero = x * 0. ;
    grad = [zero, zero] +1e-3;
    return zero ;
  }


  hc_lkT = cst_th/(T*wl) ;
  fact_exp = min(exp(min(hc_lkT,700.) ), 1e100) ;

  if ((min(fact_exp) -1.0) > 1e-20) {
    BB = 2.*SI.h*SI.c *  max(1./ ((fact_exp -1.) * wl^3), 1e-200) ; // same as bb_nu
  } else {
    BB = 1e50 ;
  }
  BB = 10.^min(a(2),100) * _fact * BB ;
  log_BB(ou_data) = log(BB(ou_data)) ;
  if (numberof(ou_upper)) log_BB(ou_upper) = BB(ou_upper) ; // On ne renvoie pas le log

  // derivee par rapport a T, teste Ok numeriquement + maxima
  if ((min(fact_exp) -1.0) > 1e-20) {
    dBB_dT(ou_data) = hc_lkT(ou_data) / T * fact_exp(ou_data) / (fact_exp(ou_data) -1.) ;
    if (numberof(ou_upper)) {
      dBB_dT(ou_upper) = BB(ou_upper) * hc_lkT(ou_upper) / T * fact_exp(ou_upper) / (fact_exp(ou_upper) -1.) ;
    }
  } else {
    dBB_dT = 1e50 ;
  }

  if (deriv) {
    grad = [dBB_dT, log(10.)] + 1e-100;
  }

  return log_BB ;
}

//*************************************************************************

func xy2sky_file(filename,coord_filename) {
/* DOCUMENT xy2sky_file(filename,coord_filename)
   Write an ASCII file with the coordinates of the pixels in degrees
   SEE ALSO:
 */

  extern _wcs_bin_dir ;

  fh = cfitsio_open(filename) ;
  naxis = cfitsio_get(fh,"NAXIS") ;
  nx = cfitsio_get(fh,"NAXIS1") ;
  ny = cfitsio_get(fh,"NAXIS2") ;
  cfitsio_close, fh ;

  /*
  if (naxis != 2) {
    "Error: Fits file is not an image" ;
    "Exiting" ;
    return ;
  }
  */

  f=open("pixel_list.tmp","w") ;
  ix = indgen(nx)(,-:1:ny) ;
  iy = indgen(ny)(-:1:nx,:) ;
  write, f, ix, iy ;
  close, f ;

  // Extraction des coordonnees en degrees
  if (strchar(strchar(filename)(-2:0))=="gz") {
    filename2 = "tmp.fits" ;
    system, "gunzip -c "+filename+" > "+filename2
  } else {
    filename2 = filename ;
  }
  system, _wcs_bin_dir+"/xy2sky -j "+filename2+" @pixel_list.tmp | awk -F \" \" '{print $1 \" \" $2}' > "+coord_filename ;
  system, "rm -rf pixel_list.tmp tmp.fits" ;

  return ;
}


//*************************************************************************

func remap_fits(filename,nx,ny,coord_filename,image=) {
/* DOCUMENT
   Remap a fits file on the coordinates provided in a ASCII file
   The dims of the new coordinates is nx * ny
   Image KEYWORD is if you want to replace the image in the fits file
   must have the same dims as the fits file
   SEE ALSO:
 */

  extern _wcs_bin_dir ;

  // Verification image
  if (is_void(image)) {
    image=cfitsioRead(filename) ;
  } else {
    fh = cfitsio_open(filename);
    nim_x = cfitsio_get(fh, (key = "NAXIS1")) ;
    nim_y = cfitsio_get(fh, (key = "NAXIS2")) ;
    cfitsio_close, fh ;

    N=dimsof(image) ;
    if ( (nim_x != N(2)) || (nim_y != N(3)) ) {
      "*** Error in remap_fits: incorrect number of pixel in provided image" ;
      write, filename, " size", nim_x, nim_y ;
      write, N(2), N(3) ;
      error ;
      return [] ;
    }
  }

  // Calcul des indices pixels
  system, "rm -rf pixel_list.tmp" ;
  if (strchar(strchar(filename)(-2:0))=="gz") {
    filename2 = "tmp.fits" ;
    system, "gunzip -c "+filename+" > "+filename2
  } else {
    filename2 = filename ;
  }
  system, _wcs_bin_dir+"/sky2xy -j "+filename2+" @"+coord_filename+" | awk -F \" \" '{print $5 \" \" $6}' > pixel_list.tmp" ;
  system, "rm -rf tmp.fits" ;

  // Extraction des flux des pixels
  pixels=OpenASCII("pixel_list.tmp",valtype="float",prompt=0) ;

  N_pixels = numberof(pixels) ;
  //write, N_pixels, "were found"

  if (N_pixels != nx*ny) {
    "*** Error in remap_fits: incorrect number of pixel in coord map" ;
    write, N_pixels ;
    write, nx, ny, "->", nx*ny ;
    error ;
    return [] ;
  }


  interp = bilinear(image,pixels.X1,pixels.X2,minus_one=0,outside=0) ;

  //ou = where(interp < -1e10) ;
  //if (numberof(ou)) interp(ou) = NaN ;

  im = array(float,nx,ny) ;
  im(*) = interp(*) ;



  return im ;
}

//------------------------------------------------------------

func pixel_sed(win=,color=,pixel=) {
  // TODO : 9/11/2010 : bug le fit lineaire n'est pas le meme que dans le routine proncipale
  // pixels 329,358 et 328,358 pour m16


  extern map, map_sigma, map_upper, wl, distance, Temperature, density, chi2 ;

  if (is_void(color)) color="black" ;
  if (is_void(win)) win=window() ;

  if (is_void(pixel)) {

    window, win ;
    m = mouse(1) ;
    i = int(m(1)) ;
    j = int(m(2)) ;

    if (i==0 && (j==0)) {
      write, "Pixel error" ;
      return ;
    }
  } else {
    i = pixel(1) ;
    j = pixel(2) ;
  }

  write, "Pixel", i, j ;


  // define flux and error in the pixel
  flux = map(i,j,) ;
  sigma = map_sigma(i,j,) ;
  upper = map_upper(i,j,)

  write, "lambda (micron)   Flux (Jy)   Error (Jy)" ;
  write, wl*1e6, flux, sigma ;

  window, win+1; fma ; lxy, 1, 1 ;
  ou_data = where(!upper) ;
  ou_upper = where(upper) ;

  W = array(0.,numberof(flux)) ;

  if (numberof(ou_data)) plp, flux(ou_data), wl(ou_data) * 1e6, symbol = 6, dy =sigma(ou_data)  ;
  if (numberof(ou_upper)) plp, 3*sigma(ou_upper), wl(ou_upper) * 1e6, symbol = 7 ;
  limits ; relimits_log ;

  if (numberof(ou_data)) {


    T = Temperature(i,j) ;
    dens = density(i,j) ;
    write, "Map parameters:    T=", T, "dens=", dens, "chi2=", chi2(i,j) ;

    wl2 = spanl(50.,1000.,100) * 1e-6 ;
    plg, _kappa0 * (wl2/_l0)^-2 *  SI.Msun / ((distance*SI.pc)^2)  * bb_nu(T,wl2) * dens / _density_factor, wl2* 1e6, color=color ;

    /*

    // do the fitting (linear)
    if (numberof(ou_data) > 2) {
      W(ou_data) = 1.0/ sigma(ou_data)^2 ;

      // initial guess:
      // compute the density from the last point given a chosen temperature
      T0 = 10. ;//
      l =  ou_data(0) ;
      dens0 = log10( flux(l) / (_fact(l) * bb_nu(T0,wl(l))) ) ;
      a = [T0,dens0] ;

      result= lmfit(black_body, wl, a, flux, W, deriv=1, stdev=0, monte_carlo=0, correl=0) ;

      plg, _kappa0 * (wl2/_l0)^-2 *  SI.Msun / ((distance*SI.pc)^2)  * bb_nu(a(1),wl2) * 10^a(2), wl2* 1e6, color="green", type=2 ;

      write, "linear fit (green) T=", a(1), "dens=", 10^a(2) * _density_factor ;
    }

    // do the fitting (log)
    if (numberof(ou_data) > 2) {
      // initial guess:
      // compute the density from the last point given a chosen temperature
      T0 = 10. ;//
      l =  ou_data(0) ;
      dens0 = log10( flux(l) / (_fact(l) * bb_nu(T0,wl(l))) ) ;
      a = [T0,dens0] ;

      W(ou_data) = 1.0/ sigma(ou_data)^2 ;

      flux(ou_data) = log(flux(ou_data)) ;
      W = array(0.,numberof(flux)) ;
      W(ou_data) = 1.0 ;

      result= lmfit(log_black_body, wl, a, flux, W, deriv=0, stdev=0, monte_carlo=0, correl=0) ;

      plg, _kappa0 * (wl2/_l0)^-2 *  SI.Msun / ((distance*SI.pc)^2)  * bb_nu(a(1),wl2) * 10^a(2), wl2* 1e6, color="blue", type=2 ;


      write, "log fit (blue)    T=", a(1), "dens=", 10^a(2) * _density_factor ;
    }
    */
    window, win ;

    return ;
  }
}

//------------------------------------------------------------

func log_normal(x,a) {

  mu = a(1) ;
  sigma = a(2) ;
  fact = a(3) ;

  y = fact/(x*sigma*sqrt(2*pi)) * exp( -(log(x) - mu)^2/(2*sigma^2) ) ;
  return y ;
}

//------------------------------------------------------------

func read_scale_offset(file) {
  f = open(file) ;
  rdline, f, 17 ;

  instru = filter = unite = norm = array(string,8) ;
  wav = scale = offset = beam_surf_sr = beam_surf_arcsec2 = ratio = array(float,8) ;
  read, f, instru, filter, wav, unite, norm, scale, offset, beam_surf_sr, beam_surf_arcsec2, ratio ;
  close, f ;


  // PACS and SPIRE fluxes
  select = [1,3,4,5,6] ;
  select = [] ;

  offset = offset(select) ;
  scale = scale(select) ;

  /*
  // Test to find out what is ratio
  ratio = ratio(select) ;
  beam_surf_sr = beam_surf_sr(select) ;
  beam_surf_arcsec2 = beam_surf_arcsec2(select) ;
  // write, beam_surf_sr / beam_surf_arcsec2, "" ;
  write, scale * ratio * beam_surf_sr * 1e6, "" ; // this is equal to 1
  */


  return [offset,scale] ;
}


//------------------------------------------------------------

func read_scale_offset2(file,wl) {

  n_lines = int() ;
  f = popen("wc -l "+file,0) ;
  read, f, n_lines ;
  close, f ;


  f = open(file) ;
  rdline, f, 17 ;

  instru = filter = unite = norm = array(string,n_lines-17) ;
  wav = scale = offset = beam_surf_sr = beam_surf_arcsec2 = ratio = array(float,n_lines-17) ;
  read, f, instru, filter, wav, unite, norm, scale, offset, beam_surf_sr, beam_surf_arcsec2, ratio ;
  close, f ;

  Offset = Scale = array(float,numberof(wl)) ;
  for (i=1 ; i<=numberof(wl) ; i++) {
    ou = where(abs(wav - wl(i)*1e6) < 1e-3) ;
    if (numberof(ou) != 1) {
      write, "ERROR : scale and offset are not available (or multiple) for this wavelength", wl(i) ;
      error ;
    } else {
      Offset(i) = offset(ou) ;
      Scale(i) = scale(ou) ;
    }
  }

  return [Offset,Scale] ;
}

//------------------------------------------------------------

func read_median_offset(file) {

  f = open(file) ;
  rdline, f, 19 ;

  instru = filter = unite = norm = array(string,5) ;
  wav = scale = offset = beam_surf_sr = beam_surf_arcsec2 = ratio =  scale_error = offset_error = rchi2 = medianp = medianh = array(float,5) ;
  read, f, instru, filter, wav, unite, norm, scale, offset, beam_surf_sr, beam_surf_arcsec2, ratio, scale_error, offset_error, rchi2, medianp, medianh ;
  close, f ;

  return [medianp,medianh] ;
}

//------------------------------------------------------------

func read_median_offset2(file,wl) {

  n_lines = int() ;
  f = popen("wc -l "+file,0) ;
  read, f, n_lines ;
  close, f ;

  f = open(file) ;
  rdline, f, 19 ;

  instru = filter = unite = norm = array(string,n_lines-19) ;
  wav = scale = offset = beam_surf_sr = beam_surf_arcsec2 = ratio =  scale_error = offset_error = rchi2 = medianp = medianh = array(float,n_lines-19) ;
  read, f, instru, filter, wav, unite, norm, scale, offset, beam_surf_sr, beam_surf_arcsec2, ratio, scale_error, offset_error, rchi2, medianp, medianh ;
  close, f ;


  Medianp = Medianh = array(float,numberof(wl)) ;
  for (i=1 ; i<=numberof(wl) ; i++) {
    ou = where(abs(wav - wl(i)*1e6) < 1e-3) ;
    if (numberof(ou) != 1) {
      write, "ERROR : scale and offset are not available (or multiple) for this wavelength", wl(i) ;
      error ;
    } else {
      Medianp(i) = medianp(ou) ;
      Medianh(i) = medianh(ou) ;
    }
  }

  return [Medianp,Medianh] ;
}

//------------------------------------------------------------

func make_PDF(data,xtitle=,ytitle=,log=,color=,mask=,N_bins=,plot_error=,type=,width=,cumul=) {
/* DOCUMENT make_PDF(data,xtitle=,ytitle=,log=,color=,mask=,N_bins=,plot_error=,type=,width=,cumul=)

   SEE ALSO:
 */

  if (is_void(N_bins)) N_bins=200 ;
  if (is_void(log)) log=0 ;
  if (is_void(plot_error)) plot_error=0 ;
  if (is_void(color)) color="black" ;
  if (is_void(type)) type=1 ;
  if (is_void(cumul)) cumul=0 ;


  dims = dimsof(data) ;
  nx = dims(2) ; ny = dims(3) ;

  if (is_void(mask)) {
    ou = where( !ieee_test(data) & data > 0) ;
  } else {
    ou = where( mask & !ieee_test(data) & data > 0) ;
  }

  if (numberof(ou)) {

    // On redefinit x et y
    if (log) {
      y = log10(data(ou)) ;
    } else {
      y = data(ou) ;
    }

    threshold = 1e-3 ;
    bins = span(floor(min(y)),floor(max(y))+2,N_bins+1) ;
    x = bins(zcen) ; // bin center

    histoData = digitize(y,bins) ;
    N = histogram(histoData)  ;

    // Checking dimensions
    if(numberof(where(histoData==N_bins+1))!=0)
      N = N(:-1) ;
    while(numberof(N)<numberof(x))
      grow, N, 0.0 ;

    dy = 1.0 * sqrt(N) ;

    norme = sum(N) ;
    N = 1.0 * N / norme ;
    dy = dy / norme ;

    if (log) {
      x_plot=  10^x ;
    } else {
      x_plot = x ;
    }

    // info on the PDF
    write, " peak at =", x_plot(N(mxx)) ;
    write, " median =", median(data(ou)) ;
    write, " 1st quartile =", percentile(data(ou),0.25) ;
    write, " 3rd quartile =", percentile(data(ou),0.75) ;


    if (cumul) N = N(psum) ;

    plh, N, x_plot, color=color, type=type ;
    if (plot_error) plp, N, x_plot, dy = dy, symbol=0, fill=1, size=0.1, color=color ;

    limits ;
    limits, , , 0.5 * threshold * max(N), 1.2 * max(N) ;
    xytitles, xtitle, ytitle, [0.0,0.02] ;

  } else { // Polygon is empty
    write, "Polygon is empty, skipping PDF" ;
  }

}

//------------------------------------------------------------

func make_CMF(data,color=,N_bins=,cumul=) {

  if (is_void(color)) color="black" ;
  if (is_void(N_bins)) N_bins = 15 ;
  if (is_void(cumul)) cumul = 0 ;


  y = log10(data) ;

  bins = span(min(y),max(y)+1,N_bins+1) ;


  x_h = 10^bins  ;
  x = x_h / sqrt(x_h(2)/x_h(1)) ;// bin center

  if (!cumul) {
    N = histogram(digitize(y,bins),top=numberof(bins)) ;
    dy = sqrt(N) ;

    plh, N, x_h, color=color, just=3 ;
    plp, N, x, dy = 10 * dy, color=color, symbol=20, fill=1, size=0.5 ;

    ou = where(dy > 0) ;

    Y = N(ou) ;
    X = x(ou) ;
    W = 1./dy(ou)^2 ;

    a = [1,log10(Y(mxx)),100000.] ;
    result= lmfit(log_normal, X, a, Y, W, deriv=0, stdev=0, monte_carlo=0, correl=0) ;

    y = log_normal(X,a) ;

    x = spanl(1e-3,1e7,100) ;
    y = log_normal(x,a) ;
    plg, y, x, color=color ;

    write, "Log-normal fit: peak=", x(y(mxx)), "Msun, deviation=",a(2)/log(10)," in log10(M)" ;


    pltitle_height=14 ;
    xytitles, "Mass [M_sun]", "# of sources: !DN/!Dlog(M)" ;
    limits, , , 0.3, 1.5 * max(N) ;
  } else {
    plot_histo(y,x,cumul=1,color=color) ;
  }
}

//-----------------------------------------------------------

func calc_Lum(T,M) {

  extern _kappa0, _l0 ;

  lambda = spanl(1,10000,100) * 1e-6 ;

  alambda = lambda(2) / lambda(1) ;
  alambda = sqrt(alambda) - 1.0/sqrt(alambda) ;

  delta_lambda = alambda * lambda ;

  Lsun = SI.sigma * 5780.^4 * SI.Rsun^2 ;// Sans le 4pi pour etre conherent avec les autres luminosites

  kappa = _kappa0 * (lambda/_l0)^-2 ;


  BB = SI.Msun * M * bb(T,lambda) * kappa ;// 1 pour F_lambda
  L = sum(BB  * delta_lambda) / Lsun;

  return L ;
}

//------------------------------------------------------------

func neighbour_distance(clumps,stage,color=) {

  if (is_void(color)) color="black" ;

  N_clumps = numberof(clumps) ;
  x = clumps.XCO_A ;
  y = clumps.YCO_A ;

  neighbour = array(float,N_clumps) ;
  for(i=1 ; i<=N_clumps ; i++) {
    x0 = x(i) ; y0 = y(i) ;
    dist2 = (x-x0)^2 + (y-y0)^2 ;

    neighbour(i) = min(dist2(where(dist2 > 0))) ;
  }
  neighbour = sqrt(neighbour) ;

  //difference between prestellar and protostars
  neighbour1 = neighbour(where(stage==1)) ;
  neighbour2 = neighbour(where(stage==2)) ;

  if (numberof(neighbour2) > 1) write, "median neighbour distance for protostars = ",median(neighbour2), "arcsec" ;
  if (numberof(neighbour1) > 1) write, "median neighbour distance for starless = ",median(neighbour1), "arcsec" ;


  //plot nearest neighbour comparison pre/proto stellar
  bins = span(1.,1300,1000) ;
  plot_histo, neighbour2,bins,color=color, type=1, cumul=1 ;
  plot_histo, neighbour1,bins,color=color, type=2, cumul=1 ;
  limits, 1, 1000, 0, 1 ;
  plg, [0.9,0.9],[2,4],   width=3, type=1, color="black";
  plt, "protostars", 5, 0.9, tosys=1 ;

  plg, [0.85,0.85],[2,4],   width=3, type=2, color="black";
  plt, "starless", 5, 0.85, tosys=1 ;

  xytitles, "Distance to closest neighbour (arcsec)", "Cumulative distribution" ;

  /*
    pdf, region+"_neighbour.pdf" ;
    eps, region+"_neighbour.eps" ;
    system, "mv -f "+region+"_neighbour.pdf "+dir ;
    system, "mv -f "+region+"_neighbour.eps "+dir ;
  */

  if (numberof(neighbour1) * numberof(neighnour2) ) print, "KS test =", kstwo(neighbour1,neighbour2) ;

  return neighbour ;
}

//------------------------------------------------------------

func write_catalog(cat_filename, cat) {

  f = open(cat_filename,"w") ;
  write, f, "#   NO   XCO_P  YCO_P   XCO_A   YCO_A   SNR_SSX  FS  SN_RATIO1  FXP_BEST1  FXP_ERRO1  FXT_BEST1  FXT_ERRO1 AFWHM1 BFWHM1 THEPA1  SN_RATIO2  FXP_BEST2  FXP_ERRO2  FXT_BEST2  FXT_ERRO2 AFWHM2 BFWHM2 THEPA2  SN_RATIO3  FXP_BEST3  FXP_ERRO3  FXT_BEST3  FXT_ERRO3 AFWHM3 BFWHM3 THEPA3  SN_RATIO4  FXP_BEST4  FXP_ERRO4  FXT_BEST4  FXT_ERRO4 AFWHM4 BFWHM4 THEPA4  SN_RATIO5  FXP_BEST5  FXP_ERRO5  FXT_BEST5  FXT_ERRO5 AFWHM5 BFWHM5 THEPA5" ;
  write, f, indgen(numberof(cat)), cat.XCO_P, cat.YCO_P, cat.XCO_A, cat.YCO_A, cat.SNR_SSX, cat.FS, cat.SN_RATIO1, cat.FXP_BEST1, cat.FXP_ERRO1, cat.FXT_BEST1, cat.FXT_ERRO1,  cat.AFWHM1, cat.BFWHM1, cat.THEPA1, cat.SN_RATIO2, cat.FXP_BEST2, cat.FXP_ERRO2, cat.FXT_BEST2, cat.FXT_ERRO2, cat.AFWHM2, cat.BFWHM2, cat.THEPA2, cat.SN_RATIO3, cat.FXP_BEST3, cat.FXP_ERRO3, cat.FXT_BEST3, cat.FXT_ERRO3, cat.AFWHM3, cat.BFWHM3, cat.THEPA3, cat.SN_RATIO4, cat.FXP_BEST4, cat.FXP_ERRO4, cat.FXT_BEST4, cat.FXT_ERRO4, cat.AFWHM4, cat.BFWHM4, cat.THEPA4, cat.SN_RATIO5, cat.FXP_BEST5, cat.FXP_ERRO5, cat.FXT_BEST5, cat.FXT_ERRO5, cat.AFWHM5, cat.BFWHM5, cat.THEPA5 ;
  close, f ;
}

//------------------------------------------------------------

func remap_cat(cat,ref_file) {

  // New call to OpenASCII to update the structure format ;
  ref = OpenASCII(ref_file,prompt=0) ;

  new_cat = array(OPEN_ASCII_TMP,numberof(cat)) ;

  new_cat.FS       = cat.FS       ;
  //new_cat.SNR_SSX  = cat.SNR_SSX  ;
  new_cat.XCO_A    = cat.XCO_A    ;
  new_cat.XCO_P    = cat.XCO_P    ;
  new_cat.YCO_A    = cat.YCO_A    ;
  new_cat.YCO_P    = cat.YCO_P    ;


  new_cat.AFWHM2      =   cat.AFWHM1       ;
  new_cat.AFWHM3      =   cat.AFWHM2       ;
  new_cat.AFWHM4      =   cat.AFWHM3       ;
  new_cat.AFWHM5      =   cat.AFWHM4       ;
  new_cat.BFWHM2      =   cat.BFWHM1       ;
  new_cat.BFWHM3      =   cat.BFWHM2       ;
  new_cat.BFWHM4      =   cat.BFWHM3       ;
  new_cat.BFWHM5      =   cat.BFWHM4       ;
  new_cat.FXP_BEST2   =   cat.FXP_BEST1    ;
  new_cat.FXP_BEST3   =   cat.FXP_BEST2    ;
  new_cat.FXP_BEST4   =   cat.FXP_BEST3    ;
  new_cat.FXP_BEST5   =   cat.FXP_BEST4    ;
  new_cat.FXP_ERRO2   =   cat.FXP_ERRO1    ;
  new_cat.FXP_ERRO3   =   cat.FXP_ERRO2    ;
  new_cat.FXP_ERRO4   =   cat.FXP_ERRO3    ;
  new_cat.FXP_ERRO5   =   cat.FXP_ERRO4    ;
  new_cat.FXT_BEST2   =   cat.FXT_BEST1    ;
  new_cat.FXT_BEST3   =   cat.FXT_BEST2    ;
  new_cat.FXT_BEST4   =   cat.FXT_BEST3    ;
  new_cat.FXT_BEST5   =   cat.FXT_BEST4    ;
  new_cat.FXT_ERRO2   =   cat.FXT_ERRO1    ;
  new_cat.FXT_ERRO3   =   cat.FXT_ERRO2    ;
  new_cat.FXT_ERRO4   =   cat.FXT_ERRO3    ;
  new_cat.FXT_ERRO5   =   cat.FXT_ERRO4    ;
  new_cat.SN_RATIO2   =   cat.SN_RATIO1    ;
  new_cat.SN_RATIO3   =   cat.SN_RATIO2    ;
  new_cat.SN_RATIO4   =   cat.SN_RATIO3    ;
  new_cat.SN_RATIO5   =   cat.SN_RATIO4    ;
  new_cat.THEPA2      =   cat.THEPA1       ;
  new_cat.THEPA3      =   cat.THEPA2       ;
  new_cat.THEPA4      =   cat.THEPA3       ;
  new_cat.THEPA5      =   cat.THEPA4       ;

  return new_cat ;
}

//------------------------------------------------------------

func secToHMS(time)
/* DOCUMENT secToHMS(time)
 * Convert from time (float in sec) to string in the form hh:mm:ss.s
 * AUTHOR : F.Rigaut, June 13 2002.
 * Modify to account for < 0    C. Pinte
 * SEE ALSO:
 */
{
  if (anyof(time < 0)) {
    offset = 1 ;
    time += 100 * 3600 ;
  }

  lt    = long(time);
  hh    = lt/3600;
  lt    = lt-hh*3600;
  mm    = lt/60;
  sec   = float(time)-hh*3600-mm*60;

  if (!is_void(offset)) hh = hh - 100 ;

  ts    = swrite(hh,mm,sec,format="%02i:%02i:%04.2f");
  return ts;
}

//------------------------------------------------------------

func make_power_spectrum(data,pix_size,xtitle=,ytitle=,log=,color=,type=,width=,mask=) {
/* DOCUMENT make_power_spectrum(data,pix_size,xtitle=,ytitle=,log=,color=,type=,width=,mask=)
     pix_size en arcsec
   SEE ALSO:
 */

  if (is_void(color)) color="black" ;
  if (is_void(type)) type=1 ;


  // Amplitude

  if (is_void(mask)) {
    D = data ;
  } else {
    D = data * mask ;
  }
  NaN_to_zero, D ;

  A = abs(fft(D,1))^2 ;
  A = roll(A) ;

  dims = dimsof(A) ;
  nx = dims(2) ; ny = dims(3) ;

  center = indexof(A,A(*)(mxx)) ;


  x = (indgen(nx)-center(1))(,-:1:ny) ; y = (indgen(ny)-center(2))(-:1:nx,) ;
  dist = abs(x,y) ;

  bins = span(1,min(nx,ny)/2.,100) ;
  n_bins = numberof(bins)-1 ;

  PS = PS_std = array(double,n_bins) ;

  for (i=1 ; i<=n_bins ; i++) {
    ou = where( (dist > bins(i)) & (dist < bins(i+1)) ) ;
    if (numberof(ou) > 1) {
      PS(i) = median(A(ou)) ;
      PS_std(i) = (A(ou))(rms) ;
      //plp, A(ou), dist(ou), color=color ;
    }
  }

  pix_fft = 1.0/pix_size * 60;
  size = max(nx,ny) ;
  dx_fft = pix_fft/size ;

  X = bins(zcen) * dx_fft ;

  plg, PS, X, color=color ;
  pldj, X, PS * (1 +PS_std/PS), X, PS / (1 +PS_std/PS), color=color ;
  lxy, 1, 1 ;
  if (!is_void(xtitle)) xytitles, xtitle, "" ;
  if (!is_void(ytitle)) xytitles, "", ytitle ;
  limits, 0.01, , ,  ;

  //plg, 3.6e54 * (x/0.026)^(-2.65), x, color="red"

}

//------------------------------------------------------------

func set_model_info(filename, P, conf) {

  fh = cfitsio_open(filename,"a") ;

  cfitsio_set, fh,"Author", P.author ;
  cfitsio_set, fh,"Herschel_program", P.Herschel_program ;
  cfitsio_set, fh,"scan_speed", P.scan_speed ;
  cfitsio_set, fh,"region_name", P.region_name ;
  cfitsio_set, fh,"map_resolution", P.beam_size ;
  cfitsio_set, fh,"distance", P.distance ;
  cfitsio_set, fh,"parameter_file", P.parameter_file ;
  cfitsio_set, fh,"log_fit", P.logfit ;
  cfitsio_set, fh,"min_noise", P.min_noise ;
  cfitsio_set, fh,"Planck_correct", P.Planck_correct ;
  if (P.Planck_correct) cfitsio_set, fh,"Planck_correct_method", P.Planck_correct_method ;

  cfitsio_set, fh,"nu0", P.nu0/1e9, "[GHz]" ;
  cfitsio_set, fh,"kappa0", P.kappa0 * 10, "[cm^2/g]" ;


  for (i=1 ; i<=numberof(conf) ; i++) {
    cfitsio_set, fh,"file_"+swrite(format="%1i",i), conf(i).file ;
    //if (conf(i).ref == 1) cfitsio_set, fh, "ref_file", "file_"+swrite(format="%1i",i) ;
  }

  cfitsio_close, fh ;
}

//------------------------------------------------------------

func save_fits(filename, image, WCS, P, conf, unit=) {

  fh = cfitsio_open(filename,"w",overwrite=1);
  cfitsio_create_img, fh, cfitsio_bitpix_of(image), dimsof(image);
  dim       = dimsof(image);
  fpixels   = array(1,dim(1));
  nelements = numberof(image);
  cfitsio_write_pix, fh, fpixels, nelements, image;

  if (WCS.type==1) {
    cfitsio_set, fh,"CDELT1", WCS.CDELT1 ;
    cfitsio_set, fh,"CDELT2", WCS.CDELT2 ;
  } else {
    cfitsio_set, fh,"CD1_1", WCS.CD_1 ;
    cfitsio_set, fh,"CD1_2", WCS.CD_2 ;
    cfitsio_set, fh,"CD2_1", WCS.CD_1 ;
    cfitsio_set, fh,"CD2_2", WCS.CD_2 ;
  }

  cfitsio_set, fh,"CRVAL1", WCS.CRVAL1 ;
  cfitsio_set, fh,"CRVAL2", WCS.CRVAL2 ;
  cfitsio_set, fh,"CRPIX1", WCS.CRPIX1 ;
  cfitsio_set, fh,"CRPIX2", WCS.CRPIX2 ;
  cfitsio_set, fh,"CTYPE1", WCS.CTYPE1 ;
  cfitsio_set, fh,"CTYPE2", WCS.CTYPE2 ;
  cfitsio_set, fh,"CROTA1", WCS.CROTA1 ;
  cfitsio_set, fh,"CROTA2", WCS.CROTA2 ;
  cfitsio_set, fh,"EQUINOX", WCS.EQUINOX ;

  cfitsio_set, fh,"Author", P.author ;
  cfitsio_set, fh,"Software", "Yorick column density / temperature map maker" ;
  cfitsio_set, fh,"Herschel_program", P.Herschel_program ;
  cfitsio_set, fh,"scan_speed", P.scan_speed ;
  cfitsio_set, fh,"region_name", P.region_name ;
  cfitsio_set, fh,"map_resolution", P.beam_size ;
  cfitsio_set, fh,"distance", P.distance ;
  cfitsio_set, fh,"parameter_file", P.parameter_file ;
  cfitsio_set, fh,"log_fit", P.logfit ;
  cfitsio_set, fh,"min_noise", P.min_noise ;
  cfitsio_set, fh,"Planck_correct", P.Planck_correct ;
  if (P.Planck_correct) cfitsio_set, fh,"Planck_correct_method", P.Planck_correct_method ;

  cfitsio_set, fh,"nu0", P.nu0/1e9, "[GHz]" ;
  cfitsio_set, fh,"kappa0", P.kappa0 * 10, "[cm^2/g]" ;

  if (!is_void(unit)) cfitsio_set, fh,"BUNIT", unit ;

  for (i=1 ; i<=numberof(conf) ; i++) {
    cfitsio_set, fh,"file_"+swrite(format="%1i",i), conf(i).file ;
    //if (conf(i).ref == 1) cfitsio_set, fh, "ref_file", "file_"+swrite(format="%1i",i) ;
  }

  cfitsio_write_date,fh;
  cfitsio_close,fh;
  return 1;
}
