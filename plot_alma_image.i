func plot_alma_image(filename,win=,dx=,x_beam=, y_beam=) {

  if (is_void(win)) win=0 ;
  window, win, style="boxed.gs" ;

  fh = cfitsio_open(filename) ;
  image = cfitsio_read_image(fh) ;

  CDELT1 = cfitsio_get(fh,"CDELT1") ;
  CDELT2 = cfitsio_get(fh,"CDELT2") ;
  pixel_scale = abs(CDELT1) * 3600;
  distance = cfitsio_get(fh,"distance") ;

  a = cfitsio_get(fh,"BMAJ",comments) *3600. ;
  b = cfitsio_get(fh,"BMIN",comments) *3600.;
  theta = cfitsio_get(fh,"BPA",comments) ;

  write, "BEAM ", a, "\" x", b, "\"  PA=", theta ;

  if (!is_void(theta)) theta *= pi/180 ;

  img   = cfitsio_read_image(fh);
  cfitsio_close, fh ;

  // Plot image
  nx = dimsof(img)(2) ;
  ny = nx ;
  img_size = nx * pixel_scale ;
  img = img(,,1) ;


  if (is_void(cmax)) {
    Cmax = max(img) ;
  } else {
    Cmax = cmax ;
  }

  if (is_void(cmin)) {
    Cmin = min(img) ;
  } else {
    Cmin = cmin ;
  }
  pli, img, -img_size/2, -img_size/2, img_size/2, img_size/2, cmin=Cmin, cmax=Cmax ;

  // Plot beam
  if (!is_void(a)) {
    phi = span(0.,2*pi,360);
    x0 = -img_size/2 + max(nx / 10 * pixel_scale, 0.6*max(a,b)) ;
    y0 = x0 ;

    if (!is_void(x_beam)) x0 = x_beam ;
    if (!is_void(y_beam)) y0 = y_beam ;

    x= - a/2* cos(phi) * sin(theta) + b/2* sin(phi) * cos(theta)  ;
    y= + b/2* sin(phi) * sin(theta) + a/2* cos(phi) * cos(theta);

    for (l=1 ; l<=100 ; l++) {
      plg, y0 + 0.01*l*y, x0+ 0.01*l*x, width=1, color="white" ;
    }
  }
  xytitles, "RA offset (arcsec)", "DEC offset (sercsec)", [0.01,0.02];

}
