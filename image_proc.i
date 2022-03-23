func convol2d(A,b) {
/* DOCUMENT
     Convolution 2d
   SEE ALSO:
 */
  dimsA=dimsof(A);
  dimsb=dimsof(b);


  if ((dimsA(1)!=2) || dimsb(1)!=2) {
    write("Error : convol2d works only in images");
    return [];
  }

  nx=dimsA(2);
  ny=dimsA(3);
  nx_mask=dimsb(2);
  ny_mask=dimsb(3);

  xstart = (nx_mask+1)/2;
  ystart = (ny_mask+1)/2;

  C=array(double,nx,ny);

  A_0 =array(double,nx+nx_mask,ny+ny_mask) ;
  /*
  for (i=1+nx_mask/2 ; i<=nx-nx_mask/2 ; i++) {
    for (j=1+ny_mask/2 ; j<=ny-ny_mask/2 ; j++) {
      C(i,j) = sum(A(i-nx_mask/2:i+nx_mask/2,j-ny_mask/2:j+ny_mask/2)*b)
        }
  }
  */

  A_0(xstart+1:xstart+nx,xstart+1:xstart+ny) = A;
  for (i=1 ; i<=nx ; i++) {
    for (j=1 ; j<=ny ; j++) {
      C(i,j) = sum(A_0(i+1:i+nx_mask,j+1:j+ny_mask)*b);
      }
    }


  return C;
}


func convol2df(im,psf){

  psf2 = im * 0.;

  dx=dimsof(im)(2);
  dy=dimsof(im)(3);

  dx2=dimsof(psf)(2)  ;
  dy2=dimsof(psf)(2)  ;

  startx = (dx-dx2)/2;
  starty = (dy-dy2)/2;

  //psf2(startx+1:startx+dx2,starty+1:starty+dy2) = psf;
  psf2(startx+1:startx+dx2,starty+1:starty+dy2) = psf;

  fim = fft(im,1);
  fpsf = fft(psf2,1);

  fim = fim * fpsf;

  im = roll(fft(fim,-1));
  im = im.re / (dx*dy);

  return im
}


func convol2df_noise(im,psf,noise=){

  if (is_void(noise)) noise=0.0 ;

  psf2 = im * 0.;

  dx=dimsof(im)(2);
  dy=dimsof(im)(3);

  dx2=dimsof(psf)(2)  ;
  dy2=dimsof(psf)(2)  ;

  startx = (dx-dx2)/2;
  starty = (dy-dy2)/2;

  //psf2(startx+1:startx+dx2,starty+1:starty+dy2) = psf;
  psf2(startx+1:startx+dx2,starty+1:starty+dy2) = psf;

  fim = fft(im,1) ;
  fpsf = fft(psf2,1);

  fim = fim * fpsf;

  fim += abs(fim) * noise * random(dx,dy) ;


  im = roll(fft(fim,-1));
  im = abs(im) / (dx*dy);

  return im
}


func gauss_kernel_1d(npix, sigma) {
/* DOCUMENT
     Cree un noyau gaussien dont l'intégrale vaut 1
     Marche pour npix pair et impair
     sigma en pixel
     SEE ALSO:
 */
  if (sigma < 1.0e-30) sigma=1.0e-30;

  centre = npix/2. + 0.5;

  x = indgen(npix) - centre;

  tmp = exp(-0.5*x^2/sigma^2);

  return tmp/sum(tmp);
}


func gauss_kernel(npix, sigma) {
/* DOCUMENT
     Cree un noyau gaussien dont l'intégrale vaut 1
     Marche pour npix pair et impair
     sigma en pixel
     SEE ALSO:
 */
  if (sigma < 1.0e-30) sigma=1.0e-30;

  centre = npix/2. + 0.5;

  tmp1 = indgen(npix)(,-:1:npix) - centre;
  tmp2 = indgen(npix)(-:1:npix,) - centre;
  pixel_map = [tmp1, tmp2];

  dist2 =  pixel_map(,,1)^2+pixel_map(,,2)^2;

  tmp = exp(-0.5*dist2/sigma^2);

  return tmp/sum(tmp);
}


func gauss_kernel2(npix, sigma_x, sigma_y, PA) {
/* DOCUMENT
    gauss_kernel2(npix, sigma_x, sigma_y, PA)
    sigma_x et y en pix PA en degres
    PA defini depuis la direction verticale dans le sens horaire
    Cree un noyau gaussien dont l'intégrale vaut 1
     Marche pour npix pair et impair
   SEE ALSO:
 */
  if (sigma_x < 1.0e-30) sigma_x=1.0e-30;
  if (sigma_y < 1.0e-30) sigma_y=1.0e-30;

  PA = PA * pi / 180.;

  centre = npix/2. + 0.5;

  tmp1 = indgen(npix)(,-:1:npix) - centre;
  tmp2 = indgen(npix)(-:1:npix,) - centre;

  x = tmp1 * cos(PA) + tmp2 * sin(PA) ;
  y = tmp1 * sin(PA) - tmp2 * cos(PA) ;

  x = x / sigma_x ;
  y = y / sigma_y ;

  dist2 =  x^2 + y^2 ;

  tmp = exp(-0.5*dist2);

  return tmp/sum(tmp);

}

include, "bessel.i";

func airy(x) {
  /* DOCUMENT
     fonction d'Airy
     1er zero en 1.22
  SEE ALSO:
  */

  y=array(double,dimsof(x));

  xx = x*pi ;

  zero =  where(abs(xx)<1.0e-10);
  not_zero = where(abs(xx)>1.0e-10);

  if (numberof(zero)) y(zero) = 1.0;
  y(not_zero) = 4.0*(bessj1(xx(not_zero))/xx(not_zero))^2;

  return y;

}

func airy_kernel(npix, sigma) {

/* DOCUMENT
     Cree un noyau gaussien dont l'intégrale vaut 1
     Marche pour npix pair et impair
   SEE ALSO:
 */
  centre = npix/2. + 0.5;

  tmp1 = indgen(npix)(,-:1:npix) - centre;
  tmp2 = indgen(npix)(-:1:npix,) - centre;
  pixel_map = [tmp1, tmp2];

  dist =  sqrt(pixel_map(,,1)^2+pixel_map(,,2)^2);

  tmp = airy(dist/sqrt(sigma));

  return tmp/sum(tmp);

}

func deriv(A) {
/* DOCUMENT
     Derivée d'un tableau 1d
   SEE ALSO:
 */

  dims=dimsof(A);
  if (dims(1) != 1) {
    write("Error : not 1d table");
  }

  n=dims(2);

  B=array(double,n);
  B(2:n-1) = 0.5*(A(3:n)-A(1:n-2));

  return B;
}

func deriv2(A) {
/* DOCUMENT
     Derivée seconde d'un tableau 1d
   SEE ALSO:
 */

  dims=dimsof(A);
  if (dims(1) != 1) {
    write("Error : not 1d table");
  }

  n=dims(2);

  B=array(double,n);
  B(2:n-1) = A(3:n)+A(1:n-2)-2*A(2:n-1);

  return B;
}


func get_pix(win) {
/* DOCUMENT
     Coordonnées des pixels cliqués
     Click droit pour quitter
   SEE ALSO:
 */
  if (is_void(win)) win=0;
  window, win;

  x=[];
  y=[];

  do {
    m=mouse();
    grow, x, m(1);
    grow, y, m(2);
  } while(m(-1)==1);

  return [x,y];

}


func telescope_beam(lambda,D,verbose=) {
/* DOCUMENT  telescope_beam(lambda,D,verbose=)
   - lambda and D in m
   - return FWHM  in as
   SEE ALSO:
 */

  if (is_void(verbose)) verbose=0 ;

  lambda_D = lambda/D * 648000 / pi ;

  x=span(-3,3,1000000) ;
  y = airy(x) ;

  if (verbose) {
    write, "resolution (1st zero of Airy):", 1.22 * lambda_D, "arcsec" ;
    window, 1 ; fma ;
    plg, y, x * lambda_D, color="red" ;
  }

  // Calcule la FWHM et sigma d'une finction d'Airy
  // pour fitter une gaussienne dessus
  // 1er zero fct airy 1.22 * lambda/D
  // http://en.wikipedia.org/wiki/Airy_disc

  ou = where(y > 0.5) ;
  FWHM1 = (x(max(ou)) - x(min(ou))) ;
  FWHM2 = (x(max(ou)+1) - x(min(ou)-1)) ;
  FWHM = 0.5 * (FWHM1 + FWHM2) * lambda_D;

  if (verbose) write, "FWHM", FWHM, "arcsec" ;

  sigma =  FWHM / (2.*sqrt(2*log(2.)) ) ;
  //sigma = 0.42 ; // meilleur d'apres-wikipedia

  if (verbose) {
    write, "sigma", sigma, "arcsec";
    plg, exp(-x^2/(2*(sigma/lambda_D)^2)), x * lambda_D   ;
    xytitles, "Dist (arcsec)", "Beam" ;
  }

  // --> pour une gaussienne, la FHWM est presque egale a lambda/D
  // FWHM = 0.989 * lambda/D ;
  // sigma = 0.42 * lambda/D ;

  return FWHM ;
}


func FWHM_airy(n) {

  x=span(-10,10,n) ;
  y = airy(x) ;


}

func median_image(img,n) {
/* DOCUMENT median_image(img,N)
     Apply a median filter on img on a box N x N
     N must be an odd number
     if (N==0) --> return img ;
   SEE ALSO:
 */

  if (n==0) {
    return img ;
  } else {
    write, "Median filter image ..." ;
    dims = dimsof(img) ;
    nx = dims(2) ;
    ny = dims(3) ;

    img2 = array(double,nx,ny) ;

    dy = dx = n/2 ;

    // NaN
    mask = ieee_test(img) ;
    ou = where(mask) ;
    if (numberof(ou)) img2(ou) = img(ou) ;

    // Pixel avec flux
    mask = !mask ;

    // Milieu de l'image
    for (i=1+dx ; i<= nx-dx ; i++) {
      for (j=1+dy ; j<= ny-dy ; j++) {
        if (mask(i,j)) {
          square = img(i-dx:i+dx,j-dy:j+dy) ;
          ou = where(!ieee_test(square)) ;
          img2(i,j) = median(square(ou)) ;
        }
      }
    }

    // Bord de l'image
    img2(1:1+dx-1,:) = img(1:1+dx-1,:) ;
    img2(nx-dx+1:nx,:) = img(nx-dx+1:nx,:) ;
    img2(:,1:1+dy-1) = img(:,1:1+dy-1) ;
    img2(:,ny-dy+1:ny) = img(:,ny-dy+1:ny) ;

    write, "Done" ;
    return img2 ;
  }
}
