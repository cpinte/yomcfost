func relimits(rapx,rapy) {
  if (is_void(rapx)) rapx=8;
  if (is_void(rapy)) rapy=rapx;
  x=limits()(2)-limits()(1);
  y=limits()(4)-limits()(3);
  limits, (limits()(1)-x*rapx/100.0) , (limits()(2)+x*rapx/100.0) , (limits()(3)-y*rapy/100.0) , (limits()(4)+y*rapy/100.0);
}

func relimits_log(rapx,rapy) {
  if (is_void(rapx)) rapx=50;
  if (is_void(rapy)) rapy=50;
  limits, (limits()(1)*rapx*0.01) , (limits()(2)/(rapx*0.01)) , (limits()(3)*rapy*0.01) , (limits()(4)/(rapy*0.01));
}


func density(x,y,xgrid,ygrid,hx,hy) {
/* DOCUMENT
     convertie un nuage de points en densite par convolution
   SEE ALSO:
 */
  if(is_void(hy)) hy=hx

  if (anyof(dimsof(x)!=dimsof(y))) {
    write("Error : dim");
    return;
  }

  if (anyof(dimsof(xgrid)!=dimsof(ygrid))) {
    write("Error : dim");
    return;
  }

  z= exp(-0.5*((xgrid(-,) - x)^2/hx^2+(ygrid(-,) - y)^2/hy^2))(sum,,);

  return z;
}


func indexof(x,w) {
/* DOCUMENT  indexof(x,w)
     converti les indices 1d w de x en 2d (marche en nd)
   SEE ALSO:
 */
 if (!is_array(w)) return [];
 d= dimsof(x);
 n= d(1);
 if (!n) return w;  /* catcall for passing a scalar */
 d= d(2:);
  o= orgsof(x)(2:);
  w2= w(-:1:n,);
  w-= o(1);
  for (i=1 ; i<=n ; i++) {
    w2(i,)= w%d(i) + o(i);
    w/= d(i);
  }
  return w2;
}


func conv_ind(x,w) {
/* DOCUMENT
     converti les indices 1d w de x en nd
   SEE ALSO:
 */
 if (!is_array(w)) return [];
 d= dimsof(x);
 n= d(1);
 n
 if (!n) return w;  /* catcall for passing a scalar */

 d= d(2:);
 o= orgsof(x)(2:);
 w2= w(-:1:n,);
 w-= o(1);
 for (i=1 ; i<=n ; i++) {
   w2(i,)= w%d(i) + o(i);
   w/= d(i);
 }
 return w2;
}



func fitgauss(x, y) {
  /* DOCUMENT fitgauss(x, y)

      return [A,x0,sigma] where y = A exp( -(x-x0)^2/(2*sigma^2))

      X and Y must be arrays of 3 or more points.

    SEE ALSO: poly2
  */

  ou = where(y> 0);
  y=y(ou);
  x=x(ou);

  pol = fitpoly(2,x,log(y));

  sig2=-1.0/(2.0*pol(3));
  x0 = pol(2)*sig2;
  A = exp(pol(1) + x0*x0/(2*sig2));

  return [A,x0,sqrt(sig2)];

}

func corr(f1,f2) {
  dim1=dimsof(f1);
  dim2=dimsof(f1);
  dim=dim1(2)+dim2(2)-1;
  f = array(double,);

  for (i=1 ; i<=dim ; i++) {
    f(i)= f() ;
  }

}

func correl(f1,f2) {
  ff1 = fft(f1,1);
  ff2 = fft(f2,2);
  ff = ff1 * conj(ff2);

  f= fft(ff1,-1);
}


func xinv {l=limits() ; limits, l(2), l(1);}
func yinv {l=limits() ; limits, , ,l(4), l(3);}
func hhmm(value) {gridxy, hhmm=value;}

func lsdir2(expr) {
  tab_out = [];
  f=popen("ls "+expr, 0);
  for (i=0;line=rdline(f); i++)
    _,tab_out,line;
  if (_PROMPT==2) write,format="%d files",i;
  close, f;
  return tab_out;
}

func strmatch2 (str, match) {
  size = strlen(match);
  test = array(1 , dimsof(str));
  for (i=1; i<=size; i++)
    test = test & strmatch(str, strpart(match , i:i));

  return test;
}

func prmat(mat,t=,format=)
{
/* DOCUMENT prmat,matrix,format=,t=
	print matrices neatly. format= is format string for elements.
	If t= is used then does a transpose before printing.
 */
	typ=typeof(mat);
	if(is_void(format)){
		w=where(typ==["complex","float","double","short", "int", "long","string","char"]);
		if(numberof(w)!=1) error,"matrix is not of a recognized type";
		format=["%g+%gi","%6.3f ","%6.3f ","%6d ","%6d ","%6d ","%s ","%c "](w(1));
		}
	d=dimsof(mat);
	if(!is_void(t)) mat=transpose(mat);
	if(typ=="complex") {
		for(i=1;i<=d(2);i++){
			for(j=1;j<=d(3);j++){
				s=swrite(format=format,mat(i,j).re,mat(i,j).im);
				write,format="%10s ",s;
				}
			write,format="%s","\n";
			}
		}
	else {
		for(i=1;i<=d(2);i++){
			for(j=1;j<=d(3);j++) write,format=format,mat(i,j);
			write,format="%s","\n";
			}
		}
}


func undersample2d(a, nsub, &sigma, op=) {

/* DOCUMENT : undersample2d(a, nsub, op=)
      undersample2d(a, nsub, &sigma, op=)
      op=avg par defaut
   SEE ALSO:
 */
  nsub= int(nsub);
  if (nsub==1) "WARNING : nsub=1, cannot estimate uncertainties";

  if (is_void(op)) op=avg;
  tfunc = is_func(op);
  dim = dimsof(a);
  Nx = dim(2) / nsub + ((dim(2)/float(nsub)- dim(2)/nsub)>1e-30 ? 1 :0);
  Ny = dim(3) / nsub + ((dim(3)/float(nsub)- dim(3)/nsub)>1e-30 ? 1 :0);
  out = sigma = array(structof(a), Nx, Ny);
  N = Nx*Ny;

  for (i=1; i<=Nx ; i++) {
    for (j=1; j<=Ny ; j++) {
      i1 = (i-1)*nsub+1;
      i2 = i*nsub; i2 = (i2<=dim(2) ? i2 : dim(2));

      j1 = (j-1)*nsub+1;
      j2 = j*nsub; j2 = (j2<=dim(3) ? j2 : dim(3));

      if (tfunc) {
        out(i,j) = op(a(i1:i2,j1:j2)(*));
        sigma(i,j) = stddev(a(i1:i2,j1:j2))*nsub ;
      } else {
        out(i,j) = a(i1:i2,j1:j2)(*)(op);
        sigma(i,j) = stddev(a(i1:i2,j1:j2))*nsub ;
      }
    }
  }

  return out;
}


func convol1d(A,b) {
/* DOCUMENT
     Convolution 1d
   SEE ALSO:
 */
  dimsA=dimsof(A);
  dimsb=dimsof(b);


  if ((dimsA(1)!=1) || dimsb(1)!=1) {
    write("Error : convol2d works with 1d tables");
    return [];
  }

  nx=dimsA(2);
  nx_mask=dimsb(2);

  C=array(double,nx);

  C=A;
  b=b*1.0/sum(b);


  A_0 =array(double,nx+nx_mask) ;
    /*  for (i=1+nx_mask/2 ; i<=nx-nx_mask/2 ; i++) {
    for (j=1+ny_mask/2 ; j<=ny-ny_mask/2 ; j++) {
      C(i,j) = sum(A(i-nx_mask/2:i+nx_mask/2,j-ny_mask/2:j+ny_mask/2)*b)
        }
        }*/

  A_0(nx_mask/2+1:nx_mask/2+nx) = A;

  for (i=1 ; i<=nx ; i++) {
    C(i,j) = sum(A_0(i+1:i+nx_mask)*b);
  }

  return C;
}

c=3.0e8;
hp = 6.626e-34;
kb = 1.38e-23;
sigma = 5.6697e-8;
cst_th = hp*c/kb ;

func bb(T,lambda) {
/* DOCUMENT  bb(T,lambda)
   B_lambda(lambda) en W.m^-2.m^-1
   T en K
   lambda en m
   SEE ALSO:
 */
  extern cst_th ;
  wl=lambda;
  if (dimsof(T)(1)==0) {
      return 2.*hp*SI.c^2 *  max(1./ ( ((exp(min(cst_th/(T*wl),700.)) -1.)+1.e-30) * (wl^5)), 1e-200) ;
      //  cst_wl = cst_th/(T*wl)
      //return  max(1./ ( wl^5 * ( exp(min(cst_wl,700.)) -1.)+1.0e-30), 1.0e-20) ;
    } else {
      return  2.*hp*SI.c^2 * max(1./ ( ((exp(min(cst_th/(T(-,)*wl),700.)) -1.)+1.e-30) * (wl^5)), 1e-200) ;
    }
  // La somme * dlambda vaut bien  SI.sigma * T^4 / pi
}


func bb_nu_old(T,lambda,Jy=) {
/* DOCUMENT  bb_nu(T,lambda)
   B_nu(lambda) en W.m^-2.Hz^-1 or Jy if Jy=1
   T en K
   lambda en m
   This is the version of the routine used by Tracey's code
   SEE ALSO:
 */
  //F_nu = lambda^2 / c * F_lambda

  extern cst_th ;

  if (is_void(Jy)) Jy=0 ;

  if (Jy) {
    factor = 1e26 ;
  } else {
    factor = 1 ;
  }

  wl=lambda;
  if (dimsof(T)(1)==0) {
    return factor * 2.*hp*SI.c *  max(1./ ( ((exp(min(cst_th/(T*wl),700.)) -1.)+1.e-30) * (wl^3)), 1e-200) ;
  } else {
    return  factor * 2.*hp*SI.c * max(1./ ( ((exp(min(cst_th/(T(-,)*wl),700.)) -1.)+1.e-30) * (wl^3)), 1e-200) ;
  }
}


func bb_nu(T,lambda,Jy=) {
/* DOCUMENT  bb_nu(T,lambda)
   B_nu(lambda) en W.m^-2.Hz^-1 or Jy if Jy=1
   T en K
   lambda en m
   update 29/08/2016 to check surface brightness temperature
   SEE ALSO:
 */
  extern cst_th ;

  if (is_void(Jy)) Jy=0 ;

  if (Jy) {
    factor = 1e26 ;
  } else {
    factor = 1 ;
  }

  nu=SI.c/lambda;

  //if (dimsof(T)(1)==0) {
  return factor * 2.*SI.h/SI.c^2 *  max(nu^3/ ( (exp(min(SI.h * nu / (SI.kB * T),700.)) -1.)+1.e-30), 1e-200) ;
  //} else {
  //return factor * 2.*hp/SI.c^2 *  max(nu^3/ ( (exp(min(SI.h * nu(,-) / (SI.kB * T(-,)),700.)) -1.)+1.e-30), 1e-200) ;
  //}
}


func bb_etoile(T,R,dist,lamb){
/* DOCUMENT bb_etoile(T,R,dist,lamb)
     renvoie lambda.F_lambda en W.m-2
     T en K
     R en rayon solaire
     dist en pc
     lamb en m
  SEE ALSO:
 */
  // Verif
  //tmp = bb(T,lamb)
  //alambda = sqrt(lamb(2)/lamb(1)) ;
  //delta_lambda = lamb * alambda - lamb / alambda ;
  // A = SI.sigma * T^4 / pi ;
  //sum(tmp * delta_lambda ) / A ;// plus precis avec un nombre de points faibles : OK
  //sum(tmp(zcen) * lamb(dif) ) / A ;

  return pi * bb(T,lamb) * lamb * ((R*SI.Rsun)/(dist * SI.pc))^2

}


func Tbrightness_to_Flux(T, wl, BMAJ, BMIN) {
  /* DOCUMENT
     Convert brightness temperature to Flux density in Jy

     T [K]
     wl [m]
     BMAJ, BMIN in [deg], ie as in fits header

     Flux in Jy

     SEE ALSO:

  */
  conversion_factor = (BMIN * BMAJ * (3600*SI.as)^2 * pi/4./log(2.)) ;

  return bb_nu(T,wl,Jy=1) * conversion_factor ;
}


func Flux_to_Tbrightness(F, wl, BMAJ, BMIN) {
  /* DOCUMENT
     Convert Flux density in Jy/beam to brightness temperature [K]
     Flux [Jy]
     wl [m]
     BMAJ, BMIN in [deg], ie as in fits header

     T [K]

     SEE ALSO:
  */

  nu = SI.c/wl ;

  factor = 1e26 ;
  conversion_factor = (BMIN * BMAJ * (3600*SI.as)^2 * pi/4./log(2.)) ;

  //F = factor *  conversion_factor * 2.*hp/SI.c^2 * nu^3/ (exp(SI.h * nu / (SI.kB * T)) -1.) ;

  exp_m1 = factor *  conversion_factor * 2.*hp/SI.c^2 * nu^3/F ;
  hnu_kT =  log(max(exp_m1,1e-10) + 1) ;

  T = SI.h * nu / (hnu_kT * SI.kB) ;

  return T ;
}


func stddev(x) {

/* DOCUMENT
     computes standard deviation of an array
     similar as IDL function

     different from rms : 1/(N-1) and not 1/N
   SEE ALSO:
 */
  if (numberof(x) > 1)
    return sqrt(sum((x-avg(x))^2)/(numberof(x)-1.));
  else
    return 0;
}

func retourne(A) {
  n=numberof(A);
  return A(int(span(n,1,n)));
}

struct WCS_struct {
  int type ;

  double CDELT1 ;
  double CDELT2 ;
  double CRVAL1 ;
  double CRVAL2 ;
  double CRPIX1 ;
  double CRPIX2 ;
  double CD1_1 ;
  double CD1_2 ;
  double CD2_1 ;
  double CD2_2 ;
  double CROTA1 ;
  double CROTA2 ;
  string CTYPE1 ;
  string CTYPE2 ;
  double EQUINOX ;
} ;

/*
CD1_1 = CDELT1 * cos (CROTA2)
CD1_2 = -CDELT2 * sin (CROTA2)
CD2_1 = CDELT1 * sin (CROTA2)
CD2_2 = CDELT2 * cos (CROTA2)
*/

func get_WCS(filename,hdu=) {

  WCS = array(WCS_struct) ;
  buffer = float() ;

  fh = cfitsio_open(filename,"r") ;

  CRVAL1 = cfitsio_get(fh,"CRVAL1") ;
  if (typeof(CRVAL1)=="string") {
    sread, CRVAL1, buffer ; WCS.CRVAL1 = buffer ;
   } else {
    WCS.CRVAL1 = CRVAL1 ;
  }
  CRVAL2 = cfitsio_get(fh,"CRVAL2") ;
  if (typeof(CRVAL2)=="string") {
    sread, CRVAL2, buffer ; WCS.CRVAL2 = buffer ;
   } else {
    WCS.CRVAL2 = CRVAL2 ;
  }
  CRPIX1 = cfitsio_get(fh,"CRPIX1") ;
  if (typeof(CRPIX1)=="string") {
    sread, CRPIX1, buffer ; WCS.CRPIX1 = buffer ;
   } else {
    WCS.CRPIX1 = CRPIX1 ;
  }

  CRPIX2 = cfitsio_get(fh,"CRPIX2") ;
  if (typeof(CRPIX2)=="string") {
    sread, CRPIX2, buffer ; WCS.CRPIX2 = buffer ;
   } else {
    WCS.CRPIX2 = CRPIX2 ;
  }

  WCS.CTYPE1 = cfitsio_get(fh,"CTYPE1") ;
  WCS.CTYPE2 = cfitsio_get(fh,"CTYPE2") ;

  test = cfitsio_get(fh,"CDELT1") ;
  if (!is_void(test)) {
    WCS.type = 1 ;
    WCS.CDELT1 = test ;
    WCS.CDELT2 = cfitsio_get(fh,"CDELT2") ;
  } else {
    WCS.type = 2 ;
    WCS.CD1_1 =  cfitsio_get(fh,"CD1_1") ;
    WCS.CD1_2 =  cfitsio_get(fh,"CD1_2") ;
    WCS.CD2_1 =  cfitsio_get(fh,"CD2_1") ;
    WCS.CD2_2 =  cfitsio_get(fh,"CD2_2") ;
  }

  CROTA1 = cfitsio_get(fh,"CROTA1") ; if (is_void(CROTA1)) CROTA1 = 0.0 ;
  CROTA2 = cfitsio_get(fh,"CROTA2") ; if (is_void(CROTA2)) CROTA2 = 0.0 ;
  equinox = cfitsio_get(fh,"EQUINOX") ; if (is_void(equinox)) equinox = cfitsio_get(fh,"EPOCH") ;
  if (typeof(CROTA1)=="string") {
    sread, CROTA1, buffer ; WCS.CROTA1 = buffer ;
    sread, CROTA2, buffer ; WCS.CROTA2 = buffer ;
    sread, equinox, buffer ; WCS.EQUINOX = buffer ;
  } else {
    WCS.CROTA1 = CROTA1 ;
    WCS.CROTA2 = CROTA2 ;
    if (!is_void(equinox)) WCS.EQUINOX = equinox ;
  }

  cfitsio_close, fh ;

  return WCS ;
}

func set_WCS(filename,WCS) {

  fh = cfitsio_open(filename,"a") ;

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

  cfitsio_close, fh ;

}

func get_WCS2(filename) {
  // use fits.i instead of cfitsio
  WCS = array(WCS_struct) ;

  fh = fits_open(filename,"r") ;

  test = fits_get(fh,"CDELT1") ;
  if (!is_void(test)) {
    WCS.type = 1 ;
    WCS.CDELT1 = test ;
    WCS.CDELT2 = fits_get(fh,"CDELT2") ;
  } else {
    WCS.type = 2 ;
    WCS.CD1_1 =  fits_get(fh,"CD1_1") ;
    WCS.CD1_2 =  fits_get(fh,"CD1_2") ;
    WCS.CD2_1 =  fits_get(fh,"CD2_1") ;
    WCS.CD2_2 =  fits_get(fh,"CD2_2") ;
  }

  WCS.CRVAL1 = fits_get(fh,"CRVAL1") ;
  WCS.CRVAL2 = fits_get(fh,"CRVAL2") ;
  WCS.CRPIX1 = fits_get(fh,"CRPIX1") ;
  WCS.CRPIX2 = fits_get(fh,"CRPIX2") ;
  WCS.CTYPE1 = fits_get(fh,"CTYPE1") ;
  WCS.CTYPE2 = fits_get(fh,"CTYPE2") ;
  equinox = fits_get(fh,"EQUINOX") ;
  if (is_void(equinox)) equinox = fits_get(fh,"EPOCH") ;
  WCS.EQUINOX = equinox ;
  fits_close, fh ;

  return WCS ;
}

func set_WCS2(filename,WCS) {
  // use fits.i instead of cfitsio
  filename2 = filename+".tmp" ;

  system, "rm -rf "+filename2 ;

  fh = fits_open(filename,"r") ;
  fh2 = fits_open(filename2,"w") ;

  keywords_name = fits_get_keywords(fh) ;
  n_key = numberof(keywords_name) ;

  // Copy the header from the initial file
  for (i=1 ; i<= n_key ; i++) {
    fits_set, fh2, keywords_name(i), fits_get(fh, keywords_name(i)) ;
  }

  if (WCS.type==1) {
    fits_set, fh2,"CDELT1", WCS.CDELT1 ;
    fits_set, fh2,"CDELT2", WCS.CDELT2 ;
  } else {
    fits_set, fh2,"CD1_1", WCS.CD1_1 ;
    fits_set, fh2,"CD1_2", WCS.CD1_2 ;
    fits_set, fh2,"CD2_1", WCS.CD2_1 ;
    fits_set, fh2,"CD2_2", WCS.CD2_2 ;
  }
  fits_set, fh2,"CRVAL1", WCS.CRVAL1 ;
  fits_set, fh2,"CRVAL2", WCS.CRVAL2 ;
  fits_set, fh2,"CRPIX1", WCS.CRPIX1 ;
  fits_set, fh2,"CRPIX2", WCS.CRPIX2 ;
  fits_set, fh2,"CTYPE1", WCS.CTYPE1 ;
  fits_set, fh2,"CTYPE2", WCS.CTYPE2 ;
  fits_set, fh2,"EQUINOX", WCS.EQUINOX ;

  fits_write_header, fh2;

  // Copy the data from the initial file
  data = fits_read_array(fh) ;
  fits_write_array, fh2, data ;

  fits_close, fh ;
  fits_close, fh2 ;

  system, "mv -f "+filename2+" "+filename ;

}

//*************************************************************************

NaN = array(float,1) ; ieee_set, NaN,2 ; NaN = NaN(1) ;

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

  ou = where(ieee_test(img)) ;

  if (numberof(ou)) {
    img(ou) = 0. ;
  }
}

//*************************************************************************

func percentile(x, percent, which)
/* DOCUMENT percentile(x, percent)
         or percentile(x, percent, which)
     returns the median of the array X.  The search for the median takes
     place along the dimension of X specified by WHICH.  WHICH defaults
     to 1, meaning the first index of X.  The median function returns an
     array with one fewer dimension than its argument X (the WHICH
     dimension of X is missing in the result), in exact analogy with
     rank reducing index range functions.  If dimsof(X)(WHICH) is
     odd, the result will have the same data type as X; if even, the
     result will be a float or a double, since the median is defined
     as the arithmetic mean between the two central values in that
     case.
   SEE ALSO: sort
 */
{
  if (is_void(which)) which= 1;
  list= sort(x, which);
  dims= dimsof(x);
  if (which<1) which= dims(1)-which;
  n= dims(1+which);
  odd= n%2;
  n *= percent; n = int(n) ;   /* index with half above, half below... */
  n+= 1;         /* ...corrected for 1-origin */
  stride= 1;
  for (i=1 ; i<which ; i++) stride*= dims(1+i);
  ldims= dims(1)-which+1;
  /**/ local l;
  reshape, l, &list, long, stride, grow(ldims, dims(1+which:));
  lm= l(,n,..);
  if (which<dims(1)) dims(1+which:-1)= dims(2+which:0);
  --dims(1);
  reshape, lm, long, dims;
  xm= x(lm);
  if (!odd) {     /* even length dimensions have more complicated median */
    reshape, lm;  /* undefine the LValue lm so following define works */
    lm= l(,n-1,..);
    reshape, lm, long, dims;
    xm= 0.5f*(xm+x(lm));
  }
  return xm;
}

//*************************************************************************

func int_random(N,n=) {
/* DOCUMENT
     Return an uniformly distributed integer between 1 and N
   SEE ALSO:
 */

  r = int( random(n) * N + 1.0)  ;
  return r ;
}

//*************************************************************************

func Gauss_random(n) {

  r = array(double,n) ;

  for (i=1 ; i<=n ; i=i+2) {
      do {
      x = 2.0 * random(2) - 1.0;
      w = x(1) * x(1) + x(2) * x(2);
    } while ( w >= 1.0 );

    w = sqrt( (-2.0 * log(w)) / w );

    r(i)   = x(1) * w;
    if (i+1 <= n) r(i+1) = x(2) * w;
  }
  return r ;
}


//*************************************************************************

func popen(command, mode) {

  //write, "WARNING: popen does not work, using temporary file instead" ;

  system, "rm -rf yomcfost.tmp ; "+command+" > yomcfost.tmp " ;
  if (mode) {
    f = open("yomcfost.tmp","w") ;
  } else {
    f = open("yomcfost.tmp","r") ;
  }
  return f ;
}
