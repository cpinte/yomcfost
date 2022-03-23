func precompute_fft_speed {

  N = 2740 ;
  speed = array(float,N) ;

  write, "Analysing fft speed ..." ;
  for (i=1 ; i<= N ; i++) {
    write, i, "/", N ; 
    im = array(float,i,i) ;
    tic ;
    fim = fft(im, 1) ;
    speed(i) = tac() ;
    write, "speed=", speed(i) ;
  }
  
  cfitsWrite, "init/fft_speed.fits.gz",  speed ;
}


func fft_best_size(N) {

  fft_speed = cfitsRead("~/yorick/init/fft_speed.fits.gz") ;

  N_extra = 30 ;
  speed = array(float,N_extra) ;
  
  for (i=1 ; i<= N_extra ; i++) {
    speed(i) = fft_speed(N+2*(i-1)) ;
  }
  ou = speed(mnx) ;

  return N+2*(ou-1) ;  
}

  

// TEST 
/*
n = 1301 ;
im = array(float,n,n) ;
im(*) = 1 ;
im = gauss_kernel(n, n/100.) ;
n2 = fft_best_size(n) ;

write, n, n2 ; 

tic

n2 = 1323 ;
im2 = array(float,n2,n2) ; n_extra = (n2 -n)/2 ;
im2(n_extra+1:-n_extra,n_extra+1:-n_extra) = im ;
tac()

tic ; fim = fft(im, 1) ; tac() ;
tic ; fim2= fft(im2, 1) ; tac() ;

window, 10 ; plim, im2 ;

window, 11 ; fma ; 
plg, abs(fim(,1)), span(0,1,n) ;
plg, abs(fim2(,1)), span(0,1,n2), color="red" ;
limits, 0, 0.05 ; 
*/
