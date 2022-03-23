struct polygon {
  int n_points ;
  float x(100) ;
  float y(100) ;
}

//------------------------------------------------------------

func PolygonMask(map, polygon) {

  nx = dimsof(map)(2) ;
  ny = dimsof(map)(3) ;

  // Mask
  mask = array(0,nx,ny) ;

  // Elimate pixels not around the polygon
  imin = int( floor(min(polygon.x(1:polygon.n_points))) )  ;
  imax = int( floor(max(polygon.x(1:polygon.n_points)))+1 );
  jmin = int( floor(min(polygon.y(1:polygon.n_points))) ) ;
  jmax = int( floor(max(polygon.y(1:polygon.n_points)))+1 ) ;


  // Loop over remaining pixels
  for (i=imin ; i<imax ; i++) {
    for (j=jmin ; j<jmax ; j++) {
      mask(i,j) = inside_polygon(i,j,polygon) ;
    } // j
  } // i

  return mask ;
}

//------------------------------------------------------------

func inside_polygon(x,y,polygon) {
  /* DOCUMENT inside_polygon(x,y,polygon)
     compute if a point is inside a polygon or not
     SEE ALSO: PolygonMask
  */

  // Loop over polygon's corners
  counter = 0 ;
  for (ip=1; ip<=polygon.n_points ; ip++) {
    p1x = polygon.x(ip);
    p1y = polygon.y(ip);
    p2x = polygon.x(ip%polygon.n_points + 1);
    p2y = polygon.y(ip%polygon.n_points + 1);

    if ( (y > min(p1y,p2y)) & (y <= max(p1y,p2y)) ) {
      if (x <= max(p1x,p2x)) {
        if (p1y != p2y) {
          xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x;
          if (p1x == p2x || x <= xinters) counter = !counter;
        }
      } // test on x
    } // test on y
  } // loop over polygon points

  return counter
}

//------------------------------------------------------------

func read_polygons(polygon_file) {

  n_polygons = int() ;
  system, "grep polygon "+polygon_file+" | sed s/,/' '/g  | sed s/'polygon('// | sed s/') # color=black'// > polygons.txt" ;
  system, "wc -l polygons.txt > lignes.txt" ;
  f = open("lignes.txt","r") ;
  read, f, n_polygons ;
  close, f ;

  polygons = array(polygon,n_polygons) ;
  f = open("polygons.txt","r") ;
  for (i=1;i<=n_polygons;i++) {
    value=array(float,100) ;
    s = rdline(f)  ;
    sread, s, value ;
    value = value(where(value > 0)) ;
    N = numberof(value) ;
    iy = indgen(N/2) * 2;
    ix = iy -1 ;

    polygons(i).n_points = N/2 ;
    polygons(i).x(1:N/2) = value(ix) ;
    polygons(i).y(1:N/2)  = value(iy) ;
  }
  close, f ;
  return polygons ;
}

//------------------------------------------------------------

func remap_polygon(polygon,initial_fits_file,new_fits_file) {

  coord_filename = "polygons.coord" ;

  // Convert to RA & dec (degrees)
  f=open("pixel_list.tmp","w") ;
  write, f, polygon.x(1:polygon.n_points), polygon.y(1:polygon.n_points) ;
  close, f ;
  system, _wcs_bin_dir+"/xy2sky -j "+initial_fits_file+" @pixel_list.tmp | awk -F \" \" '{print $1 \" \" $2}' > "+coord_filename ;
  system, "rm -rf pixel_list.tmp" ;

  // Convert to the new pixel units
  system, _wcs_bin_dir+"/sky2xy -j "+new_fits_file+" @"+coord_filename+" | awk -F \" \" '{print $5 \" \" $6}' > pixel_list_new.tmp" ;
  pixels=OpenASCII("pixel_list_new.tmp",valtype="float",prompt=0) ;
  system, "rm -rf pixel_list_new.tmp "+coord_filename ;

  // Create new polygon
  new_polygon = polygon() ;
  new_polygon.n_points = polygon.n_points ;
  new_polygon.x(1:polygon.n_points) = pixels.X1 ;
  new_polygon.y(1:polygon.n_points) = pixels.X2 ;

  return new_polygon ;
}

//------------------------------------------------------------

func plot_polygon(polygon,color=,width=,type=) {

  if (is_void(color)) color="black" ;
  if (is_void(width)) width=1 ;
  if (is_void(type)) type=1 ;

  X = polygon.x(1:polygon.n_points) ; grow, X, polygon.x(1) ;
  Y = polygon.y(1:polygon.n_points) ; grow, Y, polygon.y(1) ;


  plg, Y, X, color=color ;
}

//------------------------------------------------------------

/*
polygons = read_polygons("getsources/ds9.reg") ;



Mask = array(0,dimsof(map)) ;
for (i=1 ; i<=numberof(polygons) ; i++) {
  print, "Polygon", i ;
  polygon = polygons(i) ;
  map = cfitsRead("getsources/Vela_density_zo_wo70um.fits") ;
  mask = PolygonMask(map, polygon) ;
  Mask = Mask | mask ;
 }

// Verif
window, 1 ; fma ;
plim, Mask  ;
for (i=1 ; i<=numberof(polygons) ; i++) {
  polygon = polygons(i) ;
  plg, polygon.y(1:polygon.n_points), polygon.x(1:polygon.n_points), color="green" ;
 }


window, 2; fma ;
polygon = polygons(1) ;
plg, polygon.y(1:polygon.n_points), polygon.x(1:polygon.n_points), color="green" ;

new_polygon = remap_region(polygons(1),"getsources/Vela_density_zo_wo70um.fits","getsources/spire350_bot_2_thirds.resamp.m.fits") ;
plg, new_polygon.y(1:polygon.n_points), new_polygon.x(1:polygon.n_points), color="red" ;
*/
