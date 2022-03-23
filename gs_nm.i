
/* Author : LeBouquin, jean-baptiste.lebouquin@obs.ujf-grenoble.fr
 * $Log: gs_nm.i,v $
 * Revision 1.4  2004/07/10 12:41:18  jblebou
 * Modify the limits_nm function. Now limits_nm without argument will set
 * the same limits at all the system of the window.
 *
 * Revision 1.3  2004/07/09 14:26:07  jblebou
 * global revision of all the files :
 *  - replace the old file name (toto_v1.i) by new one (tot.i)
 *  - add some author name
 *  - clean the CVS history
 *
 * Revision 1.2  2004/07/09 10:27:46  fmillour
 *
 */

if(_PRINT) write,"#include \"gs_nm.i\"";
require,"style.i";

/*---------------------------------------------------------------------------
This package contain some function to plot in N by M viewport :

 - gs_nm
 - plTitle
 - pltitle_nm
 - xytitles_nm
 - limits_nm
 - range_nm
 - logxy_nm
 - rangex_nm

---------------------------------------------------------------------------*/
func gs_nm(win, n, m, d= ,dx=, dy=, fy=, fx=, rx=, ry=, V=, style=, file=, square=, landscape=)
/* DOCUMENT gs_nm, win, n, m,
                   d= ,dx=, dy=, fy=, fx=,
                   rx=, ry=, V=, style=, file=, square=

   ***********
   arguments :
   
   'win' is the number of the window to style (default is current window)
   'n' is the number of column (default is n of the current window)
   'm' is the number of line (default is m of the current window)

   ********************
   optional arguments :
   
   V= : global viewport containing all the systems
   style= : default system type, valid gist file, default is "boxed.gs"
   square= : 1, force the systems to be square
   
   fx= : 1/0 set/unset the x-ticks                       (default [1,1,1,..])
   fy= : 1/0, set/unset the y-ticks                      (default [1,1,1,..])
   rx= : relative x-width of the different column        (default [1,1,1...])
   ry= : relative y-width of the different lines         (default [1,1,1...])
   
   dx= : x offset between the column
   dy= : y offset between the lines
   d=  : equivalent to dy=dx=d, default is 0.01

   file= : optional file output name (Gist style sheet format)

   ******
   call :
   
   gs_nm, win, n, m    : set a n.x.m window
   gs_nm, win, nwin    : set a n.x.m window with n.m>nwin
   gs_nm, win          : set a window with the same n.x.m than the curent one
   gs_nm,, dx=1.       : change the dx of the current window without destroy
                         the plot !

   **********
   history :
   - 2004/03/15 : LeBouquin, release 1.0

   SEE ALSO: gs_nm, limits_nm, range_nm, rangex_nm, logx_nm, xytitles_nm
             pltitle_nm, plTitle
 */
{
  local n,m,nn,mm,num,Dvx,Dvy,land,systems,leg,cleg;
  if(is_void(V) &&  is_void(landscape))      V = [0.12,0.73,0.1,0.9];
  if(is_void(V) && !is_void(landscape))      V = [0.1,0.9,0.12,0.73];
  if(is_void(style))  style = "boxed.gs";
  if(is_void(win))    win = window();
  if(is_void(d))      d = 0.01;
  
  if(is_void(n) && is_void(m)) read_nm,win,n,m;
  if(is_void(n) && !is_void(m)) read_nm,win,n;
  if(!is_void(n) && is_void(m)) {m = int(ceil(n/ceil(sqrt(n)))); n=int(ceil(sqrt(n))); n0=n;n=m;m=n0;}
    
  if(is_void(dx)) dx = double(d); if(n==1) dx=[dx(1)]; else {if(numberof(dx)==1) dx=dx(1)(-:1:n-1); if(numberof(dx)>n-1) dx=dx(1:n-1); if((num=numberof(dx))<n-1) grow,dx,dx(0)(-:1:n-1-num);}
  if(is_void(dy)) dy = double(d); if(m==1) dy=[dy(1)]; else {if(numberof(dy)==1) dy=dy(1)(-:1:m-1); if(numberof(dy)>m-1) dy=dy(1:m-1); if((num=numberof(dy))<m-1) grow,dy,dy(0)(-:1:m-1-num);}
  if(is_void(rx)) rx = 1.;        if(n==1) rx=[rx(1)]; else {if(numberof(rx)==1) rx=rx(1)(-:1:n);   if(numberof(rx)>n) rx=rx(1:n);     if((num=numberof(rx))<n)   rx = grow(double(rx(1:-1)),double(rx(0))/(n-num+1)(-:1:n-num+1));}
  if(is_void(ry)) ry = 1.;        if(m==1) ry=[ry(1)]; else {if(numberof(ry)==1) ry=ry(1)(-:1:m);   if(numberof(ry)>m) ry=ry(1:m);     if((num=numberof(ry))<m)   ry = grow(double(ry(1:-1)),double(ry(0))/(m-num+1)(-:1:m-num+1));}

  /* invert the order of the y plot (from the top to the bottom) */
  dy = dy(::-1);
  ry = ry(::-1);
                        
  /* read the default style */
  read_style,style,land,systems,leg,cleg;
  if(!is_void(landscape)) land = landscape;
  
  /* grow the system to be an array */
  systems = array(systems(1),n,m);
  
  /* compute the system dimensions */
  Dvx = (V(2)-V(1) - dx(sum))/n;
  Dvy = (V(4)-V(3) - dy(sum))/m;

  /* normalization of the ratio */
  rx = rx * (V(2)-V(1) - dx(sum)) / rx(sum) / Dvx;
  ry = ry * (V(4)-V(3) - dy(sum)) / ry(sum) / Dvy;
  
  /* foce the system to be square */
  if(square) Dvx = Dvy = min(Dvx,Dvy);

  /* compute the limits of each viewport */
  xmin = (V(1) + grow(0.,Dvx*rx)(psum)(:-1) + grow(0.,dx)(psum))(,-:1:m); if(n==1) xmin = xmin([[1]]);
  xmax = xmin +  Dvx*rx;
  ymin = (V(3) + grow(0.,Dvy*ry)(psum)(:-1) + grow(0.,dy)(psum))(-:1:n,); if(m==1) ymin = ymin([[1]]);
  ymax = ymin + (Dvy*ry)(-,);

  /* if recentering */
  if(square) {
    delta = (V(1)+V(2))/2. - (min(xmin)+max(xmax))/2. ;
    xmin += delta; xmax += delta;
    delta = (V(3)+V(4))/2. - (min(ymin)+max(ymax))/2. ;
    ymin += delta; ymax += delta;
  }

  /* configure the labels */
  if(is_void(fx))           fx=1;
  if(numberof(fx)==1)       fx=fx(1)(-:1:n,-:1:m);
  else if(numberof(fx)==n)  fx=fx(,-:1:m);
  else if(numberof(fx)==m)  fx=fx(-:1:n,);
  if(is_void(fy))           fy=1;
  if(numberof(fy)==1)       fy=fy(1)(-:1:n,-:1:m);
  else if(numberof(fy)==m)  fy=fy(-:1:n,);  /* where f<0, we verifie the delta */
  else if(numberof(fy)==n)  fy=fy(,-:1:m);
  if(m==1) fx=fx(1); else if(is_array((fn=where(fx>0 & grow(dy>0.035,1)(-:1:n,)<0.035)))) fx(fn) = 0.;
  if(n==1) fy=fy(1); else if(is_array((fn=where(fy>0 & grow(1,dx>0.035)(,-:1:m)<0.035)))) fy(fn) = 0.;
  
  /* set the systems */
  systems.viewport = transpose([xmin,xmax,ymin(,::-1),ymax(,::-1)],2);
  systems.ticks.horiz.flags -= 0x020*!fx;
  systems.ticks.vert.flags  -= 0x020*!fy;

  /* set this style in the 'win' window */ 
  window,win;
  set_style,land,systems(*),leg,cleg;

  /* if need, write the style */
  if(file) write_style,file,land,systems,legends,clegends;

  return win;
}

/*---------------------------------------------------------------------------*/

func plTitle(title)
/* DOCUMENT plTitle,title
   Plot a Global title to a sheet.
   SEE ALSO: pltitle, xytitles, pltitle_nm, xytitles_nm
 */
{
  local xt,yt,n,m;
  extern pltitle_height;
  extern pltitle_font;

  read_nm,,n,m,sys;
  xt = (sys(1,1,1)+sys(2,n)) / 2.;
  yt = sys(4,1,1);
  
  plt, title, xt, yt + 0.06,
    font=pltitle_font, justify="CB", height=pltitle_height,tosys=0;  
}


/*---------------------------------------------------------------------------*/

func pltitle_nm(titles)
/* DOCUMENT pltitle_nm, titles
   The same as pltitle, but for n by m window. 'titles' could be array.
   SEE ALSO: plt, xytitles_nm
*/
{
  local n,m,sys,xmin,xmax,tops,imax;
  extern pltitle_height;
  extern pltitle_font;

  read_nm,,n,m,sys;
  xmin = sys(1,)(*); xmax = sys(2,)(*);
  ymin = sys(3,)(*); ymax = sys(4,)(*);
  tpos = [ ((xmax+xmin)/2. )(1:n) , (max(ymax))(-:1:n) ];  
  
  imax = max(numberof(titles),numberof(tpos(,1)));
  if (is_array(titles)) for(i=1; i<=imax; i++) {
    ptitles= numberof(titles)>=1 ? &titles(i%numberof(titles)) : &titles;
    px= numberof(tpos(,1))>=1 ? &tpos(,1)(i%numberof(tpos(,1))) : &tpos(,1);  
    py= numberof(tpos(,2))>=1 ? &tpos(,2)(i%numberof(tpos(,2))) : &tpos(,2);    
    plt, *ptitles, *px, *py + 0.02,
      font=pltitle_font, justify="CB", height=pltitle_height,tosys=0;
  }
}

/*---------------------------------------------------------------------------*/

func xytitles_nm(xtitles, ytitles, adjust)
/* DOCUMENT xytitles_nm, xtitles, ytitles;
            xytitles_nm, xtitles, ytitles, adjust;
   The same as xytitles, but for n by m window. 'xtitles' and/or
   'ytitles' could be array.           
   SEE ALSO: gs_nm, limits_nm, logx_nm, xytitles_nm
 */
{
  local imax,n,m;
  extern pltitle_height;
  extern pltitle_font;
  if (is_void(adjust)) adjust= [0.,0.];

  read_nm,,n,m,sys;
  xmin = sys(1,)(*); xmax = sys(2,)(*);
  ymin = sys(3,)(*); ymax = sys(4,)(*);
  
  xpos = [ ((xmax+xmin)/2. )(1:n) , (min(ymin))(-:1:n) ];
  ypos = [ (min(xmin))(-:1:m) , ((ymax+ymin)/2.)(::n)  ];

  imax = max(numberof(xtitles),numberof(xpos(,1)));
  if (is_array(xtitles)) for(i=1; i<=imax; i++) { 
    pxtitles= numberof(xtitles)>=1 ? &xtitles(i%numberof(xtitles)) : &xtitles;
    px= numberof(xpos(,1))>=1 ? &xpos(,1)(i%numberof(xpos(,1))) : &xpos(,1);  
    py= numberof(xpos(,2))>=1 ? &xpos(,2)(i%numberof(xpos(,2))) : &xpos(,2);  
    plt, *pxtitles, *px, *py + adjust(2) - 0.05,
      font=pltitle_font, justify="CT", height=pltitle_height, tosys=0;
  }
  
  imax = max(numberof(ytitles),numberof(ypos(,1)));
  if (is_array(ytitles)) for(i=1; i<=imax; i++) { 
    pytitles= numberof(ytitles)>=1 ? &ytitles(i%numberof(ytitles)) : &ytitles;
    px= numberof(ypos(,1))>=1 ? &ypos(,1)(i%numberof(ypos(,1))) : &ypos(,1);  
    py= numberof(ypos(,1))>=1 ? &ypos(,2)(i%numberof(ypos(,2))) : &ypos(,2);  
    plt, *pytitles, *px + adjust(1) - 0.05, *py,
      font=pltitle_font, justify="CB", height=pltitle_height, orient=1, tosys=0;
  }
}

/*---------------------------------------------------------------------------*/

func logxy_nm(xflag, yflag)
/* DOCUMENT  logxy_nm, xflag, yflag
   The same as logxy but for n by m window;
   SEE ALSO: gs_nm, limits_nm, logx_nm, xytitles_nm
 */
{
  local nsys,Xmin,Xmax,Ymin,Ymax;  
  get_style,landscape,systems,legends,clegends;
  nsys = numberof(systems);
  for(i=1;i<=nsys;i++) {
    plsys,i;
    logxy, xflag, yflag;
  }
}

/*---------------------------------------------------------------------------*/

func limits_nm(xmin, xmax, ymin, ymax, square=, nice=, restrict=, individual=)
/* DOCUMENT limits_nm, xmin, xmax, ymin, ymax, square=, nice=, restrict=
   The same as limits() but for n by m window;
   SEE ALSO: gs_nm, limits_nm, logx_nm, xylitles_nm
 */
{
  local n,m;
  
  /* if individual, all limits are restore independantly */
  if(individual) _set_limits_nm,;
  /* set the global limits */
  else if(is_void(xmin) && is_void(xmax) && is_void(ymin) && is_void(ymax) && am_subroutine())
    {_set_limits_nm,_get_limits_nm(1);}
  /* set the limits store in xmin array*/
  else if(numberof(xmin)>1) _set_limits_nm,xmin;
  /* else, set the inputs as limits */
  else {
    read_nm,,n,m;
    for(i=1;i<=n*m;i++) {
      plsys,i;limits,xmin,xmax,ymin,ymax,square=square,nice=nice,restrict=restrict;   
    }
  }
  /* return the previous global limits */
  return _get_limits_nm();    
}

func _get_limits_nm(mode)
{
  /* if mode<0 return the current limits, else return the max of this systems */
  /* if abs(mode)=1 return the max over all systems */
  local n,m,l;
  if(is_void(mode)) mode=0;
  read_nm,,n,m,sys;
  l = array(double,5,n,m);
  for(i=1;i<=n;i++)
    for(j=1;j<=m;j++) {
      plsys_nm,i,j,n;
      if(!(mode<0)) {old_l=limits(); limits;}
      l(,i,j) = limits();
      if(!(mode<0)) limits,old_l;
    }
  if(abs(mode)==1) l=[l(1,min,min),l(2,max,max),l(3,min,min),l(4,max,max),l(5,1,1)];
  return l;
}

func _set_limits_nm(l)
{
  local n,m;
  read_nm,,n,m,sys;
  if(is_array(l)) {dim=dimsof(l); if(dim(1)==1) dim=[3,dim(2),1,1];}
  for(i=1;i<=n;i++)
    for(j=1;j<=m;j++) {
      plsys_nm,i,j,n;
      if(is_array(l)) {limits,l(1,i%dim(3),j%dim(4)),l(2,i%dim(3),j%dim(4)),l(3,i%dim(3),j%dim(4)),l(4,i%dim(3),j%dim(4));}
      else limits;
    }
}

/*---------------------------------------------------------------------------*/

func read_nm(win,&n,&m,&sys)
/* DOCUMENT read_nm(win)
            read_nm,win,n,m,sys
   Return the number of horizontal (n) and vertical (m) system in a NxM
   window.
   SEE ALSO: plsys_nm, limits_nm, gs_nm
 */
{
  local ymin, nsys,_land,_sys,_leg,_cleg;
  window,win;
  get_style, _land, _sys, _leg, _cleg;

  nsys = numberof(_sys);
  ymin = _sys.viewport(3,);
  
  /* find n */
  n = where(ymin!=ymin(1));
  if(is_array(n)) n = n(1)-1;
  else n = nsys;

  /* find m */
  m = nsys/n;

  sys    = array(double,4,n,m);
  sys(,*) = _sys.viewport;
  
  return [n,m]
}

/*---------------------------------------------------------------------------*/

func plsys_nm(col,line,n)
/* DOCUMENT plsys_nm,col,line,n;
            plsys_nm,col,line;
   Switch to the colxline system in a NxM (horyzontalxvertical)
   window. n, the number of horyzontal system, is an optional argument
   SEE ALSO: plsys_nm, limits_nm, gs_nm
*/ 
{
  local n;
  if(is_void(n)) n = read_nm()(1);
  plsys,(line-1)*n + col;
}

/*---------------------------------------------------------------------------*/

func range_nm(ymin, ymax, line)
/* DOCUMENT range_nm, ymin, ymax, line
            range_nm, ymin, ymax
            range_nm

   Set the range for a n.m ploting window. 'line' is a vector of the line
   number you want to range. 'ymin' and 'ymax' could be vector of the
   same dimension.
   The 'line' is nil, the default is indgen(m) (all line).
   If 'ymin' and/or 'ymax' are nil, the default is value is the extremum
   of the considered line.
   
   SEE ALSO: logxy_nm, pltitle_nm, rangex_nm ...
 */
{
  local n,m;
  read_nm,,n,m;
  if(is_void(line)) line = indgen(m);
  if(is_void(ymin)) fmin = 1; else fmin = 0;
  if(is_void(ymax)) fmax = 1; else fmax = 0;

  for(l=1;l<=numberof(line);l++) {
    if(fmin) {plsys,(line(l)-1)*n+1;limits;ymin=limits()(3); for(c=2;c<=n;c++) {plsys,(line(l)-1)*n+c; limits;ymin=min(ymin,limits()(3));}}
    if(fmax) {plsys,(line(l)-1)*n+1;limits;ymax=limits()(4); for(c=2;c<=n;c++) {plsys,(line(l)-1)*n+c; limits;ymax=max(ymax,limits()(4));}}
    for(c=1;c<=n;c++) {
      plsys,(line(l)-1)*n + c;
      pymin= numberof(ymin)>=1 ? &ymin(c%numberof(ymin)) : &ymin;  
      pymax= numberof(ymax)>=1 ? &ymax(c%numberof(ymax)) : &ymax;  
      range,*pymin,*pymax;
    }
  }
}

/*---------------------------------------------------------------------------*/

func rangex_nm(xmin, xmax, column)
/* DOCUMENT rangex_nm, xmin, xmax, line
            rangex_nm, xmin, xmax
            rangex_nm

   Set the X-range for a n.m ploting window. 'column' is a vector of the column
   number you want to X-range. 'yxin' and 'xmax' could be vector of the
   same dimension.
   If 'column' is nil, the default is indgen(n) (all column).
   If 'xmin' and/or 'xmax' are nil, the default is value is the extremum
   of the considered column.
   
   SEE ALSO: logxy_nm, pltitle_nm, rangex_nm ...
 */
{
  local n,m;
  read_nm,,n,m;
  if(is_void(column)) column = indgen(n);
  if(is_void(xmin)) fmin = 1; else fmin = 0;
  if(is_void(xmax)) fmax = 1; else fmax = 0;
  
  for(c=1;c<=numberof(column);c++) {
    if(fmin) {plsys,column(c); limits;xmin=limits()(1); for(l=2;l<=m;l++) {plsys,(l-1)*n + column(c); _l=limits();limits;xmin=min(xmin,limits()(1));limits,_l; }} 
    if(fmax) {plsys,column(c); limits;xmax=limits()(2); for(l=2;l<=m;l++) {plsys,(l-1)*n + column(c); _l=limits();limits;xmax=max(xmax,limits()(2));limits,_l; }}
    for(l=1;l<=m;l++) {
      plsys,(l-1)*n + column(c);
      pxmin= numberof(xmin)>=1 ? &xmin(l%numberof(xmin)) : &xmin;  
      pxmax= numberof(xmax)>=1 ? &xmax(l%numberof(xmax)) : &xmax;      
      limits, *pxmin,*pxmax;
    }
  }
}

/*---------------------------------------------------------------------------*/

func test_nm(win)
{
  local nsys,x;
  window,win;
  get_style,landscape,systems,legends,clegends;
  nsys = numberof(systems);

  x = span(-1,1.1,100);

  for(i=1;i<=nsys;i++) {
    plsys,i;
    plg,x^i,x,marks=0,type="solid";
  }
  
  pltitle_nm,"titles";
  xytitles_nm,"xtitles","ytitles";
  plTitle,"GROS TITRE";
  
  return nsys;
}
