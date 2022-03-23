if(_PRINT) write,"#include \"plot.i\"";
local  plot;
/* DOCUMENT plot.i -- Additional routines for plotting in Yorick.
   Provides routines:
 *	  - pla: plot several curves at the same time;
 *	  - plp: plot points/error bars;
 *        - plpa: plot several set of points/error bars at the same times
 *	  - plh: plot in an "histogram" style;
 *	  - pls: plot surface as a filled mesh with contours;
 *	  - pl3s: plot 3D surface;
 *	  - pl3dj: plot disjoint lines in 3D space;
 *	  - pl3t: plot text in 3D space;
 *
 * Copyright (c) 1996-1999, Eric THIEBAUT.
 *	See the file "LICENSE" for information on usage and redistribution
 *	of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 *
 * History:
 *	Jan, 2003 LeBouquin: add the function plpa
 *	Dec, 2003 LeBouquin: add contour plot in color_bar,add keyword edges
 *	Nov, 2003 LeBouquin: add vline option in plh;
 *	Nov, 2003 LeBouquin: add multidimentional option in pla
 *	Nov, 2003 LeBouquin: add dX and dY in pla
 *	$Id: plot.i,v 1.10 2002/05/14 10:03:49 eric Exp $
 *	$Log: plot.i,v $
 *	Revision 1.10  2002/05/14 10:03:49  eric
 *	 - plp: use directly plfp instead of plmk to plot symbols and ticks.
 *	 - plp: draw each symbol inside a unit circle so that they look the
 *	   same size.
 *	 - plp: new "star" symbol (SYMBOL=8) and draw polygon for SYMBOL>=9.
 *	 - plp: new keywords XLO, XHI, YLO, YHI to draw non-symmetrical error
 *	   bars and FILL to draw filled symbols.
 *	 - plp: removed unused keyword HIDE.
 *
 *	Revision 1.9  2001/11/15 09:46:19  eric
 *	 - plp: use plmk to plot symbols and ticks.
 *	 - plp: change usage of keywords SYMBOL and SIZE.
 *	 - plp: remove keyword ASPECT.
 *
 *	Revision 1.8  2001/06/21 12:45:17  eric
 *	- fix indentation and floating point representation.
 *
 *	Revision 1.7  1997/04/03 15:52:39  eric
 *	Added routines color_bar and pl_fc.
 *
 *	Revision 1.6  1996/11/22 11:42:02  eric
 *	 - Add keyword ASPECT for plp.
 *
 *	02/20/96 Eric THIEBAUT: plp;
 *	02/20/96 release 1.1.
 *	02/21/96 Eric THIEBAUT: plh;
 *	02/21/96 release 1.2.
 *	02/21/96 Christophe PICHON and Eric THIEBAUT: pla;
 *	03/04/96 Eric THIEBAUT (from routines written by Christophe PICHON
 *	  and David MUNRO): pl3s, pl3dj pl3t;
 *	03/04/96 release 1.3.
 *	03/04/96 Eric THIEBAUT: added more symbols in plp;
 *	24/03/96 Eric THIEBAUT: added "box" keyword in pl3s and fix some bugs;
 *	24/03/96 release 1.4.
 *	May 3, 1996 by Eric THIEBAUT: added point symbol in routine plp;

   SEE ALSO:
 */
/*---------------------------------------------------------------------------*/

func pl_fc(z, y, x, ireg, levs=, legend=, hide=, type=, width=, color=,
	   colors=, smooth=, marks=, marker=, mspace=, mphase=,
	   triangle=, region=)
{
  d= dimsof(z);
  if (d(1) != 2)
    error, "expecting a 2D array for Z";
  nx= d(2);
  ny= d(3);
  if (is_void(x))
    x= span(1, nx, nx)(,-:1:ny);
  else if (dimsof(x)(1) == 1)
    x= x(,-:1:ny);
  if (is_void(y))
    y= span(1, ny, ny)(-:1:nx,);
  else if (dimsof(y)(1) == 1)
    y= y(-:1:nx,);
  if (is_void(levs)) {
    zmin= min(z);
    zmax= max(z);
    if (zmin >= zmax) levs= zmin;
    else levs= zmin+indgen(9)*0.1*(zmax-zmin);
  }
  plfc, z, y, x, ireg, levs=levs, colors=colors,
    triangle=triangle, region=region;
  plc, z, y, x, ireg, levs=levs, legend=legend, hide=hide, type=type,
    width=width, color=color, smooth=smooth, marks=marks, marker=marker,
    mspace=mspace, mphase=mphase, triangle=triangle, region=region;
}

/*----------------------------------------------------------------------*/

func color_bar(levs, colors, vert=, labs=, adjust=, ecolor=, edges=,
               height=, vport=, format=, font=, flip=)
/* DOCUMENT color_bar
         or color_bar, levs, colors
     Draw a color bar below the current coordinate system.  If LEVS is
     not specified uses plfc_levs (set by previous call to plfc).  If
     COLORS is specified, it should have one more value than LEVS,
     otherwise equally spaced colors are chosen, or plfc_colors if
     plfc_levs was used.  With the vert=1 keyword the color bar appears
     to the left of the current coordinate system (vert=0 is default).
     By default, color_bar will attempt to label some of the color
     interfaces.  With the labs= keyword, you can force the labelling
     algorithm as follows: labs=0 supresses all labels, labs=n forces
     a label at every nth interface, labs=[i,n] forces a label at every
     nth interface starting from interface i (0<=i<=numberof(LEVS)).

     You can specify the viewport coordinates by keyword
     VPORT=[xmin,xmax,ymin,ymax]; by default the colorbar is drawn next to
     the current viewport. You can use the ADJUST keyword to move the bar
     closer to (adjust<0) or further from (adjust>0) the viewport.

     You can specify the string format for labels with keyword FORMAT (default
     "%g"), the font type with keyword FONT and the font
     height with keyword HEIGHT.

   SEE ALSO: plfc
 */
{
  if(is_void(edges)) edges=0;
  if(is_void(flip)) flip=0;
  if (is_void(levs)) {
    if (is_void(plfc_levs)) error, "no levels specified";
    levs= plfc_levs;
    n= numberof(levs)+1;
    if (is_void(colors)) colors= plfc_colors;
  }
    else {
      if(numberof(levs)==2) levs = span(levs(1),levs(2),200);
    n= numberof(levs)+1;
    if (is_void(colors)) colors= bytscl(span(1,n,n),cmin=0.5,cmax=n+0.5);
  }
  if (n != numberof(colors))
    error, "numberof(colors) must be one more than numberof(levs)";

  if (is_void(vport)) vport= viewport();
  if (is_void(adjust)) adjust= 0.0;
  dx= dy= 0.0;
  if (vert) {
    x= (vport(2)+adjust+[0.022,0.042])(-:1:n+1,); x0=x(1,1);x1=x(1,0);
    dx= 0.005;
    y= span(vport(3),vport(4),n+1)(,-:1:2);       y0=y(1,1);y1=y(0,1);
  } else {
    y= (vport(3)-adjust-[0.045,0.065])(-:1:n+1,); y0=y(1,1);y1=y(1,0);
    dy= -0.005;
    x= span(vport(1),vport(2),n+1)(,-:1:2);       x0=x(1,1);x1=x(0,1);
  }
  sys= plsys(0);
  plf,[colors],y,x,edges=edges,ecolor=ecolor, legend="";

  plg,[y0,y0,y1,y1],[x0,x1,x1,x0], closed=1,
    marks=0,color=ecolor,width=1,type=1,legend="";

  plsys, sys;

  if (is_void(labs) || labs(0)>0) {
    if (numberof(levs)>1) {
      dz= levs(dif);
      if (numberof(dz)!=numberof(levs)-1 ||
          anyof((dz>0.0)!=(dz(1)>0.0)) || !dz(1))
        error, "levs must be monotone 1D";
      levs= levs(1:0);
      levs= grow([2*levs(1)-levs(2)],levs,[2*levs(0)-levs(-1)]);
    } else {
      levs= double(levs(1));
      if (!levs) levs= [-1.0,levs,1.0];
      else levs= [0.0,levs,2*levs];
    }
    if (numberof(labs)<2) {
      if (is_void(labs)) labs= (n-1)/4 + 1;
      orig= where(levs<1.0e-9*max(levs(dif)));
      if (numberof(orig)==1) labs= [orig(1)%labs,labs];
      else labs= [(n%labs)/2,labs];
    }

    list= where(indgen(0:n)%labs(2)==labs(1));

    // CP : update in case we don't want to plot the 0 tick
    if (!numberof(list)) {
      list= where(indgen(0:n)%labs(2)==0);
      list = list(2:) ;
    }

    info, x ;

    x= x(list,);
    y= y(list,);
    if (is_void(format)) format= "%3.2g";
    labs= strtrim(swrite(format=format,levs(list)));

    plsys, 0;
    if (flip) {
      pldj, x(,2)-4*dx,y(,1),x(,2)-6*dx,y(,1)-dy, legend=""; // flip either x or y : x for HD163_pola
    } else {
      pldj, x(,2),y(,2),x(,2)+dx,y(,2)+dy, legend="";
    }
    plsys, sys;
    if (flip) {
      //plt1, labs,x(,2)+2*dx,y(,1)-5*dy, justify=(vert?"CH":"CT"),height=height,font=font;
      //plt1, labs,x(,2),y(,1)-5*dy, justify=(vert?"CH":"CT"),height=height,font=font; // for HD163 planet
      plt1, labs,x(,2)-9*dx,y(,1)-5*dy, justify=(vert?"CH":"CT"),height=height,font=font; // for HD163 pola
    } else {
      plt1, labs,x(,2)+3*dx,y(,2)+2*dy, justify=(vert?"CH":"CT"),height=height,font=font;
    }
  }
}

/*----------------------------------------------------------------------*/

func pla(y, x, every=, legend=, hide=, type=, width=, color=, closed=, smooth=,
	 marks=, marker=, mspace=, mphase=, msize=, rays=, arrowl=, arroww=,
         rspace=, rphase=, dx=, dy=, dtype=, tosys=)
{
/* DOCUMENT pla, y, x
         or pla, y
	Plot the  buddle of curves Y versus X labelled by the last indice.
	Y must be 2-dimensional, and X may be 2-dimensional, 1-dimensional
	or omitted.  If X is 2-dimensional, it must have the same dimensions
	as Y and Y(,i) versus X(,i) is plotted for each last indice i.  If
	X is 1-dimensional, it must have the same length as the 1st dimension
	of Y and Y(,i) versus X is plotted for each last indice i.  If X is
	omitted, it defaults to [1, 2, ..., numberof(Y(,1))].

	The plotting keywords of plg are accepted. Most of then could be
        array of the same lenght as the secund dimension of y

        he optional keyword
	every=N which can be used to plot every N curves in the bundle
	(default N=1).

        The dx= and dy= keyword command X and Y errors bars
        (could be 2D array). The type of the bars are command by dtype=

   EXAMPLE
	x =span(0,1,25)(,-:1:25);
	pla, x*transpose(x), marks=0, every=3;
*/
  if(typeof(y) == "complex") y=[y.re,y.im];
  if (is_void(every)) {
    n= 1;
  } else if (numberof(every)!=1 || (n= long(every(1))) <= 0) {
    error, "EVERY must be a scalar >= 1";
  }
  if (dimsof(y)(1) == 1) {y = [y];}
  y = y(,*);
  imax= dimsof(y)(3);
  if (is_void(x)) {
    x2d= 0N;
  } else {
    x2d= dimsof(x)(1) >= 2;
  }

  // default for x
  if(is_void(x)) x=indgen(1:numberof(y(, 1)));
  // loop on i
  for(i= (n+1)/2; i <= imax; i+= n) {
    // test on x
    px= x2d ? &x(,i) : &x;
    // test on scalar keyword
    ptype= numberof(type)(1)>=1 ? &type(i%numberof(type)) : &type;
    pmarker= numberof(marker)(1)>=1 ? &marker(i%numberof(marker)) : &marker;
    pmarks= numberof(marks)(1)>=1 ? &marks(i%numberof(marks)) : &marks;
    pwidth= numberof(width)(1)>=1 ? &width(i%numberof(width)) : &width;
    if (numberof(color) && dimsof(color)(1)==1)
      pcolor= numberof(color)(1)>=1  ? &color(i%numberof(color)) : &color;
    else if (!numberof(color)) pcolor=&[];
    else
      pcolor =numberof(color)(1,)>=1 ? &color(,i%numberof(color)) : &color;

    plegend=numberof(legend)(1)>=1 ? &legend(i%numberof(legend)) : &legend;
    pmsize= numberof(msize)(1)>=1 ? &msize(i%numberof(msize)) : &msize;
    pdtype= numberof(dtype)(1)>=1 ? &dtype(,i%numberof(dtype)) : &dtype;
    ptosys= numberof(tosys) ? tosys(i%numberof(tosys)) : [];
    // test on array keyword
    pdx = (!is_void(dx)  && dimsof(dx)(1)>=2)  ? &dx(,i) : &dx;
    pdy = (!is_void(dy)  && dimsof(dy)(1)>=2)  ? &dy(,i) : &dy;

    plsys, ptosys;
    plg, y(, i), *px, legend=*plegend, hide=hide, type=*ptype, width=*pwidth,
      color=*pcolor, closed=closed, smooth=smooth, marks=*pmarks,
      marker=*pmarker,
      mspace=mspace, mphase=mphase, rays=rays, arrowl=arrowl, arroww=arroww,
      rspace=rspace, rphase=rphase, msize=*pmsize;

    if(!is_void(*pdy))
      pldj, *px, y(, i)-*pdy, *px, y(, i)+*pdy, color=*pcolor, width=*pwidth,
        type=*pdtype;

    if(!is_void(*pdx))
      pldj, *px-*pdx, y(, i), *px+*pdx, y(, i), color=*pcolor, width=*pwidth,
        type=*pdtype;

  }
}

/*----------------------------------------------------------------------*/

func pls_mesh(&x, &xx, d, which=, inhibit=)
/* DOCUMENT err_msg= pls_mesh(x, xx, dimsof(z), which=1/2, inhibit=1/2)

     build X and/or XX arrays of coordinates (abscissa if last argument is
     0/nil; otherwise ordinate) for 2-D array Z.  Normally, the returned
     value is string(0) otherwise it is an error message.

     X is input and output, it will have the same shape as Z and will be
     suitable for contour plots.  XX is purely output, it will have 1 more
     element than Z in each dimension and will be suitable for mesh plots.
     In other words, X(i,j) will furnish the coordinate of the centre of
     cell Z(i,j) whereas XX(i,j), XX(i,j+1), XX(i+1,j) and XX(i+1,j+1)
     will give the coordinates of the corners of cell Z(i,j).

     Assuming the length of Z along the considered dimension is N
     (N must be >= 2) there are 3 possibilities:
       (1) if X is a vector with N elements or has the same shape as Z,
           then X is considered to give the coordinates at the centre of Z
           cells: X is unchanged and output XX is build by interpolating
	   (and extrapolating at the edges) X ;
       (2) if X is a vector with N+1 elements or has 1 more element than Z
           in each dimension, then X is considered to give the coordinates
           at the corners of Z cells: output XX is set to input X and
	   output X is build by interpolating output XX;
       (3) if X is nil, it defaults to [0.5, 1.5, ..., N-0.5] and XX
           defaults to [0, 1, ..., N] along the considered dimension.
     Finally, if X is 1-D, it is expanded in the other direction.

     If keyword WHICH is 1 (the default), abscissa is the dimension of
     interest; otherwise WHICH must be 2 and ordinate is the dimension
     of interest.

     If keyword INHIBIT is 1, then only X output is computed; if INHIBIT
     is 2 then only XX output is computed.

   SEE ALSO: pls, pl3s, plmesh.
 */
{
  xx= [];
  if (is_void(which))
    which= 1;
  do_x= inhibit != 1;
  do_xx= inhibit != 2;
  expand=1;
  if (d(1) != 2 || anyof(d < 2))
    return "Z must be 2-dimensional and have at least 2-by-2 elements";
  n1= d(2);
  n2= d(3);
  n= d(which+1);
  if (is_void((dx= dimsof(x)))) {
    if (do_x)
      x= span(0.5, n-0.5, n);
    if (do_xx)
      xx= span(0, n, n+1);
  } else if (dx(1) == 1) {
    if (dx(2) == n) {
      if (do_xx) {
	xx= x(pcen);
	xx(1)= 2.0 * x(1) - x(2);
	xx(0)= 2.0 * x(0) - x(-1);
      }
    } else if (dx(2) == n+1) {
      xx= x;
      x= do_x ? xx(zcen) : [];
    }
  } else if (dx(1) == 2) {
    expand= 0;
    if (allof(dx == d)) {
      if (do_xx) {
	t= x(pcen,);
	t(1,)= 2.0 * x(1,) - x(2,);
	t(0,)= 2.0 * x(0,) - x(-1,);
	xx= t(,pcen);
	xx(,1)= 2.0 * t(,1) - t(,2);
	xx(,0)= 2.0 * t(,0) - t(,-1);
	t= [];
      }
    } else if (allof(dx == d + [0,1,1])) {
      xx= x;
      x= do_x ? xx(zcen,zcen) : [];
    }
  }
  if (is_void(xx) && is_void(x)) {
    return "X, Y and Z are not compatible";
  }
  if (expand) {
    if (which == 1) {
      if (do_x)
	x= x(,-:1:n2);
      if (do_xx)
	xx= xx(,-:1:n2+1);
    } else {
      if (do_x)
	x= x(-:1:n1,);
      if (do_xx)
	xx= xx(-:1:n1+1,);
    }
  }
  return string(0);
}

func pls(z, y, x, cbar=, viewport=, title=, xtitle=, ytitle=,
	 legend=, hide=, top=, cmin=, cmax=, edges=, ecolor=, ewidth=,
         height=, font=, levs=, nlevs=, type=, width=, color=,
         marks=, marker=, mspace=, mphase=, smooth=)
/* DOCUMENT pls, z, y, x
         or pls, z
     draws surface plot of Z versus (X,Y) as a filled mesh with
     optional contours.  The Z array must be a 2-dimensional array,
     see documentation of pls_mesh for the meaning of X and Y.

     If keyword CBAR is set to non-zero, a color bar is drawn on the
     right of the plot.  The current viewport (in NDC) may be
     specified with keyword VIEWPORT, default is:
       [0.19, 0.60, 0.44, 0.85].

     The appearance of the filled mesh can be modified by means of
     keywords: LEGEND, HIDE, TOP, CMIN, CMAX, EDGES, ECOLOR and EWIDTH
     (see plf documentation).

     Optional contour plot of Z may be superimposed by either keyword
     NLEVS to set the number of contours or by with keyword LEVS to
     specify the level values.  The appearance of the contour plot can
     be modified by means of keywords: LEGEND, HIDE, TYPE, WIDTH,
     COLOR, MARKS, MARKER, MSPACE, MPHASE and SMOOTH (see plc
     documentation).

   SEE ALSO: pls_mesh, pl3s, plc, plf, plmesh.
*/
{
  local r, g, b;	// these variables are used to query colors
  local xx, yy;

  /*
   * Set some defaults.
   */
  if (is_void(edges))
    edges= 0;
  if (is_void(height)) {
    height= 12;
    small = 10;
  } else {
    s= [8,10,12,14,18,24];
    i= where(height == s);
    if (numberof(i) != 1)
      error, "bad font HEIGHT";
    i= i(1);
    small= i > 1 ? s(i-1) : height;
  }
  if (numberof(levs)) {
    nlevs= numberof(levs);
  } else if (is_void(nlevs)) {
      nlevs= 8;
  }

  /*
   * Compute mesh coordinates.
   */
  i= nlevs >= 1 ? 0 : 1;
  if ((msg= pls_mesh(x, xx, dimsof(z), which=1, inhibit=i)) != string(0) ||
      (msg= pls_mesh(y, yy, dimsof(z), which=2, inhibit=i)) != string(0))
    error, msg;

  /*
   * Plot color bar and titles.
   */
  vpmax= [0.127, 0.672, 0.363, 0.908];
  if (numberof(viewport) != 4)
    viewport= [0.19, 0.60, 0.44, 0.85];		// standard viewport
  if (cbar) {
    local r, g, b;
    plsys, 0;
    margin= vpmax(2)-viewport(2);
    x0= viewport(2) + 0.7 * margin;
    x1= viewport(2) + 0.9 * margin;
    y0= viewport(3);
    y1= viewport(4);
    palette, r, g, b, query=1;
    n= numberof(r);
    r= g= b= [];
    pli, char(indgen(n)-1)(-,), legend=string(0), x0, y0, x1, y1;
    plg, [y0,y0,y1,y1,y0], [x0,x1,x1,x0,x0], legend=string(0), marks=0,
      width=1;
    plsys, 1;
  }
  xc= 0.5*(viewport(1)+viewport(2));
  yc= 0.5*(viewport(3)+viewport(4));
  if (!is_void(title)) {
    plt, title, xc, viewport(4) + 0.9 * (vpmax(4) - viewport(4)), tosys=0,
      legend=string(0), justify="CT", path=0,
      font=font, height=height, opaque=opaque;
  }
  if (!is_void(xtitle)) {
    plt, xtitle, xc, vpmax(3) + 0.05 * (viewport(3) - vpmax(3)), tosys=0,
      legend=string(0), justify="CB", path=0,
      font=font, height=small, opaque=opaque;
  }
  if (!is_void(ytitle)) {
    plt, ytitle, vpmax(1) + 0.05 * (viewport(1) - vpmax(1)), yc, tosys=0,
      legend=string(0), justify="LH", path=1,
      font=font, height=small, opaque=opaque;
  }

  /*
   * Plot filled mesh.
   */
  plf, z, yy, xx, legend=legend, hide=hide,
    top=top, cmin=cmin, cmax=cmax, edges=edges, ecolor=ecolor, ewidth=ewidth;
  xx= yy= [];

  /*
   * Plot contours.
   */
  if (nlevs) {
    if (is_void(levs)) {
      zmax= double(max(z));
      zmin= double(min(z));
      levs= zmin + (zmax-zmin) / double(nlevs+1) * indgen(nlevs);
    }
    plc, z, y, x, levs=levs, legend=legend, hide=hide, type=type, width=width,
      color=color, marks=marks, marker=marker, mspace=mspace, mphase=mphase,
      smooth=smooth;
  }
}

/*----------------------------------------------------------------------*/

func plh(y, x, just=, legend=, hide=, type=, width=, color=, marks=, marker=,
	 mspace=, mphase=, vline=, ecolor=, bwidth=)
/* DOCUMENT plh, y, x
         or plh, y
	plots a graph of Y versus X in an "histogram" style (i.e., with
	steps).  Y and X must be 1-D arrays of equal length; if X is
	omitted, it defaults to [1, 2, ..., numberof(Y)].

	The optional keyword JUST set justification of the histogram:
	JUST=1, 2 or 3 makes the graph be left justified, centered or
	right justified respectively along X axis.  Default is centered.
        JUST= 4 --> justified with  bars whose width is given by keyword bwidth.
        Default is 0.5. 1.0 means touching bars.

        The optional keyword VLINE=0 print vertical line betwen the steps.

	Other plotting keywords (legend, hide, type, width, color, marks,
	marker, mspace, and mphase) are passed to the plg routine.

   SEE ALSO: plg, plm, plc, plv, plf, pli, plt, pldj, plfp
             limits, logxy, range, fma, hcp
*/
{
  if(is_void(marks)) marks=0;
  if(is_void(type))  type="solid";
  if(is_void(vline)) vline=0;

  // parse/check arguments
  if (!is_array(y) || dimsof(y)(1)!=1 || (n= numberof(y)) < 2)
    error, "Y must be a vector of at least 2 elements";
  if (is_void(x))
    x= double(indgen(numberof(y)));
  else if (!is_array(x) || dimsof(x)(1)!=1 || numberof(x) != n)
    error, "X must be a vector of same length as Y";
  if (is_void(just))
    just= 2;

  // build new X vector
  n2= 2 * n;
  n4 = 2*n2;
  x2= array(double, n2);
  x4= array(double, n4);
  if (just == 1) {
    // left justify
    x2(1::2)= x;
    x2(2:-1:2)= x(2:);
    x2(0)= 2 * x(0) - x(-1);
  } else if (just == 2) {
    // center
    d= 0.5 * x(dif);
    dx= d(1);
    grow, dx, d, d(0);
    d= [];
    x2(1::2)= x - dx(:-1);
    x2(2::2)= x + dx(2:);
    dx= [];
  } else if (just == 3) {
    // right justify
    x2(1)= 2 * x(1) - x(2);
    x2(2::2)= x;
    x2(3::2)= x(:-1);
  } else if (just == 4) {
    // center
    if (is_void(bwidth)) bwidth=0.5;
    d= 0.5 * bwidth * x(dif);
    dx= d(1);
    grow, dx, d, d(0);
    d= [];
    x2(1::2)= x - dx(:-1);
    x2(2::2)= x
    x4(1::4) = x4(2::4) = x - dx(:-1);
    x4(3::4) = x4(4::4) = x + dx(2:);
    dx= [];
  } else {
    error, "bad value for JUST";
  }

  // build new Y vector
  y2= array(double, n2);
  y2(1::2)= y2(2::2)= y;

  y4= array(double, n4);
  y4(2::4) = y4(3::4) = y

  //Filled graph :
  if (!is_void(ecolor)) {
    l = limits();
    if (structof(ecolor) == string) {
      n = where(ecolor == __pl_color_list);
      if (numberof(n)!=1) error, "unrecognized color name: "+ecolor;
      ecolor = char(-n(1));
    } else if (structof(ecolor) != char) {
      ecolor = char(ecolor);
    }
    if (just==4) {
      plfp, array(ecolor,1) , _(y4,l(3),l(3)) , _(x4,x4(0),x4(1)), numberof(y4)+2;
    } else {
      plfp, array(ecolor,1) , _(y2,l(3),l(3)) , _(x2,x2(0),x2(1)), numberof(y2)+2;
    }
  }


  // plot the graph
  if (just==4) {
    plg, y4, x4,
      legend=legend, hide=hide, type=type, width=width, color=color,
      marks=marks, marker=marker, mspace=mspace, mphase=mphase;
  } else {
    plg, y2, x2,
      legend=legend, hide=hide, type=type, width=width, color=color,
      marks=marks, marker=marker, mspace=mspace, mphase=mphase;
  }

    //================================
  if(vline==1) {
    l = limits();
    pldj,x2,y2,x2,y2*0.+l(3),type=type, width=width, color=color;
  }

  if (0) {
    dens = 3;
    nb_point = 40;
    dimx = numberof(x2-1)*dens;

    test = array(int, nb_point , dimx);


    ty  = span(min(y2), max(y2), nb_point);
    tx  = y2(-:1:dens,)(*);
    test = ty(-,) < tx;  ou = where(test);
    xp  = span(min(x2), max(x2), dimx)(,-:1:nb_point)(*)(ou);
    yp  = ty(-:1:dimx,)(*)(ou);


    plp, yp , xp , symbol=15, size=.3;
  }

}

/*---------------------------------------------------------------------------*/

func plp(y, x, dx=, xlo=, xhi=, dy=, ylo=, yhi=, size=, symbol=, ticks=,
	 legend=, type=, width=, color=, fill=)
/* DOCUMENT plp, y, x
       -or- plp, y, x, dx=sigma_x, dy=sigma_y
     Plots points (X,Y) with symbols and/or  error bars.  X, and Y may have
     any dimensionality, but  must have the same number  of elements.  If X
     is nil, it defaults to indgen(numberof(Y)).

     Keyword SYMBOL may be used to choose the shape of each symbol:
       0    nothing (just draw error bars if any)
       1    square
       2    cross (+ sign)
       3    triangle
       4    circle (hexagon)
       5    diamond
       6    cross (rotated 45 degrees)   <- this is the default
       7    triangle (upside down)
       8    star
     >=9    SYMBOL-side polygon

     Keyword  SIZE may be used to  change the size of the  symbols and tick
     marks (SIZE acts as a multiplier, default value is 1.0).

     If value of  keyword FILL is true (non-nil  and non-zero), symbols are
     filled with COLOR (default is to draw open symbols).

     Keywords XLO, XHI, YLO, and/or YHI  can be used to indicate the bounds
     of the optional  error bars (default is to draw  no error bars).  Only
     specified bounds get plotted as  error bars. If value of keyword TICKS
     is true  (non-nil and non-zero), ticks  get drawn at  the endpoints of
     the error bars.   Alternatively, keywords DX and/or DY  can be used to
     plot  error bars  as segments  from XLO=X-DX  to XHI=X+DX  and/or from
     YLO=Y-DY to  YHI=Y+DY.  If keyword  DX (respectively DY) is  used, any
     value of XLO and XHI (respectively YLO and YHI) is ignored.

     The other keywords are the same as for pldj (TYPE is only used to draw
     error bars):
   KEYWORDS: legend, type, width, color.

   SEE ALSO: pldj, plg, plm, plc, plv, plf, pli, plt, pldj, plfp, plmk,
             limits, logxy, range, fma, hcp. */
{
  /* NDC units for symbols/ticks (one pixel = 0.00125268 NDC at 75 DPI) */
  u0 = 0.0;       // zero
  u1 = 0.00500214; // radius of about 5 pixels at 75 DPI
  if (! is_void(size)) u1 *= size;

  /* parse color */
  if (is_void(color)) color = char(-2); /* fg */
  if (structof(color) == string) {
    n = where(color == __pl_color_list);
    if (numberof(n)!=1) error, "unrecognized color name: "+color;
    color = char(-n(1));
  } else if (structof(color) != char) {
    color = char(color);
  }

  /* default X */
  if (is_void(x)) (x = array(double, dimsof(y)))(*) = indgen(numberof(y));

  /* error bars */
  if (is_void(dx)) {
    err = (! is_void(xlo)) + 2*(! is_void(xhi));
  } else {
    xlo = x - dx;
    xhi = x + dx;
    err = 3;
  }
  if (err) {
    pldj, (is_void(xlo) ? x : xlo), y, (is_void(xhi) ? x : xhi), y,
      type=type, width=width, color=color;
    if (ticks) {
      xm = [ u0, u0];
      ym = [-u1, u1];
      if      (err == 1) __plp,   y,       xlo;
      else if (err == 2) __plp,   y,       xhi;
      else               __plp, [y, y], [xlo, xhi];
    }
    xhi = xlo = [];
  }
  if (is_void(dy)) {
    err = (! is_void(ylo)) + 2*(! is_void(yhi));
  } else {
    ylo = y - dy;
    yhi = y + dy;
    err = 3;
  }
  if (err) {
    pldj, x, (is_void(ylo) ? y : ylo), x, (is_void(yhi) ? y : yhi),
      type=type, width=width, color=color;
    if (ticks) {
      xm = [-u1, u1];
      ym = [ u0, u0];
      if      (err == 1) __plp,    ylo,       x;
      else if (err == 2) __plp,    yhi,       x;
      else               __plp, [ylo, yhi], [x, x];
    }
    yhi = ylo = [];
  }

  /* symbols */
  if (! symbol) {
    if (is_void(symbol)) symbol = 6;
    else return;
  }
  if (symbol == 1) {
    /* square */
    u2 = u1*sqrt(0.5);
    xm = [-u2, u2, u2,-u2];
    ym = [ u2, u2,-u2,-u2];
  } else if (symbol == 2) {
    /* + cross */
    xm = [-u1, u1, u0, u0, u0, u0];
    ym = [ u0, u0, u0, u1,-u1, u0];
    fill = 0;
  } else if (symbol == 3) {
    /* triangle */
    u2 = u1*0.5;
    u3 = u1*sqrt(0.75);
    xm = [u0, u3,-u3];
    ym = [u1,-u2,-u2];
  } else if (symbol == 4) {
    /* hexagon */
    u2 = u1*0.5;
    u3 = u1*sqrt(0.75);
    xm = [ u1, u2,-u2,-u1,-u2, u2];
    ym = [ u0, u3, u3, u0,-u3,-u3];
  } else if (symbol == 5) {
    /* diamond */
    xm = [u1, u0,-u1, u0];
    ym = [u0, u1, u0,-u1];
  } else if (symbol == 6) {
    /* x cross (rotated 45 degrees) */
    u2 = u1*sqrt(0.5);
    xm = [u2,-u2, u0, u2,-u2, u0];
    ym = [u2,-u2, u0,-u2, u2, u0];
    fill = 0;
  } else if (symbol == 7) {
    /* triangle (upside down) */
    u2 = u1*0.5;
    u3 = u1*sqrt(0.75);
    xm = [ u0, u3,-u3];
    ym = [-u1, u2, u2];
  } else if (symbol == 8) {
    /* 5 branch star */
    /* Notations: C18 = cos(18*ONE_DEGREE)
     *            S18 = sin(18*ONE_DEGREE)
     *            C54 = cos(54*ONE_DEGREE)
     *            S54 = sin(54*ONE_DEGREE)
     */
    u2 = 0.224514*u1; // C54*S18/S54
    u3 = 0.309017*u1; // S18
    u4 = 0.951057*u1; // C18
    u5 = 0.363271*u1; // C18*S18/S54
    u6 = 0.118034*u1; // S18*S18/S54
    u7 = 0.587785*u1; // C54
    u8 = 0.809017*u1; // S54
    u9 = 0.381966*u1; // S18/S54
    xm = [ u0, u2, u4, u5, u7, u0,-u7,-u5,-u4,-u2];
    ym = [ u1, u3, u3,-u6,-u8,-u9,-u8,-u6, u3, u3];
  } else if (symbol > 0) {
    /* N-side polygon in unit circle */
    PI = 3.141592653589793238462643383279503;
    a = (2.0*PI/symbol)*indgen(0:symbol-1);
    xm = u1*cos(a);
    ym = u1*sin(a);
  } else {
    error, "bad SYMBOL value";
  }
  __plp, y, x;
  //fx=abs(limits()(1)-limits()(2))*0.05;
  //fy=abs(limits()(3)-limits()(4))*0.05;
  //limits,limits()(1)-fx,limits()(2)+fx,limits()(3)-fy,limits()(4)+fy;
}

func __plp(y, x)
/* DOCUMENT __plp, x, y;
     Private routine used by plp. */
{
  extern xm, ym, color, fill, legend, width;
  local z;
  n = array(1, 1 + numberof(y));
  n(1) = numberof(ym);
  if (fill && n(1) > 2) {
    if (numberof(color) == 3) {
      z = array(char, 3, numberof(n));
      z(,) = color;
    } else {
      z = array(color, numberof(n));
    }
  }
  plfp, z, grow(ym,y(*)), grow(xm,x(*)), n,
    legend=legend, edges=1, ewidth=width, ecolor=color;
}

/*----------------------------------------------------------------------*/

func plpa(y, x, dx=, xlo=, xhi=, dy=, ylo=, yhi=, size=, symbol=, ticks=,
	 legend=, type=, width=, color=, fill=, every=, tosys=)
{
/* DOCUMENT plpa, y, x
         or plpa, y
	Plot with marker the buddle of curves Y versus X labelled by the
        last indice.Y must be 2-dimensional, and X may be 2-dimensional,
        1-dimensional
	or omitted.  If X is 2-dimensional, it must have the same dimensions
	as Y and Y(,i) versus X(,i) is plotted for each last indice i.  If
	X is 1-dimensional, it must have the same length as the 1st dimension
	of Y and Y(,i) versus X is plotted for each last indice i.  If X is
	omitted, it defaults to [1, 2, ..., numberof(Y(,1))].

	The plotting keywords of plp are accepted. Most of then could be
        array of the same lenght as the secund dimension of y

        The optional keyword
	every=N which can be used to plot every N curves in the bundle
	(default N=1).

   EXAMPLE
        n = 10;
        x = span(0,1,25);
	y = x(,-)^indgen(1:n)(-,);
	plpa, y, x,
        every=1,
        symbol=indgen(1:n),
        color=-5-indgen(1:n),
        legend=swrite(indgen(1:n));
*/
  if(typeof(y) == "complex") y=[y.re,y.im];
  if (is_void(every)) {
    n= 1;
  } else if (numberof(every)!=1 || (n= long(every(1))) <= 0) {
    error, "EVERY must be a scalar >= 1";
  }
  if (dimsof(y)(1) == 1) {y = [y];}
  y = y(,*);
  imax= dimsof(y)(3);
  if (is_void(x)) {
    x2d= 0N;
  } else {
    x2d= dimsof(x)(1) >= 2;
  }
  // default for x
  if(is_void(x)) x=indgen(1:numberof(y(, 1)));
  // loop on i
  for(i= (n+1)/2; i <= imax; i+= n) {
    // test on x
    px= x2d ? &x(,i) : &x;
    // test on scalar keyword
    ptype= numberof(type)(1)>=1 ? &type(i%numberof(type)) : &type;
    pwidth= numberof(width)(1)>=1 ? &width(i%numberof(width)) : &width;
     if (numberof(color) && dimsof(color)(1)==1)
      pcolor= numberof(color)(1)>=1  ? &color(i%numberof(color)) : &color;
     else if (!numberof(color)) pcolor=&[];
     else
       pcolor =numberof(color)(1,)>=1 ? &color(,i%numberof(color)) : &color;
    plegend=numberof(legend)(1)>=1 ? &legend(i%numberof(legend)) : &legend;
    psize= numberof(size)(1)>=1 ? &size(i%numberof(size)) : &size;
    psymbol= numberof(symbol)(1)>=1 ? &symbol(i%numberof(symbol)) : &symbol;
    pticks= numberof(ticks)(1)>=1 ? &ticks(i%numberof(ticks)) : &ticks;
    pfill= numberof(fill)(1)>=1 ? &fill(i%numberof(fill)) : &fill;
    ptosys= numberof(tosys) ? tosys(i%numberof(tosys)) : [] ;
    // test on array keyword
    pdx = (!is_void(dx)  && dimsof(dx)(1)>=2)  ? &dx(,i) : &dx;
    pdy = (!is_void(dy)  && dimsof(dy)(1)>=2)  ? &dy(,i) : &dy;
    pxlo= (!is_void(xlo) && dimsof(xlo)(1)>=2) ? &dxlo(,i) : &xlo;
    pxhi= (!is_void(xhi) && dimsof(xhi)(1)>=2) ? &dxhi(,i) : &xhi;
    pylo= (!is_void(ylo) && dimsof(ylo)(1)>=2) ? &dxlo(,i) : &ylo;
    pyhi= (!is_void(yhi) && dimsof(yhi)(1)>=2) ? &dxhi(,i) : &yhi;

    plsys, ptosys;
    plp,y(, i), *px, dx=*pdx, xlo=*pxlo, xhi=*pxhi, dy=*pdy, ylo=*pylo,
      yhi=*pyhi, size=*psize, symbol=*psymbol, ticks=*pticks,
      legend=*plegend, type=*ptype, width=*pwidth, color=*pcolor,
      fill=*pfill;
  }
}

func plta(y, x, text, legend=, hide=, color=, font=, height=, opaque=,
          orient=, justify=, tosys=)

{
/* DOCUMENT plta, y, x, text
         or plta, y, text
*/
  if (is_void(tosys)) tosys=1;
  if(typeof(y) == "complex") y=[y.re,y.im];
  if (is_void(every)) {
    n= 1;
  } else if (numberof(every)!=1 || (n= long(every(1))) <= 0) {
    error, "EVERY must be a scalar >= 1";
  }
  if (is_void(y) && is_void(x) && !is_void(text) && tosys(1)==0) {
    y=span(0.9,0.1, numberof(text(1,))) ;
    x=span(0.1,0.7, numberof(text(,1))) ;
    y=y(-:1:numberof(x),);  x=x(,-:1:numberof(y));
  }
  if (is_void(y) && !is_void(text))
    y=indgen(numberof(text(1,)))(-:1:numberof(text(,1)),);
  else if (numberof(y)==1) y = y(-:1:numberof(x));
  if (numberof(x)==1) x = x(-:1:numberof(y));

  if (dimsof(y)(1) == 1) {y = [y];}
  y = y(,*);
  imax= dimsof(y)(3);
  if (is_void(x)) {
    x2d= 0N;
    x=indgen(1:numberof(y(, 1)));
  } else {
    x2d= dimsof(x)(1) >= 2;
  }

  if (is_void(text)) {
    text2d= 0N;
    text=indgen(1:numberof(y(, 1)));
  } else {
    text2d= dimsof(text)(1) >= 2;
  }

  // default for x



  // loop on i
  for(i= (n+1)/2; i <= imax; i+= n) {
    // test on x
    px= x2d ? &x(,i) : &x;
    ptext = text2d ? &text(,i) : &text;
    // test on scalar keyword
    pjustify= numberof(justify)(1)>=1 ? &justify(i%numberof(justify)) : &justify;
    pheight= numberof(height)(1)>=1 ? &height(i%numberof(height)) : &height;
    pcolor= numberof(color)(1)>=1 ? &color(i%numberof(color)) : &color;
    plegend=numberof(legend)(1)>=1 ? &legend(i%numberof(legend)) : &legend;
    phide= numberof(hide)(1)>=1 ? &hide(i%numberof(hide)) : &hide;
    popaque= numberof(opaque)(1)>=1 ? &opaque(i%numberof(opaque)) : &opaque;
    porient= numberof(orient)(1)>=1 ? &orient(i%numberof(orient)) : &orient;
    pfont= numberof(font)(1)>=1 ? &font(i%numberof(font)) : &font;
    ptosys= numberof(tosys)(1)>=1 ? &tosys(i%numberof(tosys)) : &tosys;
    // test on array keyword

    for (j=1; j<=numberof(y(,i)); j++)
      plt, (typeof((*ptext)(j)) == "string" ? (*ptext)(j) : pr1((*ptext)(j))) ,
        (*px)(j), y(j , i),
        hide= *phide, opaque=*popaque, orient=*porient,
        legend=*plegend, justify=*pjustify, height=*pheight,
        color=*pcolor,
        font=*pfont, tosys=*ptosys;
  }
}
//================================================================

local __pl_color_list;
__pl_color_list = ["bg","fg","black","white","red","green","blue",
                  "cyan","magenta","yellow"];
/* DOCUMENT __pl_color_list - private list of color names
 *
 *
*/


/*----------------------------------------------------------------------*/

/*
 * 3D TRANSFORM:
 * -------------
 * Let (X,Y,Z)  be  the  data  coordinates  (in  a  direct  frame) and
 * (XP,YP,ZP) be the coordinates in the view frame (also direct, XP is
 * from left to right, YP is from bottom to top and  ZP  points toward
 * the observer), then  for  an  altitude  ALT  (angle  of  view above
 * XY-plane) and an azimuth  AZ  (angle  of  view  around  Z-axis) the
 * coordinates transform is obtained  by  a  rotation  around  Oz with
 * angle AZ (azimuth), followed by a rotation around Ox with  angle AX
 * (AX = ALT - 90deg):
 *   XP = X cos(AZ) - Y sin(AZ)
 *   YP = U cos(AX) - Z sin(AX) =   U sin(ALT) + Z cos(ALT)
 *   ZP = U sin(AX) + Z cos(AX) = - U cos(ALT) + Z sin(ALT)
 * where:
 *   U  = X sin(AZ) + Y cos(AZ)
 */

func _pl3xyz(&xp, &yp, &zp, x, y, z)
{
/* DOCUMENT _pl3xyz, xp, yp, zp, x, y, z
     transform data coordinates (X,Y,Z) into viewer coordinates
     (XP,YP,ZP) for an externally defined altitude ALT and
     azimuth AZ (in degrees).
 */
  extern alt, az;
  if (is_void(az))  az= 30.0;	// angle of view around z-axis
  if (is_void(alt)) alt= 45.0;	// angle of view above xy-plane
  d2r= pi / 180.0;
  xp= x * (c= cos(az * d2r)) - y * (s= sin(az * d2r));
  zp= x * s + y * c;
  yp= z * (c= cos(alt * d2r)) + zp * (s= sin(alt * d2r));
  zp= z * s - zp * c;
}

func _pl3xy(&xp, &yp, x, y, z)
{
/* DOCUMENT _pl3xy, xp, yp, x, y, z
     transform data coordinates (X,Y,Z) into viewer coordinates
     (XP,YP) for an externally defined altitude ALT and
     azimuth AZ (in degrees).
 */
  extern alt, az;
  if (is_void(az))  az= 30.0;	// angle of view around z-axis
  if (is_void(alt)) alt= 45.0;	// angle of view above xy-plane
  d2r= pi / 180.0;
  xp= x * (c= cos(az * d2r)) - y * (s= sin(az * d2r));
  yp= z * cos(alt * d2r) + (x * s + y * c) * sin(alt * d2r);
}

/*----------------------------------------------------------------------*/

func pl3t(text, x, y, z, alt=, az=, legend=, hide=, color=, font=, height=,
	  opaque=, path=, justify=, tosys=)
{
/* DOCUMENT pl3t, text, x, y, z, alt=alt, az=az, tosys=0/1
     plots TEXT (a string) at the point (X,Y,Z) in a 3-dimensional
     graph view from altitude ALT (default 45) and azimuth AZ
     (default 30) both in degrees.   TEXT, X, Y and Z may be arrays
     with the same number of elements.

     Other optional keywords are:
       legend, hide, color, font, height, opaque, path, justify and tosys
     and have the same meaning as in plt.

   SEE ALSO: plt.
 */
  local xp, yp;
  _pl3xy, xp, yp, x, y, z;
  n= numberof(text);
  if (n == 1) {
    plt, text, xp, yp,
      legend=legend, hide=hide, color=color,font=font, height=height,
      opaque=opaque, path=path, justify=justify, tosys=tosys;
  } else {
    for (i=1; i<=n; i++) {
      plt, text(i), xp(i), yp(i),
	legend=legend, hide=hide, color=color,font=font, height=height,
	opaque=opaque, path=path, justify=justify, tosys=tosys;
    }
  }
}

/*----------------------------------------------------------------------*/

func pl3dj(x0, y0, z0, x1, y1, z1, alt=, az=,
	   legend=, hide=, type=, width=, color=)
{
/* DOCUMENT pl3dj, x0, y0, z0, x1, y1, z1, alt=alt, az=az
     plots disjoint lines from (X0,Y0,Z0) to (X1,Y1,Z1) in a 3-dimensional
     graph view from altitude ALT (default 45) and azimuth AZ (default 30)
     both in degrees.  X0, Y0, Z0, X1, Y1 and Z1 must have the same shapes.

     Additional keywords are those accepted by pldj: legend, hide, type,
     width, and color.

   SEE ALSO: pldj, pl3t, pl3s.
 */
  local x0p, y0p, x1p, y1p;
  _pl3xy, x0p, y0p, x0, y0, z0;
  _pl3xy, x1p, y1p, x1, y1, z1;
  pldj, x0p, y0p, x1p, y1p,
    legend=legend, hide=hide, type=type, width=width, color=color;
}

/*----------------------------------------------------------------------*/

func pl3s(z, y, x, alt=, az=, axis=, box=, acolor=, fill=, legend=, hide=,
          edges=, ecolor=, ewidth=, height=, font=)
/* DOCUMENT pl3s, z, y, x, fill=0/1/2
         or pl3s, z, fill=0/1/2
      draws 3-D surface plot of Z versus (X,Y).   The  Z  array  must  be a
      2-dimensional array, say NX-by-NY and X  and  Y  must  have  the same
      shape as Z or be vectors  of  length  NX  and  NY  respectively.   If
      omitted, X and Y are set to the first and second  indice  value  of Z
      respectively.

      The FILL keyword indicates the kind of plot: 0 (default) for  3D wire
      frames, 1 for 3D mesh filled with intensity,  2  for  3D  mesh shaded
      with light source aligned with observer.

      The altitude and azimuth angles (in degrees) can be set with keywords
      ALT and AZ, their default values are 30 and 45 deg.

      A solid edge can optionally be drawn around each zone by  setting the
      EDGES keyword non-zero.  ECOLOR and EWIDTH determine  the  edge color
      and width.

      Frame axis can optionally be drawn around  the  plot  by  setting the
      AXIS keyword non-zero.  The  color  of  the  axis  and  label  can be
      modified with keyword ACOLOR.

      If BOX keyword non-zero, the 3-D box borders are drawn (with the same
      color as the axis, i.e., ACOLOR).

   EXAMPLE
     It is usually better to select an eventually new window and choose
     the "nobox" style:
       window, max(0, current_window()), wait=1, style="nobox.gs";
       x= span(-3,3,50);
       y= span(-2,2,40);
       z= cos(x(,-) * y(-,));
       pl3s, z, y, x, axis=1, fill=2, edges=1, font="timesBI", height=10;

      The following keywords are legal (each has a separate help entry):
   KEYWORDS: legend, hide, region, edges, ecolor, ewidth, font, height.
   SEE ALSO: pl3dj, pl3t, plg, plm, plc, plv, plf, pli, plt, pldj, plfp,
     plmesh, limits, range, fma, hcp, palette, bytscl, ...
*/
{
  /*
   * Check dimensions of input arrays.
   */
  local tmp;
  if ((msg= pls_mesh(x, tmp, dimsof(z), which=1, inhibit=2)) != string(0) ||
      (msg= pls_mesh(y, tmp, dimsof(z), which=2, inhibit=2)) != string(0))
    error, msg;
  tmp= [];

  /*
   * Rescale arrays.
   */
  xspan= (xmax= max(x)) - (xmin= min(x));
  yspan= (ymax= max(y)) - (ymin= min(y));
  zspan= (zmax= max(z)) - (zmin= min(z));
  if (xspan <= 0.0 || yspan <= 0.0)
    error, "X and/or Y are constant";
  if (zspan <= 0.0) {
    zmin= -1.0;
    zmax= +1.0;
    zspan= 2.0;
    z(*)= 0.0;
  } else {
    z= (z - zmin) / zspan;
  }
  x= (x - xmin) / xspan;
  y= (y - ymin) / yspan;

  /*
   * Insure that angles are in the range [-180,180].
   */
  if (is_void(alt))                alt  =  45.0;
  else if ((alt %= 360.0) > 180.0) alt -= 360.0;
  else if (alt < -180.0)           alt += 360.0;
  if (is_void(az))                 az   =  30.0;
  else if ((az %= 360.0) > 180.0)  az  -= 360.0;
  else if (az < -180.0)            az  += 360.0;

  /*
   * Plot an invisible box around the plot to left some space
   * around for the axis labels.
   */
  if (axis) {
    local bxp, byp;
    if (is_void(height))
      height= 10.0;
    q= height / 75.0;
    _pl3xy, bxp, byp,
      [0, 1, 1, 0, 0, 1, 1, 0] + q * [-1, 1, 1,-1,-1, 1, 1,-1],
      [0, 0, 1, 1, 0, 0, 1, 1] + q * [-1,-1, 1, 1,-1,-1, 1, 1],
      [0, 0, 0, 0, 1, 1, 1, 1] + q * [-1,-1,-1,-1, 1, 1, 1, 1];
    bxmin= min(bxp); bxmax= max(bxp);
    bymin= min(byp); bymax= max(byp);
    pldj, bxmin, bymin, bxmax, bymax, type="none";
    pldj, bxmin, bymax, bxmax, bymin, type="none";
  }

  /*
   * Plot the rear of the 3-D box: must figure out which box faces
   * are seen, and then which box edges are seen.  The plotting order
   * is (1) rear part of the box, (2) surface mesh and (3) front
   * part of the box and axis.
   */
  if (box) {
    local bxp0, byp0, bxp1, byp1;	// end-points of box edges
    local nxp, nyp, nzp;		// vectors normal to box faces
    _pl3xy, bxp0, byp0,
      [0,1,1,0,0,1,1,0,0,1,1,0],
      [0,0,1,1,0,0,1,1,0,0,1,1],
      [0,0,0,0,1,1,1,1,1,1,1,1];
    _pl3xy, bxp1, byp1,
      [1,1,0,0,1,1,0,0,0,1,1,0],
      [0,1,1,0,0,1,1,0,0,0,1,1],
      [0,0,0,0,1,1,1,1,0,0,0,0];
    _pl3xyz, nxp, nyp, nzp,
      [ 0, 0, 0, 1, 0,-1],
      [ 0, 0,-1, 0, 1, 0],
      [-1, 1, 0, 0, 0, 0];
    face_edges=[[1,2,3,4], [5,6,7,8], [1,5,9,10], [2,6,10,11],
		[3,7,11,12], [4,8,9,12]];
    visible = array(0, 12);
    visible(face_edges(, where(nzp >= 0.0))(*)) = 1;
    fore= where(visible);
    back= where(!visible);
    pldj, bxp0(back), byp0(back), bxp1(back), byp1(back), type=1, color=acolor;
  }

  /*
   * Rotate the surface so as to have the drawing starting at back
   * end (i.e., hidden surfaces are drawn first).  To this end, the
   * first thing to do is to localize the corner of surface z=0
   * which is the farest from the observer.  Depending on the
   * position of the farest corner, X, Y and Z arrays
   * may have to be scrambled.  After what, the 3-D projection
   * of the surface can be computed.
   */
  local xp, yp, zp;
  _pl3xyz, xp, yp, zp,
    [x(0,1), x(1,0), x(0,0), x(1,1)],
    [y(0,1), y(1,0), y(0,0), y(1,1)],
    0;
  far = zp(mnx)(1);		// index of farest point
  if (far == 1) {
    x= x(::-1,);
    y= y(::-1,);
    z= z(::-1,);
  } else if (far == 2) {
    x= x(,::-1);
    y= y(,::-1);
    z= z(,::-1);
  } else if (far == 3) {
    x= x(::-1,::-1);
    y= y(::-1,::-1);
    z= z(::-1,::-1);
  }
  _pl3xyz, xp, yp, zp, x, y, z;

  /*
   * Plot surface.
   */
  if (!fill) {
    colors= [];
    edges= 1;
  } else if (fill == 1) {
    colors= bytscl(z, cmin=0.0, cmax=1.0);
  } else {
    /* compute the two median vectors for each cell */
    m0x= xp(dif,zcen);
    m0y= yp(dif,zcen);
    m0z= zp(dif,zcen);
    m1x= xp(zcen,dif);
    m1y= yp(zcen,dif);
    m1z= zp(zcen,dif);
    /* define the normal vector to be their cross product */
    nx= m0y*m1z - m0z*m1y;
    ny= m0z*m1x - m0x*m1z;
    nz= m0y*m1x - m0x*m1y;
    m0x= m0y= m0z= m1x= m1y= m1z= [];
    colors= bytscl(nz);
    //colors= bytscl(nz / abs(nx, ny, nz), cmin=0.0, cmax=1.0);
    nx= ny= nz= [];
  }
  plf, colors, yp, xp, legend=legend, hide=hide,
    edges=edges, ecolor=ecolor, ewidth=ewidth;
  xp= yp= zp= colors= [];

  /*
   * Plot the axis and the front of the box.
   */
  if (axis) {
    if (far == 1) {
      px= [ 0, 1, 0]; vx= [ 0, 1, 0];
      py= [ 0, 0, 0]; vy= [-1, 0, 0];
      pz= [ 1, 1, 0]; vz= [ 1, 0, 0];
    } else if (far == 2) {
      px= [ 0, 0, 0]; vx= [ 0,-1, 0];
      py= [ 1, 0, 0]; vy= [ 1, 0, 0];
      pz= [ 0, 0, 0]; vz= [-1, 0, 0];
    } else if (far == 3) {
      px= [ 0, 0, 0]; vx= [ 0,-1, 0];
      py= [ 0, 0, 0]; vy= [-1, 0, 0];
      pz= [ 0, 1, 0]; vz= [ 0, 1, 0];
    } else {
      px= [ 0, 1, 0]; vx= [ 0, 1, 0];
      py= [ 1, 0, 0]; vy= [ 1, 0, 0];
      pz= [ 1, 0, 0]; vz= [ 0,-1, 0];
    }
    _pl3tick, xmin, xmax, px, [1,0,0], 0.02 * vx,
      height=height, font=font, color=acolor;
    _pl3tick, ymin, ymax, py, [0,1,0], 0.02 * vy,
      height=height, font=font, color=acolor;
    _pl3tick, zmin, zmax, pz, [0,0,1], 0.02 * vz,
      height=height, font=font, color=acolor;
  }
  if (box) {
    pldj, bxp0(fore), byp0(fore), bxp1(fore), byp1(fore), type=1, color=acolor;
  }
}

/*----------------------------------------------------------------------*/

func _pl3tick(tmin, tmax, p, d, v, color=, font=, height=, opaque=, path=)
/* DOCUMENT _pl3tick, p, d, v
     draw axis between (P(1),P(2),P(3)) and ((P+D)(1),(P+D)(2),(P+D)(3))
     and ticks with vectorial length (V(1),V(2),V(3)).
 */
{
  extern alt, az;
  local px, py, dx, dy, vx, vy;

  /*
   * Compute projections, draw axis and figure out which
   * justification is the best one..
   */
  _pl3xy, px, py, p(1), p(2), p(3);
  _pl3xy, dx, dy, d(1), d(2), d(3);
  _pl3xy, vx, vy, v(1), v(2), v(3);
  pldj, px, py, px+dx, py+dy, color=color;
  if (where(d != 0)(1) == 3) {
    // horizontal ticks for z-axis
    vx= vx <= 0.0 ? -max(abs(v)) : max(abs(v));
    vy= 0.0;
  }
  justify= vx > 0.0 ? "L" : (vx ? "R" : "C");
  justify += vy > 0.0 ? "B" : (vy ? "T" : "H");

  /*
   * Compute step size in axis direction to get approximatively 6 ticks
   * and plot ticks.
   */
  tspan= tmax - tmin;
  tstep= 10.0^floor(log10(0.2 * tspan) + 0.5);
  if (tspan / tstep < 4)
    tstep *= 0.5;
  else if (tspan / tstep > 8)
    tstep *= 2.0;
  imin= ceil(tmin / tstep);
  imax= floor(tmax / tstep);
  ni= long(imax - imin + 1.0);
  t= tstep * span(imin, imax, ni);
  tn= (t - tmin) / tspan;
  x= px + dx * tn;
  y= py + dy * tn;
  pldj, x, y, x+vx, y+vy, legend=string(0), color=color;

  /*
   * Write tick labels.
   */
  x += 2 * vx;
  y += 2 * vy;
  t= swrite(format="%.3g", t);
  for (i = 0; i <= ni; i++) {
    plt, t(i), x(i), y(i), legend=string(0), justify=justify, tosys=1,
      color=color, font=font, height=height, opaque=opaque;
  }
}

/*----------------------------------------------------------------------*/
func win_copy_lim(src, dst)
/* DOCUMENT win_copy_lim, src, dst;
     Make limits of window DST the same as those in window SRC.

   SEE ALSO: current_window, limits, window. */
{
  win = current_window();
  window, src;
  l = limits();
  window, dst;
  limits,l;
  if (win >= 0) window, win; /* restore old current window */
}
/*----------------------------------------------------------------------*/

func xytitles(xtitle, ytitle, adjust)
/* DOCUMENT xytitles, xtitle, ytitle
       -or- xytitles, xtitle, ytitle, [deltax,deltay]
     Plot XTITLE horizontally under the viewport and YTITLE vertically
     to the left of the viewport.  If the tick numbers interfere with
     the labels, you can specify the [DELTAX,DELTAY] in NDC units to
     displace the labels.  (Especially for the y title, the adjustment
     may depend on how many digits the numbers on your scale actually
     have.)  Note that DELTAX moves YTITLE and DELTAY moves XTITLE.
     WARNING: There is no easy way to ensure that this type of title
              will not interfere with the tick numbering.  Interference
              may make the numbers or the title or both illegible.
   SEE ALSO: plt, pltitle
 */
{
  if (is_void(adjust)) adjust= _xytitles_adjust;
  port= viewport();
  if (xtitle && strlen(xtitle))
    plt, xtitle, port(zcen:1:2)(1), port(3)-0.050+adjust(2),
      font=pltitle_font, justify="CT", height=pltitle_height;
  if (ytitle && strlen(ytitle))
    plt, ytitle, port(1)-0.050+adjust(1), port(zcen:3:4)(1),
      font=pltitle_font, justify="CB", height=pltitle_height, orient=1;
}
_xytitles_adjust = [0.0,0.0];

/*----------------------------------------------------------------------*/

func pltr (text, rx, ry , legend=, hide=, color=, font=,
           height=, opaque=, orient=, justify=) {
/* DOCUMENT pltr, text, Rx, Ry;
   plot text on the sheet. Rx and Ry are in pourcentage unit.
   0% is the left/bottom corner 100% is the right/top corner of the
   system.
   KEYWORDS: legend, hide
              color, font, height, opaque, orient, justify
    SEE ALSO: plt, plt1, plg, plm, plc, plv, plf, pli, plt, pldj, plfp, pledit
              limits, range, fma, hcp, pltitle
    SEE ALSO:
 */
  V = viewport();
  x = (V(1)) + (V(2)-V(1))*(rx/100.);
  y = (V(3)) + (V(4)-V(3))*(ry/100.);
  plt, text, x, y, tosys=0, legend=legend, hide=hide,
    color=color, font=font, height=height, opaque=opaque,
    orient=orient, justify=justify;

}

/*----------------------------------------------------------------------*/

func plim(z,cmin=,cmax=) {
/* DOCUMENT plim, z

   SEE ALSO: pli
 */
  dims=dimsof(z);
  if(dims(1) != 2) error, "Expecting 2D array convertable to type double as argument";

  pli, z, 0.5, 0.5, dims(2)+0.5 , dims(3)+0.5, cmin=cmin, cmax=cmax;
}

/*----------------------------------------------------------------------*/

func xytitles2(xtitle, ytitle, adjust, height=, font=, orient=, pos=)
/* DOCUMENT xytitles, xtitle, ytitle
       -or- xytitles, xtitle, ytitle, [deltax,deltay]
     Plot XTITLE horizontally under the viewport and YTITLE vertically
     to the left of the viewport.  If the tick numbers interfere with
     the labels, you can specify the [DELTAX,DELTAY] in NDC units to
     displace the labels.  (Especially for the y title, the adjustment
     may depend on how many digits the numbers on your scale actually
     have.)  Note that DELTAX moves YTITLE and DELTAY moves XTITLE.
     WARNING: There is no easy way to ensure that this type of title
              will not interfere with the tick numbering.  Interference
              may make the numbers or the title or both illegible.
   SEE ALSO: plt, pltitle
 */
{
  if (is_void(height)) height = pltitle_height;
  if (is_void(adjust)) adjust = _xytitles_adjust;
  if (is_void(font))   font   = pltitle_font;
  if (is_void(orient)) orient = pltitle_orient;
  if (is_void(pos))    pos    = pltitle_pos;

  justify = [0,0];
  if (orient(1)==0 || orient(1)==2) justify(1) = pos(1)+4;
  else justify(1) = pos(1)+ [4,12,20](pos(1));

  if (orient(2)==1 || orient(2)==3) justify(2) = pos(2)+20;
  else   justify(2) = 3 + [20,12,4](pos(2));


  if (pos(1)==2) xindex = zcen; else xindex= pos(1)==1;
  if (pos(2)==2) yindex = zcen; else yindex= pos(2)==1;
  port=viewport();
  if (xtitle && strlen(xtitle))
    plt, xtitle, port(1:2)(xindex)(1), port(3)-0.050+adjust(2),
      font=font, justify=justify(1), height = height, orient=orient(1);
  if (ytitle && strlen(ytitle))
    plt, ytitle, port(1)-0.050+adjust(1), port(3:4)(yindex)(1),
      font = font, justify=justify(2), height=height, orient=orient(2);
}
_xytitles_adjust = [0.0,0.0];
pltitle_orient = [0,1];
pltitle_pos = [2,2];


func lxy(xflag,yflag) {
 /* DOCUMENT lxy, xflag, yflag
      sets the linear/log axis scaling flags for the current coordinate
      system.  XFLAG and YFLAG may be nil or omitted to leave the
      corresponding axis scaling unchanged, 0 to select linear scaling,
      or 1 to select log scaling.
    SEE ALSO: plsys, limits, range, plg, gridxy
    from logxy defined at:  LINE: 849  FILE: /home/cpinte/bin/yorick//share/yorick/2.1/i0/graph.i
 */
  logxy, xflag,yflag ;
  checklog ;
}


func plot_histo(y,bins,cumul=,color=,width=,type=) {
/* DOCUMENT plot_histo(y,bins,cumul=,color=,width=,type=)
   plot un histograme (cumule) des valeurs de "y" dans les bins "bins"
   max(bins) doit etre plus grand que max(y)
   SEE ALSO: histogram, digitize
 */
  if (is_void(cumul)) cumul = 0;
  if (is_void(color)) color="black";
  if (is_void(width)) width=1;
  if (is_void(type)) type=1;

  if (numberof(y)) {
    histo = histogram(digitize(y,bins),top=numberof(bins)) ;
    if (cumul) histo=(1.0 * histo(psum))/sum(histo) ;
    plh, histo, bins,color=color, width=width, type=type ;
  }
}

/************************************************************************/

func PlotArrow(x1, y1, x2, y2, width=, color=, type=, dx=, dy=)
    /* DOCUMENT yocoPlotArrow(x1, y1, x2, y2, width=, color=, type=, marks=, arlng=, facX=)

       DESCRIPTION
       Plots a simple arrow from one point to another.
       Works pretty much like pldj;

       PARAMETERS
       - x1   : start x
       - y1   : start y
       - x2   : end x
       - y2   : end y
       - width: see 'pldj'
       - color: see 'pldj'
       - type : see 'pldj'
       - marks: see 'pldj'
       - arlng: arrow length
       - facX : factor to plot the arrow with right proportions whatever
       the limits are. The default guess should work as soon as
       the limits are set before plotting the arrow.

       EXAMPLES
       > yocoGuiWinKill;
       > window,3;
       > yocoPlotArc, 0, 0, 10, 0, 7*pi/4, fill=1;
       > yocoPlotCircle, 30, 50, 10, fill=1;
       > yocoPlotArrow, 28, 38, 8, 8;

       SEE ALSO
    */
{
    if(is_void(facX))
    {
        vp = viewport();
        xyfac = (vp(2) - vp(1)) / (vp(4) - vp(3));
        lims = limits();
        facX = (lims(2) - lims(1)) / (lims(4) - lims(3)) / xyfac;
    }

    ang = atan(facX*(y2-y1),x2-x1);
    lng = abs(y2-y1,(x2-x1));

    if(is_void(arlng))
        arlng = 0.1*lng;

    pldj,x1,y1,x2,y2, width=width, color=color;

    pldj, x2, y2, x2-dx, y2-dy,
        width=width, color=color;
    pldj, x2, y2, x2-dx, y2+dy,
        width=width, color=color;
}

/************************************************************************/

func my_xytitles(xtitle, ytitle, adjust, color=)
/* DOCUMENT xytitles, xtitle, ytitle
       -or- xytitles, xtitle, ytitle, [deltax,deltay]
     Plot XTITLE horizontally under the viewport and YTITLE vertically
     to the left of the viewport.  If the tick numbers interfere with
     the labels, you can specify the [DELTAX,DELTAY] in NDC units to
     displace the labels.  (Especially for the y title, the adjustment
     may depend on how many digits the numbers on your scale actually
     have.)  Note that DELTAX moves YTITLE and DELTAY moves XTITLE.
     WARNING: There is no easy way to ensure that this type of title
              will not interfere with the tick numbering.  Interference
              may make the numbers or the title or both illegible.
   SEE ALSO: plt, pltitle
 */
{
  if (is_void(adjust)) adjust= _xytitles_adjust;
  if (is_void(color)) color="black" ;
  port= viewport();
  if (xtitle && strlen(xtitle))
    plt, xtitle, port(zcen:1:2)(1), port(3)-0.050+adjust(2),
      font=pltitle_font, justify="CT", height=pltitle_height, color=color;
  if (ytitle && strlen(ytitle))
    plt, ytitle, port(1)-0.050+adjust(1), port(zcen:3:4)(1),
      font=pltitle_font, justify="CB", height=pltitle_height, orient=1, color=color;
}
_xytitles_adjust = [0.0,0.0];
