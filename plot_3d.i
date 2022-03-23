local package_plot_3D;
/* DOCUMENT package_plot_3D -- plot_3D.i
   Simple package to plot 3-D lines, curves and text: pldj3d, plg3d,...
   See help, plot_3D for a complete documentation.
 */

/* Author : Jean-Baptiste LeBouwuin - jblebou@obs.ujf-grenoble.fr
 * $Log: plot_3D.i,v $
 * Revision 1.2  2005/02/11 11:11:23  jblebou
 * Major change: add the 'package_name' local and help at the front
 * of all files.
 *
 * Revision 1.1  2004/11/25 13:53:14  jblebou
 * Put in the CVS archive.
 *
 *
 */

if(_PRINT) write,"#include \"plot_3D.i\"";

/* ------------------------------------------------------------------- */
extern plot_3D;
/* DOCUMENT plot_3D

   Simple package to plot 3-D lines, curves and texts, compliant with
   the pl3d.i package.

   content:
    - pldj3d
    - plg3d
    - plt3d
    
   * $Log: plot_3D.i,v $
   * Revision 1.2  2005/02/11 11:11:23  jblebou
   * Major change: add the 'package_name' local and help at the front
   * of all files.
   *
   * Revision 1.1  2004/11/25 13:53:14  jblebou
   * Put in the CVS archive.
   *
     
   SEE ALSO: pldj, plg, plt, rot3, draw3...
 */


func pldj3d(xyz0, xyz1, legend=, hide=, type=, width=, color=, tosys=)
/* DOCUMENT pldj3d(xyz0, xyz1, legend=, hide=, type=, width=, color=, tosys=)

      Plots 3-D disjoint lines from xyz0 to xyz1. xyz0 and xyz1 should be
      3-D vectors but may have any additional dimensionality.
      
      All the keyword can be multi-dimentional if the customizations of the
      lines are differents. The 'tosys' keyword specify the Gist system.

      exemple:
      * plot a coordinate reference:
      > clear3;
      > pldj3d, 0, [[1,0,0],[0,1,0],[0,0,1]], color=["red","black","blue"];
      > orient3;
      > limits,-3,3,-3,3;
      > spin3, 40, [1,1,0], dtmin=0.2;

    The following keywords are legal (each has a separate help entry):
    KEYWORDS: legend, hide
              type, width, color, tosys
    SEE ALSO: pldj, plg3d, plt3d, rot3, range3, plwf
 */
{
  require, "pl3d.i";
  if(_draw3) {
    local arguments, lxyz0,lxyz1,llegend,lhide,ltype,lwidth,lcolor,ltosys,imax;
    arguments= xyz0;
    lxyz0    = _nxt(arguments);
    lxyz1    = _nxt(arguments);
    llegend  = _nxt(arguments);
    lhide    = _nxt(arguments);
    ltype    = _nxt(arguments);
    lwidth   = _nxt(arguments);
    lcolor   = _nxt(arguments);
    ltosys   = _nxt(arguments);
    
    imax = max(numberof(lxyz0(1,)), numberof(lxyz1(1,)),
               numberof(llegend), numberof(lhide), 
               numberof(ltype), numberof(lwidth),
               numberof(ltosys), numberof(lcolor));

    /* loop on the different realisation */
    for(i=1;i<=imax;i++) {
    local pxyz0,pxyz1,plegend,phide,ptype,pwidth,pcolor,ptosys;
    pxyz0= numberof(lxyz0(1,))>=1 ? lxyz0(,i%numberof(lxyz0(1,))) : lxyz0; 
    pxyz1= numberof(lxyz1(1,))>=1 ? lxyz1(,i%numberof(lxyz1(1,))) : lxyz1; 
    plegend= numberof(llegend)>=1 ? llegend(i%numberof(llegend)) : llegend; 
    phide= numberof(lhide)>=1 ? lhide(i%numberof(lhide)) : lhide; 
    ptype= numberof(ltype)>=1 ? ltype(i%numberof(ltype)) : ltype; 
    pwidth= numberof(lwidth)>=1 ? lwidth(i%numberof(lwidth)) : lwidth; 
    pcolor= numberof(lcolor)>=1 ? lcolor(i%numberof(lcolor)) : lcolor; 
    ptosys= numberof(ltosys)>=1 ? ltosys(i%numberof(ltosys)) : ltosys; 
    /* recover x,y,z in the viewer coordinate */
    get3_xy,pxyz0,px0,py0;
    get3_xy,pxyz1,px1,py1;

    /* plot the line */
    plsys,ptosys;
    pldj, px0, py0, px1, py1, legend=plegend, hide=phide, \
      type=ptype, width=pwidth, color=pcolor;

    }/* end loop */    
    return;  
  }
  
  /* defaut value */
  legend = string(legend);
  if(is_void(color)) color = "black";
  if(is_void(type))  type  = 1;
  if(is_void(width)) width = 1;
  if(is_void(tosys)) tosys = plsys();
  if(is_void(hide))  hide  = 0;
  
  set3_object, pldj3d, _lst(xyz0,xyz1,legend,hide,type,width,color,tosys);
  return;
}

/* ------------------------------------------------------------------- */

func plg3d(xyz, legend=, hide=, type=, width=, color=, closed=, smooth=,
           marks=, marker=, mspace=, mphase=,
           rays=, arrowl=, arroww=, rspace=, rphase=,tosys=)
/* DOCUMENT plg3d, xyz

      Plots a 3-D curves of consecutives points XYZ.
      XYZ must be an array of 3-D vectors:

      dimsof(xyz) = [.., 3 , npoints_per_curves ,  n_curves]

      The n_curves dimention is optianl, it allows to plot several 3-D curves
      at the same time. In this case, the keywords can be arrays to customize
      each curves independantly. The 'tosys' keyword specify the Gist system.

      example:
      * plot a serpentin:
      >
      > x = span(-2.5,2.5,100);
      > y = cos(2.*pi*x);
      > z = sin(2.*pi*x);
      > xyz = transpose([x,y,z]);
      >
      > window3;
      > gnomon, 1;
      > plg3d, xyz;
      >
      > limits,-5,5,-5,5;
      > spin3, 40, [1,1,0], dtmin=0.2;

      * plot a bundle of curves:
      >
      > x = array(span(-1,1,100),10);
      > y = indgen(10)(-:1:100,);
      > z = x^y;
      > y = y-5;
      > xyz = transpose([x,y,z], [1,2,3]);
      > 
      > window3;
      > clear3;
      > gnomon, 1;
      > plg3d, xyz, color=indgen(-4:-14:-1);
      > 
      > limits,-5,5,-5,5;
      > spin3, 40, [1,1,0], dtmin=0.2;
      
    The following keywords are legal (each has a separate help entry):
    KEYWORDS: legend, hide
              type, width, color, closed, smooth
              marks, marker, mspace, mphase
              rays, arrowl, arroww, rspace, rphase, tosys
    SEE ALSO:  plg, plg3d, plt3d, rot3, range3, plwf 
 */
{
  require, "pl3d.i";
  if(_draw3) {
    local arguments, lxyz, llegend, lhide, ltype, lwidth, lcolor, lclosed, lsmooth,
           lmarks, lmarker, lmspace, lmphase,
           lrays, larrowl, larroww, lrspace, lrphaseimax;
    arguments= xyz;
    lxyz    = _nxt(arguments);
    llegend = _nxt(arguments);
    lhide   = _nxt(arguments);
    ltype   = _nxt(arguments);
    lwidth  = _nxt(arguments);
    lcolor  = _nxt(arguments);
    lclosed = _nxt(arguments);
    lsmooth = _nxt(arguments);
    lmarks  = _nxt(arguments);
    lmarker = _nxt(arguments);
    lmspace = _nxt(arguments);
    lmphase = _nxt(arguments);
    lrays   = _nxt(arguments);
    larrowl = _nxt(arguments);
    larroww = _nxt(arguments);
    lrspace = _nxt(arguments);
    lrphaseimax = _nxt(arguments);    
    ltosys = _nxt(arguments);    
    
    imax = max(numberof(lxyz(1,1,)),
               numberof(llegend), numberof(lhide), numberof(ltype),
               numberof(lwidth), numberof(lcolor), numberof(lclosed),
               numberof(lsmooth), numberof(lmarks), numberof(lmarker),
               numberof(lmspace), numberof(lmphase), numberof(lrays),
               numberof(larrowl), numberof(larroww), numberof(lrspace),
               numberof(lrphaseimax), numberof(ltosys));

    /* loop on the different realisation */
    for(i=1;i<=imax;i++) {
    local pxyz, plegend, phide, ptype, pwidth, pcolor, pclosed, psmooth,
           pmarks, pmarker, pmspace, pmphase,
           prays, parrowl, parroww, prspace, prphaseimax;

    pxyz= numberof(lxyz(1,1,))>=1 ? lxyz(,,i%numberof(lxyz(1,1,))) : lxyz; 
    plegend= numberof(llegend)>=1 ? llegend(i%numberof(llegend)) : llegend; 
    phide= numberof(lhide)>=1 ? lhide(i%numberof(lhide)) : lhide; 
    ptype= numberof(ltype)>=1 ? ltype(i%numberof(ltype)) : ltype; 
    pwidth= numberof(lwidth)>=1 ? lwidth(i%numberof(lwidth)) : lwidth; 
    pcolor= numberof(lcolor)>=1 ? lcolor(i%numberof(lcolor)) : lcolor;
    pclosed= numberof(lclosed)>=1 ? lclosed(i%numberof(lclosed)) : lclosed;
    psmooth= numberof(lsmooth)>=1 ? lsmooth(i%numberof(lsmooth)) : lsmooth;
    pmarks= numberof(lmarks)>=1 ? lmarks(i%numberof(lmarks)) : lmarks;
    pmarker= numberof(lmarker)>=1 ? lmarker(i%numberof(lmarker)) : lmarker;
    pmspace= numberof(lmspace)>=1 ? lmspace(i%numberof(lmspace)) : lmspace;
    pmphase= numberof(lmphase)>=1 ? lmphase(i%numberof(lmphase)) : lmphase;
    prays= numberof(lrays)>=1 ? lrays(i%numberof(lrays)) : lrays;
    parrowl= numberof(larrowl)>=1 ? larrowl(i%numberof(larrowl)) : larrowl;
    parroww= numberof(larroww)>=1 ? larroww(i%numberof(larroww)) : larroww;
    prpsace= numberof(lrpsace)>=1 ? lrpsace(i%numberof(lrpsace)) : lrpsace;
    prphaseimax= numberof(lrphaseimax)>=1 ? lrphaseimax(i%numberof(lrphaseimax)) : lrphaseimax;
    ptosys= numberof(ltosys)>=1 ? ltosys(i%numberof(ltosys)) : ltosys;
    
    /* recover x,y,z in the viewer coordinate */
    get3_xy,pxyz,px,py;
    
    /* plot the line */
    plsys,ptosys;
    plg, py, px, legend=plegend, hide=phide, type=ptype, width=pwidth, color=pcolor, closed=pclosed, smooth=psmooth, marks=pmarks, marker=pmarker, mspace=pmspace, mphase=pmphase, rays=prays, arrowl=parrowl, arroww=parroww, rspace=prspace, rphase=prphase;

    }/* end loop */    
    return;  
  }
  
  /* defaut value */
  if(is_void(tosys)) tosys = plsys();
  
  set3_object, plg3d, _lst(xyz, legend, hide, type, width, color, closed, smooth, marks, marker, mspace, mphase, rays, arrowl, arroww, rspace, rphase,tosys);
  return;
}

/* ------------------------------------------------------------------- */

func plt3d(text, xyz, legend=,hide=,color=,font=,height=,opaque=,orient=,justify=,tosys=)
/* DOCUMENT plt3d, text, xyz
   
     Plots TEXT (a string) at the position XYZ (a 3-D vector).
     The apparent orientation of TEXT in the sheet is determine
     only by the 'orient='and 'justify=' keyword... i.d. the text orientation
     does not follow the rotations of the coordinate instead of the text
     position.

     TEXT and/or XYZ can be array of string and 3-D vectors. In this case the
     keyword can also be arrays in order to customize independantly the
     different string ouputs.

     See 'plt' for an explanation of the string format.

     * plot a coordinate reference:
     >
     > clear3;
     > pldj3d, 0, [[1,0,0],[0,1,0],[0,0,1]], color=["red","black","blue"];
     > 
     > plt3d,"Center", [0,0,0],height=16,font="helveticaB";     
     > plt3d,["X","Y","Z"],[[1,0,0],[0,1,0],[0,0,1]];
     > 
     > orient3;
     > limits,-2,2,-2,2;
     > spin3, 40, [1,1,0], dtmin=0.2;
     
   The following keywords are legal (each has a separate help entry):
   KEYWORDS: legend, hide
              color, font, height, opaque, orient, justify
              
   SEE ALSO: plt, plg3d, plt3d, rot3, range3, plwf 
 */
{
  require, "pl3d.i";
    if(_draw3) {
    local arguments,ltext,lxyz,llegend,lhide,lcolor,lfont,lheight,lopaque,lorient,ljustify,ltosys;
    arguments= text;
    ltext = _nxt(arguments);
    lxyz    = _nxt(arguments);
    llegend = _nxt(arguments);
    lhide   = _nxt(arguments);
    lcolor  = _nxt(arguments);
    lfont   = _nxt(arguments);
    lheight = _nxt(arguments);
    lopaque = _nxt(arguments);
    lorient = _nxt(arguments);
    ljustify = _nxt(arguments);
    ltosys   = _nxt(arguments);    
   
    imax = max(numberof(ltext),numberof(lxyz(1,)),
               numberof(llegend), numberof(lhide), numberof(lfont),
               numberof(lheight),numberof(lopaque),numberof(lorient),numberof(ljustify),numberof(ltosys));

    /* loop on the different realisation */
    for(i=1;i<=imax;i++) {
    local ptext,px,py,plegend,phide,pfont,pheight,popaque,porient,pjustify,ptosys;

    ptext= numberof(ltext)>=1 ? ltext(i%numberof(ltext)) : ltext; 
    pxyz= numberof(lxyz(1,))>=1 ? lxyz(,i%numberof(lxyz(1,))) : lxyz; 
    plegend= numberof(llegend)>=1 ? llegend(i%numberof(llegend)) : llegend; 
    phide= numberof(lhide)>=1 ? lhide(i%numberof(lhide)) : lhide;
    pfont= numberof(lfont)>=1 ? lfont(i%numberof(lfont)) : lfont;
    pheight= numberof(lheight)>=1 ? lheight(i%numberof(lheight)) : lheight;
    pcolor= numberof(lcolor)>=1 ? lcolor(i%numberof(lcolor)) : lcolor;
    popaque= numberof(lopaque)>=1 ? lopaque(i%numberof(lopaque)) : lopaque;
    porient= numberof(lorient)>=1 ? lorient(i%numberof(lorient)) : lorient;
    pjustify= numberof(ljustify)>=1 ? ljustify(i%numberof(ljustify)) : ljustify;
    ptosys= numberof(ltosys)>=1 ? ltosys(i%numberof(ltosys)) : ltosys;
    
    /* recover x,y,z in the viewer coordinate */
    get3_xy,pxyz,px,py;
    
    /* plot the text */
    plsys,ptosys;
    plt, ptext, px, py, legend=plegend, hide=phide, color=pcolor, font=pfont, height=pheight, opaque=popaque, orient=porient, justify=pjustify, tosys=ptosys;

    }/* end loop */    
    return;  
  }
  
  /* defaut value */
  if(is_void(tosys)) tosys = plsys();
  if(is_void(justify)) justify = "CC";
  
  set3_object, plt3d, _lst(text, xyz, legend, hide, color, font, height, opaque, orient, justify, tosys);
  return;
}





//==============================================================

func plwf2(data, xyz, legend=, cmin=, cmax=, edges=, tosys=, hide=, ncorner=, mode=)
{
  require, "pl3d.i";
  if(_draw3) {
    local ldata,lxyz,llegend,lcmin,lcmax,ledges,ltosys,lhide,imax,ldim,lncorner,lmode;
    arguments= data;
    ldata    = _nxt(arguments);
    lxyz     = _nxt(arguments);
    llegend  = _nxt(arguments);
    lcmin    = _nxt(arguments);
    lcmax    = _nxt(arguments);
    ledges   = _nxt(arguments);
    ltosys   = _nxt(arguments);
    lhide    = _nxt(arguments);
    lmode    = _nxt(arguments);

    imax = max(numberof(ldata(1,)), numberof(lxyz(1,1,1,)),
               numberof(lcmin), numberof(lcmax),
               numberof(ledges), numberof(ltosys), numberof(hide));
    /* loop on the different realisation */
    for(i=1;i<=imax;i++) {
      local pldata,plxyz,pllegend,plcmin,plcmax,pledges,pltosys,listc,vlistec;
      pldata= numberof(ldata(1,))>=1 ? ldata(,i%numberof(ldata(1,)))  : ldata;
      plxyz= numberof(lxyz(1,1,1,))>=1 ? lxyz(,,,i%numberof(lxyz(1,1,1,))) : lxyz;
      pllegend= numberof(llegend)>=1 ? llegend(i%numberof(llegend)) : llegend;
      plcmin= numberof(lcmin)>=1 ? lcmin(i%numberof(lcmin)) : lcmin;
      plcmax= numberof(lcmax)>=1 ? lcmax(i%numberof(lcmax)) : lcmax;
      pledges= numberof(ledges)>=1 ? ledges(i%numberof(ledges)) : ledges;
      pltosys= numberof(ltosys)>=1 ? ltosys(i%numberof(ltosys)) : ltosys;
      plhide= numberof(lhide)>=1 ? lhide(i%numberof(lhide)) : lhide;
      /* recover x,y,z in the viewer coordinate */
      get3_xy,plxyz,xc,yc,zc,1;
      plxyz = transpose([xc,yc,zc],2);
      /* if mode, remove the background */
      if(lmode) {
        listec = where(zc(avg,)>0.);
        plxyz  = plxyz(,,listec);
        pldata = pldata(listec);
      }
      /* sort in the back to front order */
      listec = sort(plxyz(3,avg,));
      plxyz  = plxyz(,,listec);
      pldata = pldata(listec);
      /* align the 2 dim*/
      nplxyz = array(dimsof(plxyz)(3), dimsof(plxyz)(4));
      plxyz  = ad_2_dim(plxyz,2);
      /* plot the i-eme data */
      plsys,pltosys;
      plfp, pldata, plxyz(2,), plxyz(1,), nplxyz,
        edges=pledges,cmax=plcmax,cmin=plcmin,legend=pllegend,hide=plhide;
    }
  return;
  }
  legend = string(legend);
  if(is_void(cmin))  cmin  = min(data);
  if(is_void(cmax))  cmax  = max(data);
  if(is_void(edges)) edges = 1;
  if(is_void(tosys)) tosys = plsys();
  if(is_void(hide))  hide  = 0;
  if(is_void(mode))  mode=1;
  set3_object, plwf2, _lst(data,xyz,legend,cmin,cmax,edges,tosys,hide,mode);
}
//*************************************************************//
func ad_2_dim(data,where)
{
  dim = dimsof(data);
  if(where+1>dim(1)) error,"not a confomable dimension";

  num = dim(where+1) * dim(where+2);

  if(where==1) dim1=[0];
  else {
    dim1    = dim(1:where);
    dim1(1) = where-1;
  }
  if(where+2>dim(1)) dim2=[0];
  else {
    dim2 = dim(where+3:0);
    dim2 = grow(numberof(dim2),dim2);
  }
  
  _data = array(structof(data),dim1,num,dim2);
  _data(*) = data(*);
  return _data;
}

 
func cube3 ( xyz, a) {
  if (is_void(xyz)) xyz = [0,0,0];
  if (is_void(a)) a=1;
  a=a/2.;
  return ([ [[-1,1,-1],[-1,-1,-1],[1,-1,-1],[1,1,-1]],
           [[-1,1,-1],[-1,1,1],[1,1,1],[1,1,-1]],
           [[1,1,1],[1,-1,1],[1,-1,-1],[1,1,-1]],
           [[-1,1,-1],[-1,-1,-1],[-1,-1,1],[-1,1,1]],
           [[-1,-1,1],[-1,-1,-1],[1,-1,-1],[1,-1,1]],
            [[-1,1,1],[-1,-1,1],[1,-1,1],[1,1,1] ]]*
          a + xyz(,-:1:4,-:1:6,))(,,*);
}
func carre3 ( xyz, a) {
  if (is_void(xyz)) xyz = [0,0,0];
  if (is_void(a)) a=1;
  a=a/2.;
  return ([ [[-1,1,-1],[-1,-1,-1],[1,-1,-1],[1,1,-1]]
           ]*a + xyz(,-:1:4,-:1:6,))(,,*);
          
}
