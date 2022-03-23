func eps(name, pdf=, enlarge_bbox=)
/* DOCUMENT eps, name
     writes the picture in the current graphics window to the Encapsulated
     PostScript file NAME+".eps" (i.e.- the suffix .eps is added to NAME).
     This function requires ghostscript.  Any hardcopy file associated with
     the current window is first closed, but the default hardcopy file is
     unaffected.  As a side effect, legends are turned off and color table
     dumping is turned on for the current window.
     The external variable EPSGS_CMD contains the command to start
     ghostscript.
   SEE ALSO: pdf, png, jpeg, epsi, hcps, window, fma, hcp, hcp_finish, plg
 */
{
  if (strpart(name, -3:0) == ".eps") name = strpart(name,1:-4);
  /* dump the postscript file */
  psname = hcps(name+".pseps");

  /* begin copying to the eps file */
  f = create(name+".eps");
  g = open(psname);
  write, f, format="%s\n", "%!PS-Adobe-2.0 EPSF-1.2";
  rdline, g;
  line = rdline(g);
  if (strmatch(line,"% EPSF-3.0")) line = rdline(g); /* old ps.ps bug */
  for (i=1 ; i<=4 ; i++) {  /* Title For CreationDate Creator */
    write, f, format="%s\n", line;
    line = rdline(g);
  }

  /* use ghostscript to compute true bounding box */
  bbname = name+".bbeps";
  for (;;) {
    gscmd = EPSGS_CMD+" -sDEVICE=bbox -sOutputFile=- \"%s\" >>\"%s\" 2>&1";
    system, swrite(format=gscmd, psname, bbname);
    bb = rdline(open(bbname), 20);
    bb = bb(where(bb));

    remove, bbname;
    tok = strtok(bb);
    list = where(tok(1,) == "%%HiResBoundingBox:");
    if (!numberof(list)) {
      list = where(tok(1,) == "%%BoundingBox:");
    }
    xmn = ymn = xmx = ymx = 0.;
    if (!numberof(list) || sread(tok(2,list(1)), xmn, ymn, xmx, ymx) != 4) {
      /* Ghostscript 7.07 bbox fails if -dSAFER present,
       * Ghostscript 8.61 bbox fails if -dSAFER absent
       * the 8.61 failure gives an incorrect bounding box, so will never
       *   reach this workaround code
       * the 7.07 bug produces no BoundingBox comments at all, so will
       *   reach here
       * therefore, -dSAFER should be present in EPSGS_CMD initially
       *   and be removed for a second try with code that works in 7.07
       */
      if (strpart(EPSGS_CMD, -7:0) == " -dSAFER")
        EPSGS_CMD = strpart(EPSGS_CMD, 1:-8);
      else
        error, "ghostscript sDEVICE=bbox bug workaround failed";
    }
    break;
  }

  if (!is_void(enlarge_bbox)) {
    xmn -= enlarge_bbox ;
    xmx += enlarge_bbox ;
    ymn -= enlarge_bbox ;
    ymx += enlarge_bbox ;

    bb(1) = swrite(int(xmn), int(ymn), int(xmx), int(ymx), format= "%%%%BoundingBox: %3i %3i %3i %3i") ;
    bb(2) = swrite(xmn, ymn, xmx, ymx, format= "%%%%HiResBoundingBox: %10.6f %10.6f %10.6f %10.6f") ;
  }

  if (!pdf) {
    write, f, format="%s\n", bb;
    write, f, format="%s\n", "save countdictstack mark newpath "+
      "/showpage {} def /setpagedevice {pop} def";
  } else {
    /* concept from epstopdf perl script
     * by Sebastian Rahtz and Heiko Oberdiek,
     * distributed as part of the TeTeX package, see http://www.tug.org
     */
    write, f, format="%%BoundingBox: 0 0 %f %f\n", xmx-xmn, ymx-ymn;
    write, f, format="<< /PageSize [ %f %f ] >> setpagedevice\n",
      xmx-xmn, ymx-ymn;
    write, f, format="gsave %f %f translate\n",
      -xmn, -ymn;
  }
  write, f, format="%s\n", "%%EndProlog";
  while (line) {
    if (strpart(line,1:2)!="%%")
      write, f, format="%s\n", line;
    line = rdline(g);
  }

  close, g;
  remove, psname;

  write, f, format="%s\n", "%%Trailer";
  if (!pdf) {
    write, f, format="%s\n", "cleartomark "+
      "countdictstack exch sub { end } repeat restore";
  } else {
    write, f, format="%s\n", "grestore";
  }
  write, f, format="%s\n", "%%EOF";
  close, f;
  return name+".eps";
}

func pdf(name, enlarge_bbox=)
/* DOCUMENT pdf, name
     writes the picture in the current graphics window to the Adobe PDF
     file NAME+".pdf" (i.e.- the suffix .pdf is added to NAME).  The
     pdf file is intended to be imported into MS PowerPoint or other
     commercial presentation software, or into in pdftex or pdflatex
     documents; it is cropped.  The result should be equivalent to
     running the epstopdf utility (which comes with TeX, see www.tug.org)
     on the eps file produced by the eps command.
     This function requires ghostscript.  Any hardcopy file associated with
     the current window is first closed, but the default hardcopy file is
     unaffected.  As a side effect, legends are turned off and color table
     dumping is turned on for the current window.
     The external variable EPSGS_CMD contains the command to start
     ghostscript.
   SEE ALSO: eps, png, jpeg, hcps, window, fma, hcp, hcp_finish, plg
 */
{
  if (is_void(enlarge_bbox)) enlarge_bbox=3 ;
  if (strpart(name, -3:0) == ".pdf") name = strpart(name,1:-4);
  /* first run ghostscript to produce an eps translated to (0,0) */
  psname = eps(name+".pdf", pdf=1, enlarge_bbox=enlarge_bbox);
  /* second run ghostscript to produce the cropped pdf */
  gscmd = EPSGS_CMD+" -sDEVICE=pdfwrite -sOutputFile=\"%s\" \"%s\"";
  system, swrite(format=gscmd, name+".pdf", psname);
  remove, psname;
}
