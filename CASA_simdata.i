// See ~/Software/CASA/CASA_modif.txt

func CASA_simdata(model,i,obstime=,config=,resol=,sampling_time=,pwv=,az=,decl=,distance=,phase_noise=,name=,iTrans=,rt=,only_prepare=,interferometer=,mosaic=,mapsize=,channels=,width=,correct_flux=) {
/* DOCUMENT  CASA version 4.2.2
   CASA_simdata(model,i,obstime=,config=,sampling_time=,pwv=,az=,decl=,distance=,,phase_noise=,name=,iTrans=,rt=,only_prepare=,interferometer=,mosaic=,mapsize=,channels=,width=,correct_flux=)
   prepare un model pour le simulateur ALMA de CASA :
   - taille des pixels
   - Jy/pixels
   - frequency

   - width in Ghz

   Si resol est donne (en arcsec), simdata calcule la configuration correspondante

   Genere un fichier fits et un script pour simdata

   Cycle0 : configs :   "cycle0.compact"  et "cycle0.extended"


   SEE ALSO:
 */

  is_image = ( nameof(structof(model)) == "McfostImage") ;

  workdir="CASA/" ;
  _CASA_clean, workdir ;

  if (is_void(rt)) rt=0 ;
  if (is_void(only_prepare)) only_prepare=0 ;
  if (is_void(interferometer)) interferometer = "alma" ;
  if (is_void(mosaic)) mosaic=0 ;
  if (is_void(correct_flux)) correct_flux=1 ;

  if (!is_image) {
    if (is_void(iTrans)) {
      write, "ERROR: iTrans must be specified for line data. Exiting." ;
      return [] ;
    }

    nTrans = numberof(model.freq) ;
    if (iTrans > nTrans) {
      write, "ERROR: iTrans is not in the computed range. Exiting" ;
      return [] ;
    }
  }

  // Name
  if (is_void(name)) name="simu" ;
  Name = name ;

  // PA du disque
  if (is_void(az)) {
    az = 1 ;
  }

  // Declinaison de la source
  if (is_void(decl)) {
    decl = "-22d59m59.8" ;
    write, "Forcing declination to "+decl ;
  }

  // Config
  if (is_void(config) && is_void(resol)) {
    "ERROR: config needed" ;
    error;
    return [] ;
  }

  if (is_void(config)) {
    sresol = swrite(resol,format="%2.2f") ;
    resol_name = "_resol="+sresol ;
    // JFG 20140311
    sresol = swrite(resol,format="%6.6f") ;
    resol_name2 = "alma_"+sresol+"arcsec" ;
  } else {
    if (typeof(config)=="long") {
      sconfig2 = swrite(config,format="%2i") ;
      sconfig2 = streplace(sconfig2,strgrep(" ",sconfig2),"0");
      sconfig = "out"+sconfig2 ;
    } else {
      //if ( (config=="early.125m") || (config=="early.250m")) {
      if (config=="cycle1_out9_tmp") {
        "Config cycle 1" ;
      }  else if ( (config=="cycle0.compact") || (config=="cycle0.extented")) {
        "Early config" ;
      } else if ( (config=="csv.mid") || (config=="csv.late") ) {
        "Science Verification config" ;
      }
      sconfig2 = config ;
      sconfig = sconfig2 ;
    }
    resol_name = "_config="+interferometer+"."+sconfig2 ;
    resol_name2 = sconfig2 ;

  }


  // Temps obs
  if (is_void(obstime)) {
    "ERROR: obstime needed" ;
    return [] ;
  }

  n_chiffres = floor(log10(obstime) + 1);
  obstime_format = "%"+swrite(format="%1i",int(n_chiffres))+"i";
  sobstime = swrite(obstime,format=obstime_format) ;

  if (is_void(sampling_time)) {
    "ERROR: sampling time needed" ;
    return [] ;
  }
  n_chiffres = floor(log10(sampling_time) + 1);
  sampling_format = "%"+swrite(format="%1i",int(n_chiffres))+"i";
  ssampling_time = swrite(sampling_time,format=sampling_format) ;

  if (is_void(pwv) || pwv < 1e-6) {
    "pwv not specified --> No thermal noise" ;
    th_noise = "''" ;
    lth_noise = 0 ;
  } else {
    th_noise = "'tsys-atm'" ;
    spwv = swrite(pwv, format="%4.2f") ;
    lth_noise = 1 ;
  }

  // Phase noise
  if (is_void(phase_noise)) phase_noise = 0

  // Distance
  if (is_void(distance)) {
    distance =  model.P.map.dist ;
  }

  // Pixel scale
  pixel_scale = 2.0*model.P.map.size_neb / (model.P.map.nx * distance * model.P.map.zoom) ;

  arcsec_to_deg = 1./3600. ;
  pixel_scale_x = pixel_scale * arcsec_to_deg ;
  pixel_scale_y = pixel_scale * arcsec_to_deg ;


  // frequency
  if (is_image) {
    freq = SI.c / (model.lamb * 1e-6) * 1e-9 ;
    inwidth = 8 ; // 8 Ghz par defaut en continu
    inchan = 1 ; // 1 channel pour continu
  } else {
    dv = (model.P.mol.vmax / model.P.mol.n_speed * SI.km)(1) ; // m.s^-1
    nu = (model.freq(iTrans)) ;// Hz
    freq = nu  * 1e-9 ; // GHz
    dnu = dv/SI.c * nu ; //  Hz

    inwidth =  dnu * 1e-9 ; // Ghz
    inchan = 2 * model.P.mol.n_speed + 1 ;

    if (!is_void(channels)) {
      inchan = numberof(channels) ;
    } else {
      channels = indgen(inchan) ;
    }

    if (!is_void(width)) {
      inwidth = width ;
    }

  }

  // Flux en Jy/pixel
  if (is_image) {
    if (rt) {
      image = model.image_rt(,,i,az,1) * (model.lamb*1e-6/ SI.c) * 1e26 ;
    } else {
      image = model.image(,,i,az,1) * (model.lamb*1e-6/ SI.c) * 1e26 ;
    }
    image = [[float(image)]] ; // Adding spectral & pola dimensions
  } else {
    image = model.Fline(,,channels,iTrans,i) / nu * 1e26 ;
  }

  // Correction distance eventuelle
  image = image * (model.P.map.dist/distance)^2 ;

  // Correction evntuelle flux
  image = image * correct_flux ;

  // Cleanning
  n_iter = 1000000 ;

  n_chiffres = floor(log10(max(n_iter,1)) + 1);
  n_iter= swrite(n_iter,format="%"+swrite(format="%1i",int(n_chiffres))+"i");

  // TMP : point source
  // etoile = image(101,101) ;
  // image = image * 0 ;
  // image(101,101) = etoile * 1e3;
  // FIN TMP


  n = model.P.map.nthet ;
  inc = acos((n-i+0.5)/n)*180/pi ;
  sinc = swrite(inc,format="%2.1f") ;

  n_chiffres = floor(log10(distance) + 1);
  format = "%"+swrite(format="%1i",int(n_chiffres))+"i";
  sdistance = swrite(int(distance),format=format) ;

  if (mosaic) {
    if (is_void(mapsize)) {
      write, "ERROR: mapsize required in mosaic mode" ;
      write, "Exiting" ;
      return ;
    } else {
      n_chiffres = floor(log10(mapsize) + 1);
      format = "%"+swrite(format="%1i",int(n_chiffres))+"i";
      smapsize = swrite(int(mapsize),format=format) ;
    }
    map_name = "_mapsize="+smapsize+"arcsec" ;
  } else {
    map_name = "" ;
  }

  // File names
  if (is_image) {
    slambda =  swrite(model.lamb,format="%3.0f") ;
    simu_name =  Name+"_dist="+sdistance+"_i="+sinc+"_lambda="+slambda+resol_name+"_obstime="+sobstime+"_decl="+decl+map_name ;
  } else {
    sfreq =  swrite(freq,format="%6.3f") ;
    simu_name =  Name+"_dist="+sdistance+"_i="+sinc+"_freq="+sfreq+"GHz"+resol_name+"_obstime="+sobstime+"_decl="+decl+map_name ;
  }

  // Header du fichier fits
  fh = cfitsio_open(workdir+simu_name+".raw.fits","w",overwrite=1); // create the file
  cfitsio_add_image, fh, image, "IMAGE";        // add an image

  cfitsio_write_key, fh, "CTYPE1", "RA---TAN";
  cfitsio_write_key, fh, "CRVAL1", 0. ;
  cfitsio_write_key, fh, "CRPIX1", model.P.map.nx/2+1 ;
  cfitsio_write_key, fh, "CDELT1", pixel_scale_x ;

  cfitsio_write_key, fh, "CTYPE2", "DEC--TAN";
  cfitsio_write_key, fh, "CRVAL2", 0. ;
  cfitsio_write_key, fh, "CRPIX2", model.P.map.ny/2+1 ;
  cfitsio_write_key, fh, "CDELT2", pixel_scale_y ;


  if (is_image) {
    cfitsio_write_key, fh, "RESTFREQ", float(freq)*1e9, "Hz" ;

    // 3eme axe
    cfitsio_write_key, fh, "CTYPE3","STOKES" ;
    cfitsio_write_key, fh, "CRVAL3",1.0 ;
    cfitsio_write_key, fh, "CDELT3",1.0 ;
    cfitsio_write_key, fh, "CRPIX3",1 ;

   // 4eme axe
    cfitsio_write_key, fh, "CTYPE4","FREQ" ;
    cfitsio_write_key, fh, "CRVAL4",freq*1e9,"Hz" ;
    cfitsio_write_key, fh, "CDELT4",2e9,"Hz" ;// 2GHz by default
    cfitsio_write_key, fh, "CRPIX4",0 ;

  } else {
    // TODO : 3eme axe pour raie

    cfitsio_write_key, fh, "CTYPE3", "VELO-LSR";
    cfitsio_write_key, fh, "CRVAL3", 0. ; // V au centre de la raie
    cfitsio_write_key, fh, "CRPIX3", inchan ; // centre de la raie
    cfitsio_write_key, fh, "CDELT3",  dv ; // TODO

    /*
    cfitsio_write_key, fh, "CTYPE3", "FREQ";
    cfitsio_write_key, fh, "CRVAL3", nu ;
    cfitsio_write_key, fh, "CRPIX3", model.P.mol.n_speed+1 ; // centre de la raie
    cfitsio_write_key, fh, "CDELT3",  dnu ; // TODO
    */
    cfitsio_write_key, fh, "RESTFREQ", nu ;
  }

  cfitsio_write_key, fh, "BUNIT", "JY/PIXEL";  // add a card in the first HDU
  cfitsio_write_key, fh, "BTYPE","Intensity";

  cfitsio_write_key, fh, "phase_noise", phase_noise ;

  cfitsio_write_date,fh;                       // write the date
  cfitsio_close,fh;                            // close the file

  // Fichier de config pour simdata
  f = open(workdir+simu_name+".py","w") ;
  //write, f, "default(simdata)\n", format="%s" ;

  //--------
  write, f, "project = 'DISK'\n", format="%s" ;
  write, f, "dryrun = False\n", format="%s" ;

  //--------
  write, f, "modifymodel = True\n", format="%s" ;
  write, f, "skymodel = '"+simu_name+".raw.fits'\n", format="%s" ;
  write, f, "inbright = 'unchanged'\n", format="%s" ;
  // ascension droite ne change rien --> OK
  write, f, "indirection = 'J2000 18h00m00.02 "+decl+"' # mosaic center, or list of pointings\n",format="%s";
  incell = swrite(pixel_scale,format="%5.4f") ;
  write, f, "incell = '"+incell+"arcsec'\n",format="%s";
  write, f, "incenter = '"+swrite(freq,format="%17.15e")+"GHz'\n", format="%s" ;
  write, f, "inwidth = '"+swrite(float(inwidth),format="%17.15e")+"GHz'\n",format="%s";
  write, f, "inchan = "+swrite(inchan,format="%1i")+"\n",format="%s";
  //  write, f, "ignorecoord = True\n", format="%s" ;

  //--------
  write, f, "setpointings = True\n", format="%s" ;

  if (mosaic) {
    write, f, "mapsize ='"+smapsize+"arcsec'\n",format="%s" ;
    write, f, "maptype = 'hexagonal'\n", format="%s" ;
    write, f, "pointingspacing = ''\n", format="%s" ;
  } else {
    write, f, "mapsize = ''\n", format="%s" ;
    write, f, "pointingspacing = '1.0arcmin'\n",format="%s";
  }

  write, f, "totaltime = '"+sobstime+"s'\n", format="%s" ;
  write, f, "integration = '"+ssampling_time+"s'\n", format="%s" ;

  //--------
  write, f, "predict = True\n", format="%s" ;
  write, f, "complist = ''\n", format="%s" ;
  write, f, "refdate = '2012/06/21/03:25:00'\n", format="%s" ;
  write, f, "repodir=os.getenv(\"CASAPATH\").split(\' \')[0]\n", format="%s" ;
  if (resol) {
    write, f, "antennalist = \"alma;%farcsec\" \% "+sresol+"\n", format="%s" ;
  } else {
    //write, f, "antennalist = repodir+'/data/alma/simmos/"+interferometer+"."+sconfig+".cfg'\n", format="%s" ;
    write, f, "antennalist = repodir+'/data/alma/simmos/"+sconfig+".cfg'\n", format="%s" ;
  }
  //write, f, "checkinputs = 'yes'\n", format="%s" ;

  //--------
  write, f, "thermalnoise = "+th_noise+"\n",format="%s";
  if (lth_noise) {
    write, f, "user_pwv = "+spwv+"\n",format="%s";
    write, f, "vis = project+'.noisy.ms' # clean the data with *thermal noise added*\n", format="%s";
  }

  //--------
  write, f, "image = True\n", format="%s" ;
  write, f, "cleanmode = 'clark'\n",format="%s";
  imsize = swrite(nx,format="%3i") ;
  write, f, "imsize = ["+imsize+","+imsize+"]\n",format="%s";
  write, f, "cell = ''\n", format="%s" ;
  write, f, "niter = "+n_iter+"\n",format="%s";
  write, f, "threshold = '0.0mJy'\n",format="%s";
  write, f, "weighting = 'natural'\n",format="%s";
  write, f, "outertaper = []\n",format="%s";
  write, f, "stokes = 'I'\n",format="%s";


  //--------
  write, f, "analyze = True\n", format="%s" ;

  //--------
  write, f, "graphics = 'file'\n",format="%s";
  write, f, "overwrite = True\n",format="%s";
  write, f, "verbose = False\n",format="%s";
  write, f, "async = False\n",format="%s";

  // ---- choix plots :
  /*
  write, f, "showarray = False\n",format="%s";
  write, f, "showuv = False\n",format="%s";
  write, f, "showpsf = False\n",format="%s";
  write, f, "showmodel = False\n",format="%s";
  write, f, "showconvolved = False\n",format="%s";
  write, f, "showclean = False\n",format="%s";
  write, f, "showresidual = False\n",format="%s";
  write, f, "showdifference = False\n",format="%s";
  write, f, "showfidelity = False\n",format="%s";
  */


  //write, f, "inp()\n",format="%s";
  write, f, "simalma()\n",format="%s";
  write, f, "pl.savefig('"+simu_name+".png')\n",format="%s" ;
  //  write, f, "exportfits(imagename=project+'/'+project+'."+resol_name+".image',fitsimage='ALMA_disk.fits')\n",format="%s";
  // JFG 20140311
  write, f, "exportfits(imagename=project+'/'+project+'."+resol_name2+".noisy.image',fitsimage='"+simu_name+".fits')\n",format="%s";

  if (phase_noise) {
    sdelta_pwv = swrite(0.15*pwv, format="%4.2f") ;

    write, f, "print '[yorick] Adding phase noise pwv="+spwv+" delta_pwv="+sdelta_pwv+"'\n", format="%s" ;
    // Adding phase noise
    //    write, f, "os.system('cp -r '+project+'/'+project+'.noisy.ms '+project+'/'+project+'.phase_noise.ms')\n",format="%s";
    // JFG 20140311
    write, f, "os.system('cp -r '+project+'/'+project+'."+resol_name2+".ms '+project+'/'+project+'.phase_noise.ms')\n",format="%s";
    write, f, "sm.openfromms(project+'/'+project+'.phase_noise.ms')\n",format="%s";
    write, f, "sm.settrop(mode='screen',pwv="+spwv+",deltapwv="+sdelta_pwv+",beta=1.1,windspeed=7.0)\n",format="%s";
    write, f, "sm.corrupt()\n",format="%s";
    write, f, "sm.done()\n",format="%s";
    write, f, "print '[yorick] cleaning phase-noised image'\n", format="%s" ;
    //    write, f, "execfile(project+'/'+project+'.clean.last')\n", format="%s" ;
    // JFG 20140311
    write, f, "execfile(project+'/'+project+'."+resol_name2+".clean.last')\n", format="%s" ;
    write, f, "vis=project+'/'+project+'.phase_noise.ms'\n", format="%s" ;
    write, f, "imagename=project+'/'+project+'.phase_noise'\n", format="%s" ;
    write, f, "clean()\n", format="%s" ;
    write, f, "print '[yorick] DONE'\n", format="%s" ;

    write, f, "exportfits(imagename=project+'/'+project+'.phase_noise.image',fitsimage='ALMA_disk_phase-noise.fits')\n", format="%s";
  }
  write, f, "exit\n", format="%s" ;
  close, f ;

  // Saving input files
  //system, "cp CASA/disk.fits CASA/"+simu_name+".raw.fits" ;
  //system, "cp CASA/disk.py CASA/"+simu_name+".py" ;

  if (!only_prepare) {
   simu_name =  _run_CASA(simu_name) ;
   return simu_name
  } // !only_prepare

}

//**************************************************

func _run_CASA(simu_name,nodedir=) {

  tic ;
  "Starting yoCASAPY ..." ;

  if (is_void(nodedir)) {
    nodedir="" ;
  } else {
    nodedir=nodedir+"/" ;
  }
  workdir = "CASA/"+nodedir ;
  homedir = get_cwd() ;

  _CASA_clean(workdir) ;


  // Do we run the simulator with phase noise ?
  fh = cfitsio_open(workdir+simu_name+".raw.fits","r");
  phase_noise = cfitsio_get(fh, "phase_noise") ;
  cfitsio_close,fh;

  // Running the simulator
  cd, workdir ;

  //system, "casa --nogui -c "+simu_name+".py" ;
  system, "/Applications/CASA.app/./Contents/Resources/python/regressions/admin/runcasa_from_shell.sh 0 "+simu_name+".py"
  cd, homedir ;


  // Extract beam info and add it to the fits files
  system, "rm -rf "+workdir+"/beam.txt" ;
  //system, "grep 'Beam fit' CASA/casapy.log | awk '{print $7 \" \" $9 \" \" $13}' > CASA/beam.txt" ;
  system, "grep 'Fitted beam' "+workdir+"casapy*.log | awk '{print $10 \" \" $12 \" \" $16}' > "+workdir+"beam.txt" ;
  beam = OpenASCII(workdir+"beam.txt",valtype=float,prompt=0) ;
  system, "rm -rf "+workdir+"beam.txt" ;


  system, "mv "+workdir+"ALMA_disk.png  "+workdir+simu_name+".png" ;
  system, "mv "+workdir+"ALMA_disk.fits "+workdir+simu_name+".fits" ;
  //  system, "mv "+workdir+"casapy.log "+workdir+simu_name+".log" ;
  // JFG 20140311
  system, "mv "+workdir+"casapy*.log "+workdir+simu_name+".log" ;
  if (phase_noise) system, "mv "+workdir+"ALMA_disk_phase-noise.fits "+workdir+simu_name+"_phase-noise.fits" ;


  // Adding the beam information in the output fits file
  fh = cfitsio_open(workdir+simu_name+".fits","a") ;
  cfitsio_write_key, fh,"BEAM1",beam.X1, "beam major axis [arcsec]";
  cfitsio_write_key, fh,"BEAM2",beam.X2, "beam major axis [arcsec]";
  cfitsio_write_key, fh,"BEAM_PA",beam.X3, "beam PA [degrees]";
  //cfitsio_write_key, fh,"distance",distance, "[pc]";
  cfitsio_close, fh ;

  if (phase_noise) {
    fh = cfitsio_open(""+workdir+simu_name+"_phase-noise.fits","a") ;
    cfitsio_write_key, fh,"BEAM1",beam.X1, "beam major axis [arcsec]";
    cfitsio_write_key, fh,"BEAM2",beam.X2, "beam major axis [arcsec]";
    cfitsio_write_key, fh,"BEAM_PA",beam.X3, "beam PA [degrees]";
    //cfitsio_write_key, fh,"distance",distance, "[pc]";
    cfitsio_close, fh ;
  }


  _CASA_clean(workdir) ;

  "yoCASAPY DONE" ;
  //write, "Simulation done in ",tac(), " sec" ;

  return simu_name ;
}

//**************************************************

func _CASA_clean(workdir) {
  // Nettoyage
  system, "rm -rf "+workdir+"DISK* "+workdir+"disk.fits "+workdir+"*.last "+workdir+"disk.py "+workdir+"ALMA_disk.png "+workdir+"ALMA_disk.fits "+workdir+"*.log*" ;

}

//**************************************************

func beam2config(beam,lambda) {
/* DOCUMENT
     Renvoie la config ALMA correspondant a une resolution donnee
     beam is FWHM en arcsec
     lambda en micron
     GROSSE APPROX !!
     SEE ALSO:
 */

  nu = SI.c / (lambda * 1e-6) ;

  config = floor(0.5 - 13.72 * log10(beam * nu/672e9)) ;


  if (config < 1) {
    "Warning: config < 1" ;
    config = 1 ;
  }
  if (config > 28) {
    "Warning: config > 28" ;
    config = 28 ;
  }

  return config ;
}

//**************************************************

func config2beam(config,lambda) {
/* DOCUMENT
     Renvoie le beam ALMA (FWHM en arcsec) correspondant a une config donnee
     lambda en micron
     GROSSE APPROX !!
   SEE ALSO:
 */

  nu = SI.c / (lambda * 1e-6) ;

  Beam = 672e9/nu * 10^(( - (config+0.5))/13.72) ;

  return Beam ;

}

//**************************************************

func plot_CASA_image(filename,cmin=,cmax=,stats=,title=,negative=) {

  //window, 0, style="boxed.gs" ; fma ;

//  pltitle_height = 14 ;
//  get_style, landscape,systems,legends,clegends ;
//  systems.ticks.horiz.textStyle.height=0.015;
//  systems.ticks.vert.textStyle.height=0.015;
//  set_style, landscape,systems,legends,clegends ;

  if (is_void(negative)) negative=0 ;

  // Lecture
  fh = cfitsio_open(filename) ;

  CDELT1 = cfitsio_get(fh,"CDELT1") ;
  CDELT2 = cfitsio_get(fh,"CDELT2") ;
  pixel_scale = abs(CDELT1) * 3600;
  distance = cfitsio_get(fh,"distance") ;

  a = cfitsio_get(fh,"BEAM1",comments)  ;
  b = cfitsio_get(fh,"BEAM2",comments) ;
  theta = cfitsio_get(fh,"BEAM_PA",comments) ;  if (!is_void(theta)) theta *= pi/180 ;

  img   = cfitsio_read_image(fh);
  cfitsio_close, fh ;

  // Plot image
  nx = dimsof(img)(2) ;
  ny = nx ;
  img_size = nx * pixel_scale ;
  img = img(,,1) ;

  if (!is_void(stats)) stat, img ;


  if (negative) {
    if (is_void(cmax)) {
      Cmin = -max(img) ;
    } else {
      Cmin = -cmax ;
    }

    if (is_void(cmin)) {
      Cmax = -min(img) ;
    } else {
      Cmax = -cmin ;
    }
    pli, -img, -img_size/2, -img_size/2, img_size/2, img_size/2, cmin=Cmin, cmax=Cmax ;

  } else {
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
  }


  // Plot beam
  if (!is_void(a)) {
    phi = span(0.,2*pi,360);
    x0 = -img_size/2 + max(nx / 10 * pixel_scale, 0.6*max(a,b)) ;
    y0 = x0 ;

    x= - a/2* cos(phi) * sin(theta) + b/2* sin(phi) * cos(theta)  ;
    y= + b/2* sin(phi) * sin(theta) + a/2* cos(phi) * cos(theta);

    for (l=1 ; l<=100 ; l++) {
      plg, y0 + 0.01*l*y, x0+ 0.01*l*x, width=1, color="white" ;
    }
  }
  if (!is_void(title)) xytitles, "RA offset (arcsec)", "DEC offset (sercsec)", [0.01,0.02];

}


//**************************************************

func plot_config(filename,color=, lambda=,clear=) {
/* DOCUMENT plot_config(filename,color=, lambda=)
     lambda must be given in mm
     baseline are then provided in klambda
   SEE ALSO:
 */
  if (is_void(color)) color="black" ;
  if (is_void(clear)) clear=1 ;

  system, "grep -v '#' "+filename+" > grep.tmp" ;
  config = OpenASCII("grep.tmp",prompt=0,noheader=1) ;
  system, "rm -rf grep.tmp"

  // pour les config temporaires cycle 1
  //ou = where(config.X1 > 0) ;
  //config = config(ou) ;

  x = config.X1 ; y = config.X2 ;

  x = x - avg(x) ;
  y = y - avg(y) ;

  if (clear) {
    window, 21, style="boxed.gs";
    lxy, 0, 0 ;
  } else {
    window, 21 ;
  }
  plp, y, x, color=color ;
  limits;  relimits ;
  l = limits() ; L = max(abs(l(1:4))) ;
  limits, -L, L , -L, L ;

  max_baseline = 0 ;
  min_baseline = 1e9 ;
  B = [] ;
  for (i=1 ; i<=numberof(x) ; i++) {
    baselines = abs(x - x(i), y-y(i)) ;
    ou = where(baselines > 1e-3) ;
    baselines = baselines(ou) ;
    grow, B, baselines ;

    max_baseline = max(max_baseline,max(baselines)) ;
    min_baseline = min(min_baseline,min(baselines)) ;

  }

  window, 22 ;
  x = B(sort(B)) ;
  y = span(0,1,numberof(B)) ;
  plg, y, x, color=color ;

  if (!is_void(lambda)) {
    write, "beam size (arsec)", lambda*1e-3/max_baseline / SI.as;
    write, "Max recoverable size", 0.6 * lambda*1e-3/min_baseline / SI.as ;
  }

  if (is_void(lambda)) {
    write, "min. baseline = ", min_baseline, "m" ;
    write, "max. baseline = ", max_baseline, "m" ;
    return B ;
  } else {
    write, "min. baseline = ", min_baseline/lambda, "klambda" ;
    write, "max. baseline = ", max_baseline/lambda, "klambda" ;
    return B/lambda ;
  }



}

write, "ALMA simdata interface loaded" ;
