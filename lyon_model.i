_lm_this_file    = "~/yorick/init/lyon_model.i";
_lm_default_file_format  = "lte%2d-%1.1f-0.0.AMES-%s";
_lm_default_file_format  = "lte%2d-%1.1f-0.0.%s.7.gz";
_atm_default_file_format  = "lte%d-%1.1f.%s.fits.gz";
_kurucz_default_file_format  = "%s%d-%1.1f.fits.gz";

//_lm_default_type         = "NextGen";
_lm_default_type         = "AMES-dusty";
//_lm_default_dir          = "~/tera/SPITZER/Model/";
_lm_default_dir = "~/mcfost/utils/Stellar_Spectra//";
if (is_void(_lm_size)) _lm_size = 100000;

struct Lyon_Model {
  float lamb(_lm_size);
  float flux(_lm_size);
  float bb(_lm_size);
  int T;
  float g;
  string type;
  string file;
}

// ===================MODELS TRACK ============
dusty_file     = "~/tera/Photometrie/Model/DUSTY00_models";
dusty_header   = ["mass","T","L","g","R","li","v","r","i","j","h","k","l","m"]

dusty_megacam_file     = "~/tera/Photometrie/Model/Megacam/DUSTY00_megacam";
dusty_megacam_header   = ["mass","T","L","g","U","G","i","r","z"];

dusty_CFHT12k_file     = "~/tera/Photometrie/Model/CFHT12K/DUSTY00_cfht12k";
dusty_CFHT12k_header   = ["mass","T","L","g","U","i", "G","z"];

bcah98_file    = "~/tera/Photometrie/Model/BCAH98_iso.1";
bcah98_header  = ["mass","T","g","L","v","r","i","j","h","k","l","m"];

cond_file     = "~/tera/Photometrie/Model/COND_models";
cond_header   = ["mass","T","L","R","g","v","r","i","z","j","h","k","l","m"]


megacam_dusty_file   = "~/tera/Photometrie/Model/Megacam/DUSTY00_megacam";
megacam_dusty_header = ["mass", "T","L", "g", "u", "G", "i", "r", "z"];

megacam_cond_file   = "~/tera/Photometrie/Model/Megacam/BCAH98_megacam";
megacam_cond_header = megacam_dusty_header;
//===== Default files, headers and structure for track model===========
_lm_default_header = dusty_header;
_lm_default_file  =  dusty_file;
struct LM_TRACK {
  float  mass , T, L , g , R , li ,u, v , r , i ,z , j , h ,k ,l ,m, age;
}
struct LM_ISOCHRONE {
  float age;
  pointer track;
}
//====================================================================

func _lm_open (file, lamb0=) {
  prompt = (_PROMPT ? _PROMPT : 1);
  if (file==string() || file=="") {
    if (prompt) write, "warning File name NULL";
    return Lyon_Model(type="NOFILE");
  }
  if (!open(file, "r", 1)) {
    if (prompt) write, "warning the file "+file+" does not exist";
    return Lyon_Model( type="NOFILE");
  }


  file
  system, "rm -rf /tmp/mod_lyon ; zcat "+file+" > /tmp/mod_lyon";
  tmp = OpenASCII("/tmp/mod_lyon",valtype="double",prompt=0);
  order = sort(tmp.X1);
  lambda = tmp.X1(order)/10000.0;
  ln10 = log(10);
  mod = 10^tmp.X2(order);
  bb = 10^tmp.X3(order);
  tmp=[];

  //f = open(file);
  //header = rdline(f); size = int();
  //sread, header, size, format="! %d";
  size = dimsof(lambda)(2);
  if (size!=_lm_size) {
    _lm_size = size;
    include , _lm_this_file, 1;
  }

  if (size!=_lm_size & is_void(lamb0)) {
    error, "Size of file is "+pr1(size)+" size of structure is "+pr1(_lm_size)+", you probabely try to create a array of model with different dimension";
  }

  T = int();
  g = float();
  file_name = file;
  while (((file_name=strtok(file_name, "/")))(2)) {
    file_name = file_name(2);
  }
  file_name = file_name(1);
  sread, file_name, T, g, format="lte%d-%f-";

  if (size==0) {
    if (prompt) write, "Warning the file "+file+" is probably void";
    return Lyon_Model(T=T*100, g=g, file=file, type="VOID");
  }

  if (strgrep("dusty", file)(0)>-1) type ="dusty";
  if (strgrep("cond", file)(0)>-1) type ="cond";
  if (strgrep("NextGen", file)(0)>-1) type ="NextGen";


  if (!is_void(lamb0)) return Lyon_Model(
                                         lamb=lamb0 ,
                                         flux=interp(mod, lambda, lamb0),
                                         bb= interp(bb, lambda , lamb0),
                                         T=T*100, g=g, file=file, type=type
                                         );

  return Lyon_Model(lamb=lambda, flux=mod, bb=bb, T=T*100, g=g, file=file);
}



func _atm_open(file, lamb0=) {
  // Lecture spectre
  tmp=cfitsRead(file);
  lambda = tmp(,1);
  mod = tmp(,2);
  bb = tmp(,3);
  tm=[];

    // taille
  if (is_void(lamb0)) {
    size = dimsof(lambda)(2);
  } else {
    size = dimsof(lamb0)(2);
  }

  if (size!=_lm_size) {
    _lm_size = size;
    include , _lm_this_file, 1;
  }

  // Lecture Temperature et g
  T = int();
  g = float();
  file_name = file;
  while (((file_name=strtok(file_name, "/")))(2)) {
    file_name = file_name(2);
  }
  file_name = file_name(1);

  // Lecture type
  if (strgrep("dusty", file)(0)>-1) type ="dusty";
  if (strgrep("cond", file)(0)>-1) type ="cond";
  if (strgrep("NextGen", file)(0)>-1) type ="NextGen";
  if (strgrep("Kurucz", file)(0)>-1) type ="Kurucz";

  if (type=="Kurucz") {
    sread, file_name, T, g, format="Kurucz%d-%f-";
  } else {
    sread, file_name, T, g, format="lte%d-%f-";
  }


  if (!is_void(lamb0)) {  // Ca merde, faut faire un truc mieux que ca
    flux = cn = array(double, numberof(lamb0));

    for (i=1 ; i<=numberof(lamb0) ; i++) {
      ou = where(abs(lambda-lamb0(i)) < 0.01*lamb0(i));
      if (numberof(ou) < 2) {
        flux(i)=exp(interp(log(mod+1e-300), log(lambda), log(lamb0(i))));
        cn(i)=interp(bb, lambda, lamb0(i));
      } else {
        M = mod(ou) ; // F_lambda
        L = lambda(ou) ;
        flux(i) = sum(M(zcen) * L(dif)) / sum(L(dif));
        //cn(i) = sum(bb(zcen) * L(dif)) / sum(L(dif));
      }
    }


    return Lyon_Model(lamb=lamb0 ,
                      flux=flux,
                      bb=cn,
                      T=T, g=g, file=file, type=type
                      );
  }

  return Lyon_Model(lamb=lambda, flux=mod, bb=bb, T=T, g=g, file=file);
}


func atm_open(files, lamb0=) {
  /* DOCUMENT model = atm_open(files, lamb0=lamb)
   open Model spectra in MCFOST atm spectrum format with a liste of files.
   en F_lambda
   SEE ALSO:
 */
  N = numberof(files);
  if (files(1)==string() && is_void(lamb0)) error, "The first file name is NULL";

  //definition de la taille.
  if (is_void(lamb0)) {
    tmp_model =  _atm_open (files(1));
  } else {
    tmp_model =  _atm_open (files(1), lamb0=lamb0);
  }


  out_models = array (Lyon_Model, dimsof(files));
  out_models(1) =  tmp_model;


  if (is_void(lamb0)) {
    for (i=2; i<=N ; i++) {
      out_models(i) =  _atm_open (files(i));
    }
  } else {
    for (i=2; i<=N ; i++) {
      out_models(i) =  _atm_open (files(i), lamb0=lamb0);
    }
  }

  return out_models;
}


func lm_open (files, lamb0=) {
/* DOCUMENT model = lm_open(files, lamb0=lamb)
   open Lyon Model spectra with a liste of files.

   SEE ALSO:
 */
  extern _lm_size;
  N = numberof(files);
  if (files(1)==string() && is_void(lamb0)) error, "The first file name is NULL";

  //definition de la taille.
  if (is_void(lamb0)) {
    tmp_model =  _lm_open (files(1));
    lamb0 = tmp_model.lamb;
  }
  else size = numberof(lamb0);

  out_models = array (Lyon_Model, dimsof(files));

  if (is_void(lamb0)) {
    out_models(1) =  tmp_model;
  }
  else out_models(1) =  _lm_open (files(1), lamb0=lamb0);
  for (i=2; i<=N ; i++) {
    out_models(i) =  _lm_open (files(i), lamb0=lamb0);
  }
  return out_models;
}


func atm_get_files(T, g, type, format=, dir=, prompt=) {
  if (is_void(prompt)) prompt = _PROMPT;
  if (is_void(type)) type= _lm_default_type;
  if (is_void(dir)) dir= _lm_default_dir;
  size_T = numberof(T) ; size_g = numberof(g) ; size_type= numberof(type);
  T =       T(, -:1:size_g, -:1:size_type);
  g =       g(-:1:size_T, , -:1:size_type);
  type = type(-:1:size_T, -:1:size_g,);

  T = int(T);
  N = numberof(T);
  files = array(string, dimsof(T));

  infile = lsdir(dir);

  for (i=1; i<=N ; i++) {
    files(i) = _atm_get_file(T(i), g(i), type(i), format=format, infile=infile, prompt=prompt);
  }
  return files;
}


func lm_get_files(T, g, type, format=, dir=, prompt=) {
  if (is_void(prompt)) prompt = _PROMPT;
  if (is_void(type)) type= _lm_default_type;
  if (is_void(dir)) dir= _lm_default_dir;
  size_T = numberof(T) ; size_g = numberof(g) ; size_type= numberof(type);
  T =       T(, -:1:size_g, -:1:size_type);
  g =       g(-:1:size_T, , -:1:size_type);
  type = type(-:1:size_T, -:1:size_g,);

  T = int(T/100) // les temperature sont en K/100
  N = numberof(T);
  files = array(string, dimsof(T));

  infile = lsdir(dir);

  for (i=1; i<=N ; i++) {
    files(i) = _lm_get_file(T(i), g(i), type(i), format=format, infile=infile, prompt=prompt);
  }
  return files;
}


func atm_get_list_files( T, g, type, format=, dir=, prompt=) {
  if (is_void(prompt)) prompt = _PROMPT;
  if (is_void(type)) type= _lm_default_type;
  if (is_void(dir)) dir= _lm_default_dir;
  size_T = numberof(T) ; size_g = numberof(g) ; size_type= numberof(type);
  if (size_T!=size_g || size_g!=size_type) error , "T g and type must have the same dimension";

  files = array(string,size_T);
  infile = lsdir(dir);
  for (i=1;i<=size_T;i++) {
    files(i) = _atm_get_file(T(i), g(i), type(i), format=format, infile=infile, prompt=prompt);
  }
  return files;
}

func lm_get_list_files( T, g, type, format=, dir=, prompt=) {
  if (is_void(prompt)) prompt = _PROMPT;
  if (is_void(type)) type= _lm_default_type;
  if (is_void(dir)) dir= _lm_default_dir;
  size_T = numberof(T) ; size_g = numberof(g) ; size_type= numberof(type);
  if (size_T!=size_g || size_g!=size_type) error , "T g and type must have the same dimension";

  files = array(string,size_T);
  infile = lsdir(dir);
  for (i=1;i<=size_T;i++) {
    files(i) = _lm_get_file(T(i)/100, g(i), type(i), format=format, infile=infile, prompt=prompt);
  }
  return files;
}

func lm_get_best_file (T , age) {

  _lm_cond_dusty_limit_temperature = 1500;
  type0 = ["cond","dusty"];

  type = type0((T>_lm_cond_dusty_limit_temperature)+1);

  dusty = track_interp (lm_open_track(dusty_file, dusty_header), "T",T);
  cond  = track_interp (lm_open_track(cond_file , cond_header ), "T",T);
  cond.age = 10^cond.age;

  ou = where(T>_lm_cond_dusty_limit_temperature);
  if (numberof(ou)) {
    dusty = track_interp (lm_open_track(dusty_file, dusty_header), "T",T(ou));
    Dage = abs(dusty.age - age(-,));
    prox_track = dusty.age(Dage(mnx,));
  }
}

func _atm_get_file (T, g, type, format=, infile=, prompt=) {
  if (is_void(format)) format = _atm_default_file_format;
  if (type=="Kurucz") {
    format = _kurucz_default_file_format ;
    file =  swrite (type, T, g,  format=format);
  } else {
    file =  swrite (T , g , type ,  format=format);
  }
  file_size = strlen(file);
  ou = where( strpart(infile, 1:file_size)==file);
  if (numberof(ou)) {
    return dir+"/"+infile(ou(1));
  }
  if (prompt) write, "WARNING : I can't find the file "+dir+"/"+file+"*";
  return string();
}

func _lm_get_file (T, g, type, format=, infile=, prompt=) {
  if (is_void(format)) format = _lm_default_file_format;
  file =  swrite (T , g , type ,  format=format);
  file_size = strlen(file);
  ou = where( strpart(infile, 1:file_size)==file);
  if (numberof(ou)) {
    return dir+"/"+infile(ou(1));
  }
  if (prompt) write, "WARNING : I can't find the file "+dir+"/"+file+"*";
  return string();
}


func plot_lm(T,g,type,r_star,dist,av=) {
/* DOCUMENT  plot_lm(T,g,type,r_star,dist,av=)

   SEE ALSO:
 */
  T=double(T);
  g=double(g);
  r_star=double(r_star);
  dist=double(dist);

  sigma_stephan = 5.6697e-8;
  pc_to_AU =  206264.81;
  Rsun_to_AU = 0.00466666666;

  lm = _lm_open( lm_get_files(T,g,type)(1));

  if (is_void(av)) av=0.;

  L_bol_mod = sum(lm.lamb(dif)*lm.flux(zcen));

  L_bol = sigma_stephan * (r_star/(dist*pc_to_AU/Rsun_to_AU))^2 * T^4;

  if (av > 0.) {
    extinction = OpenASCII("~/yorick/SED/extinction_law.dat",prompt=0);
    extinction.lamb = extinction.lamb;

    ext_V = abs(extinction.lamb - 5.47e-01)(mnx);
    XX =  extinction.kpa /  extinction.kpa(ext_V);
    XXm = interp( XX , extinction.lamb ,  lm.lamb);

    tau_V=0.4*log(10)*av;
    correct_Av = exp(-tau_V * XXm);
  } else {
    correct_Av=1.0;
  }


  lm.flux = lm.flux * L_bol/L_bol_mod;
  lm.bb = lm.bb * L_bol/L_bol_mod;
  plg, lm.flux*lm.lamb*correct_Av, lm.lamb

}


func plot_atm(T,g,star_type,r_star,dist,av=,rv=,color=,type=,width=,lamb0=,bb=,sm=,write=,extinction=) {
/* DOCUMENT  plot_atm(T,g,star_type,r_star,dist,av=,color=,type=,width=,sm=)
     r_star en rayon solaire
     dist en pc
     rv = 3.1, 4.0 ou 5.5
    SEE ALSO:
 */

  if (is_void(color)) color="black";
  if (is_void(type)) type=1;
  if (is_void(width)) width=1;
  if (is_void(bb)) bb=0;
  if (is_void(sm)) sm=0;
  if (is_void(rv)) rv=3.1;
  if (is_void(write)) write=1;

  srv = swrite(rv,format="%3.1f") ;

  T=double(T);
  g=double(g);
  r_star=double(r_star);
  dist=double(dist);

  sigma_stephan = 5.6697e-8;
  pc_to_AU =  206264.81;
  Rsun_to_AU = 0.00466666666;

  if (is_void(lamb0)) {
    atm = atm_open( atm_get_files(T,g,star_type)(1));
  } else {
    atm = atm_open( atm_get_files(T,g,star_type)(1), lamb0=lamb0);
  }


  if (is_void(av)) av=0.;

  if (av > 0.) {
    OPEN_ASCII_TMP = [] ;
    eval_code = [] ;
    if (is_void(extinction)) {
      extinction = OpenASCII("~/yorick/SED/kext_albedo_WD_MW_"+srv+"_D03.all",prompt=1) ;
    }
    extinction.lamb = extinction.lamb;
    ext_V = abs(extinction.lamb - 5.47e-01)(mnx);

    kext = extinction.kpa / (1-extinction.albedo) ;
    XX =  kext /  kext(ext_V);
    XXm = interp( XX , extinction.lamb ,  atm.lamb);

    tau_V=0.4*log(10)*av;
    correct_Av = exp(-tau_V * XXm);
  } else {
    correct_Av=1.0;
  }


  L = sum(atm.flux(zcen)*atm.lamb(dif)) * 4*pi * SI.pc^2 * (r_star/dist)^2 ;
  //write, "Luminosity = ", L / SI.Lsun,"Lsun" ;




  atm.flux = atm.flux * (r_star/dist)^2 * correct_Av;
  //plg, smooth(atm.flux*atm.lamb(,1),sm) , atm.lamb(,1), color=color, type=type, width=width;
  plg, atm.flux*atm.lamb(,1) , atm.lamb(,1), color=color, type=type, width=width;

  if (bb) {
    atm.bb = atm.bb * (r_star/dist)^2 * correct_Av;
    plg, atm.bb*atm.lamb(,1) , atm.lamb(,1), type=type, width=width;
  }

}



func _atm(T,g,star_type,r_star,dist,av=,lamb0=) {
/* DOCUMENT  _atm(T,g,star_type,r_star,dist,av=,lamb0=)
     r_star en rayon solaire
     dist en pc
     return les flux et lambda plotte par plot_atm
   SEE ALSO:
 */

  if (is_void(color)) color="black";
  if (is_void(type)) type=1;
  if (is_void(width)) width=1;

  T=double(T);
  g=double(g);
  r_star=double(r_star);
  dist=double(dist);

  sigma_stephan = 5.6697e-8;
  pc_to_AU =  206264.81;
  Rsun_to_AU = 0.00466666666;

  if (is_void(lamb0)) {
    atmo = atm_open( atm_get_files(T,g,star_type)(1));
  } else {
    atmo = atm_open( atm_get_files(T,g,star_type)(1), lamb0=lamb0);
  }


  if (is_void(av)) av=0.;

  if (av > 0.) {
    extinction = OpenASCII("~/yorick/SED/extinction_law.dat",prompt=0);
    extinction.lamb = extinction.lamb;

    ext_V = abs(extinction.lamb - 5.47e-01)(mnx);
    XX =  extinction.kpa /  extinction.kpa(ext_V);
    XXm = interp( XX , extinction.lamb ,  atmo.lamb);

    tau_V=0.4*log(10)*av;
    correct_Av = exp(-tau_V * XXm);
  } else {
    correct_Av=1.0;
    }


  atmo.flux = atmo.flux * (r_star/dist)^2 * correct_Av;
  atmo.bb = atmo.bb * (r_star/dist)^2 * correct_Av ;

  //  plg, atmo.flux*atmo.lamb(,1) , atmo.lamb(,1), color=color, type=type, width=width;

  atmo.flux = atmo.flux*atmo.lamb(,1);
  atmo.bb = atmo.flux*atmo.lamb(,1);
  atmo.lamb = atmo.lamb(,1);

  return atmo

}



func atm(T,g,star_type,r_star,dist,av=,lamb0=) {
/* DOCUMENT  _atm(T,g,star_type,r_star,dist,av=,color=,type=,width=,sm=)
     r_star en rayon solaire
     dist en pc
     retourne les flux et lambda plotte par plot_atm
     interpole entre les temperatures
   SEE ALSO:
 */

  T=double(T);
  g=double(g);
  r_star=double(r_star);
  dist=double(dist);

  // Selection des modeles du bon type
  liste = lsdir(_lm_default_dir);
  ou = where((strgrep(star_type, liste))(2,) > 0) ;
  liste = liste(ou) ;

  n = dimsof(liste)(2);


  tmp_T = int();
  tmp_g = float();

  Temp = array(int,n) ;
  G = array(float, n) ;

  // On lit les temperatures et les logg
  for (i=1 ; i<=n ; i++) {
    sread, liste(i), tmp_T, tmp_g, format="lte%d-%f";
    Temp(i) = tmp_T ;
    G(i) = tmp_g ;
  }

  // Selection des modeles au bon g
  ou = where(abs(g - G) < 1e-3);
  Temp = Temp(ou);

  ou = where(Temp <= T) ;
  if (!is_array(ou)) {
    atmo = _atm(T,g,star_type,r_star,dist,av=av,lamb0=);
    return atmo ;
  }
  Tinf = max(Temp(ou));

  ou = where(Temp > T) ;
  if (!is_array(ou)) {
    atmo = _atm(T,g,star_type,r_star,dist,av=av,lamb0=);
    return atmo ;
  }
  Tsup = min(Temp(ou));

  frac = (double(Tsup)^4-T^4)/(double(Tsup)^4-double(Tinf)^4);

  atm_inf = _atm(Tinf,g,star_type,r_star,dist,av=av,lamb0=lamb0);
  atm_sup = _atm(Tsup,g,star_type,r_star,dist,av=av,lamb0=lamb0);

  atmo = atm_inf ;
  atmo.flux = atm_inf.flux * frac + atm_sup.flux * (1.0-frac);
  atmo.bb = atm_inf.bb * frac + atm_sup.bb * (1.0-frac);

  return atmo ;

}


func fit_atm(T,g,type,r_star,dist,av,obs,proba=) {

  sigma_stephan = 5.6697e-8;
  pc_to_AU =  206264.81;
  Rsun_to_AU = 0.00466666666;

  // Lecture modeles
  mod = atm_open(atm_get_files(T,g,type),lamb0=obs.lamb);
  //  plg, mod.flux(,1)* (r_star/dist)^2*mod.lamb(,1) , mod.lamb(,1), color="green"

  mod_flux = (mod.flux*mod.lamb)(..,-) * ((r_star/dist)^2)(-,-,-,-,);

  chi2 = array(double, numberof(T), numberof(g), numberof(type), numberof(r_star), numberof(av))


  // Lecture extinction
  extinction = OpenASCII("~/yorick/SED/extinction_law.dat",prompt=0);
  extinction.lamb = extinction.lamb;
  ext_V = abs(extinction.lamb - 5.47e-01)(mnx);
  XX =  extinction.kpa /  extinction.kpa(ext_V);
  XXm = interp( XX , extinction.lamb ,  obs.lamb);
  correct_Av= array(double, numberof(obs.lamb), numberof(av));


  for (i=1; i<=numberof(av) ; i++) {
    // On multiplie les modeles par Av
    if (av(i) > 0.) {
      tau_V=0.4*log(10)*av(i);
      correct_Av(,i) = exp(-tau_V * XXm);
    } else {
      correct_Av(,i)=1.0;
    }
    chi2(..,i) = (((mod_flux*correct_Av(,i) - obs.flux)/(0.1*obs.flux))^2)(sum,..)
  }

  best = indexof(chi2,(chi2(*))(mnx));


  if (!is_void(proba)) {
    write, "~/tmp/proba_"+proba+".txt" ;

    f = open("~/tmp/proba_"+proba+".txt","w") ;

    //chi2 = array(double, numberof(T), numberof(g), numberof(type), numberof(r_star), numberof(av))
    N_para=0 ;
    if (numberof(T) > 1) N_para++ ;
    if (numberof(g) > 1) N_para++ ;
    if (numberof(Type) > 1) N_para++ ;
    if (numberof(r_star) > 1) N_para++ ;
    if (numberof(av) > 1) N_para++ ;

    P = exp(-0.5 * chi2 / (numberof(obs) - N_para)) ;
    P = P/sum(P) ;

    write, f, "T proba" ;
    write, f, T, P(,sum,sum,sum,sum) ;
    write, f, "";

    write, f, "g proba" ;
    write, f, g, P(sum,,sum,sum,sum) ;
    write, f, "";

    write, f, "Type proba" ;
    write, f, type, P(sum,sum,,sum,sum) ;
    write, f, "";

    write, f, "r_star proba" ;
    write, f, r_star, P(sum,sum,sum,,sum) ;
    write, f, "";

    write, f, "Av proba" ;
    write, f, av, P(sum,sum,sum,sum,) ;

    close, f;
  }
  /*plp, obs.flux, obs.lamb;
  plg, mod_flux(,best(1),best(2),best(3),best(4)) * correct_Av(,best(5)), obs.lamb, color="red";

  print, "chi2=", min(chi2);
  print, "T=", T(best(1)), "g=", g(best(2)), type(best(3)), "r_star=", r_star(best(4)), "Av=", av(best(5));
  */
  return best
}



func _lm_to_mcfost(filenames) {
/* DOCUMENT _lm_to_mcfost(filenames)
     Converti les modeles de France Allard en format pour MCFOST
     Les flux sont normalises pour R=1 rayon solaire vu de 1 pc
     Unites : W.m-2 /mum
   SEE ALSO: lm_to_mcfost
 */
  sigma_stephan = 5.6697e-8;
  pc_to_AU =  206264.81;
  Rsun_to_AU = 0.00466666666;


  if (dimsof(filenames)(1) ==0) filenames = [filenames];
  lm = lm_open(filenames);
  n = dimsof(lm)(2);

  cst = sigma_stephan * (1.0/(pc_to_AU/Rsun_to_AU))^2;

  for (i=1 ; i<=n ; i++) {
    tmp = [];

    // Normalisation des modeles a 1 rayon solaire vu de 1 pc
    L_bol_mod = sum(lm(i).lamb(dif)*lm(i).flux(zcen));
    L_bol = cst * double(lm(i).T)^4;
    correct=L_bol/L_bol_mod;
lm(i).type
    grow, tmp, [lm(i).lamb(*)], [lm(i).flux(*)*correct], [lm(i).bb(*)*correct];
    file =  swrite (lm(i).T , lm(i).g , lm(i).type ,  format="~/mcfost_utils/lte%d-%1.1f.%s.fits.gz");
    cfitsWrite, file, tmp;
  }
}




func lm_to_mcfost(T, g, type) {
/* DOCUMENT lm_to_mcfost(T, g, type)
     Converti les modeles de France Allard en format pour MCFOST
     Les flux sont normalises pour R=1 rayon solaire vu de 1 pc
     Unites : W.m-2
   SEE ALSO: _lm_to_mcfost
 */
  _lm_to_mcfost, lm_get_files(T,g,type)(1);
}


func lm_open_track (file, header) {
  if (is_void(header)) header=_lm_default_header;
  if (is_void(file))   file=_lm_default_file;

  f = open(file, "r");
  isochrone = [];
  newtrack :
  while ( (line=rdline(f)) && !regmatch( "-------------------", line))
    {pline=line;}
  if (line) {
    age = float();
    sread, strtok(pline, "=")(2) , age, format="%f";
    rdline, f;
    rdline, f;
    grow ,isochrone, LM_ISOCHRONE(  track=&_lm_read_track(f,header,age), age=age );
    goto newtrack;
  }
  return isochrone;
}

func _lm_read_track (f, header, age) {

  out_track = []; ptab =[];
  N = numberof(header);
  while ( !regmatch("-------------------",(line=rdline(f))) ) {
    tmp_tab = array(float, N);
    //ptab    = array(pointer, N);
    sread,line , tmp_tab;

    grow, ptab, [tmp_tab];
  }
  st_name = strtrim(strtok(strtok(print(LM_TRACK)(2:-1), " ")(2,),";")(1,));
  test = header(,-:1:numberof(st_name)) == st_name(-:1:numberof(header),);

  extern __LM_TMP_VAR__;
  __LM_TMP_VAR__ = array(LM_TRACK, numberof(ptab(1,)));
  eval_code = " ";
  ou = where(test(,sum))
  for (i=1;i<=numberof(ou); i++) {

    eval_code = eval_code + "__LM_TMP_VAR__."+header(ou(i)) +"="+strjoin(print(ptab(ou(i),)))+";";

  }
  eval, eval_code, debug=0;

  out_track = __LM_TMP_VAR__;
  out_track.age = age;
  return out_track;
}


func _save_lm_read_track (f, age) {
  mass = T = L = g = R = li = v = r = i = j = h =k =l = m = float();
  out_track = [];
  while ( !regmatch("-------------------",(line=rdline(f))) ) {
    sread,line ,   mass , T, L , g , R , li , v , r , i , j , h ,k ,l ,m;
    _ , out_track,  LM_TRACK(mass=mass, T=T,  L =L,  g =g,  R =R, li=li, v=v,
                             r=r,  i=i ,j=j, h=h, k=k , l=l, m=m, age=age);
  }
  return out_track;
}


func track_interp(isochrone, elem_name, xp) {
  N1 = numberof(isochrone);
  N2 =  numberof(xp);
  out_track = array(LM_TRACK, N1, N2);
  for (i=1; i<=N1; i++) {
    out_track(i,) = struct_interp ( *isochrone(i).track , elem_name, xp );
  }
  return out_track;

}
