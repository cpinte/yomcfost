require, "init/mcfost_struct.i" ;

// REINIT FUNCTION;
func reinit_all {
  extern nrad,nz,nthet,nlamb,nphi;
  extern vect_n_phi,vect_n_inc,vect_n_sed ;
  nrad = nz  = nthet = nlamb = nphi = [];
  vect_n_phi =  vect_n_inc = vect_n_sed = [];
  include , struct_file , 1;
}


// fonctions de lecture
func _read_params(file) {
  extern lambda_min, lambda_max, n_lines_ray_tracing ;

  version = val1 = val2 = val3 = val4 = val5 = val6 = val7 = val8 = val9 =float();
  com1 = com2 = string();
  string1 = string2 = string3 = string4 = string();

  f = open(file,"r");
  param = McfostParams();
  read, f, val1, com1; param.simu.version = val1;

  rdline, f;

  lambda_min = 0.1; lambda_max = 3000. ;
  param.simu.nmol = 1 ;
  param.zones(:).vert_exp = 2 ;

  if ((abs(param.simu.version -2.21) < 1.0e-6) | (abs(param.simu.version -3.0) < 1.0e-6)) {
    rdline, f;

    //photon;
    read, f, val1 ; param.phot.n2 = val1;
    read, f, val1 ; param.phot.nlambda = val1;
    read, f, val1 ; param.phot.nimage = val1;

    rdline, f;
    rdline, f;

    //Wavelength
    read, f, val1, val2, val3 ;
    param.wave.nlamb1 = val1; param.wave.lambda_min = val2; param.wave.lambda_max = val3; lambda_min=val2 ; lambda_max=val3;
    read, f, string1, string2, string3 ; param.simu.ltemp = string1; param.simu.lsed = string2; param.simu.lsed_complete = string3;
    read, f, string1 ; param.wave.file = string1;
    param.simu.lem_disk_image = "T";
    read, f, string1, string2 ; param.simu.lsepar = string1; param.simu.lsepar_pola = string2;
    param.map.n_type=1;
    if (param.simu.lsepar=="T") param.map.n_type += 4;
    if (param.simu.lsepar_pola=="T") param.map.n_type += 3;

    rdline, f;
    rdline, f;

    //Grid
    read, f, val1; param.grid.geometry = val1;
    read, f, val1, val2, val3, val4 ;
    param.grid.nrad=val1; param.grid.nz=val2; param.grid.n_az=val3 ; param.grid.nrad_in=val4;

    rdline, f;
    rdline, f;

    //map;
    read, f, val1, val2, val3 ;
    param.map.nx = val1; param.map.ny=val2; param.map.map_size=val3; param.map.zoom = 1.0 ;
    param.map.size_neb = 0.5 * param.map.map_size ;
    read, f, val1, val2, val3, string1;
    param.map.RT_imin = val1 ; param.map.RT_imax = val2 ; param.map.RT_ntheta = val3 ; param.map.lRT_i_centered = string1 ;
    read, f, val1, val2, val3 ;
    param.map.RT_az_min = val1 ; param.map.RT_az_max = val2 ; param.map.RT_naz = val3 ;
    read, f, val1 ; param.map.dist = val1;
    read, f, val1 ; param.map.PA = val1;

    rdline, f;
    rdline, f;

    // Scattering method
    read, f, val1 ; param.simu.scatt_method = val1;
    read, f, val1 ; param.simu.aniso_method = val1;

    rdline, f;
    rdline, f;

    // Symetries
    read, f, string1; param.simu.lsym_ima = string1 ;
    read, f, string1; param.simu.lsym_centrale = string1 ;
    read, f, string1; param.simu.lsym_axiale = string1 ;

    rdline, f;
    rdline, f;

    // Dust global properties
    read, f, val1, val2, val3, com2; param.dust.strat_type = val1 ;  param.dust.strat = val2; param.dust.a_strat = val3;
    if (val1 > 0) {
      param.dust.lstrat = "T" ;
    } else {
      param.dust.lstrat = "F" ;
    }
    read, f, string1, val1; param.dust.lmigration = string1 ;
    read, f, string1, val1; param.simu.ldust_sublimation = string1 ; param.simu.correct_Rsub = val1 ;
    read, f, string1, val1; param.simu.lhydro_eq = string1 ;
    read, f, string1 , val1 ; param.simu.lchauff_int = string1 ; param.simu.viscosity = val1;

    read, f;
    read, f;

    //Nbre de zones
    read,f, val1 ; param.simu.n_zones = val1;

    // On arrete ici si on n'a pas le bon nombre de zones et on recharche la structure
    if (n_zones != param.simu.n_zones) return param

    read, f;
    read, f;

    //Disk
    for (i=1 ; i <=  param.simu.n_zones ; i++) {
      read, f, val1 ; param.zones(i).geometry = val1 ;
      read, f, val1, val2; param.zones(i).mass = val1; param.zones(i).gas_to_dust = val2;
      read, f, val1, val2, val3; param.zones(i).Ho = val1; param.zones(i).Ro = val2; param.zones(i).vert_exp = val3 ;
      read, f, val1, val2, val3, val4;  param.zones(i).rin = val1; param.zones(i).edge = val2; param.zones(i).rout = val3;  param.zones(i).Rc = val4;
      read, f, val1; param.zones(i).beta = val1;
      read, f, val1, val2; param.zones(i).dens = val1; param.zones(i).gamma_exp = val2;
    }
    rdline, f;
    rdline, f;

    //Grain
    param.simu.n_especes_max = 0 ;
    param.simu.n_components_max = 0 ;

    for (i=1 ; i <=  param.simu.n_zones ; i++) {
      read, f, val1 ; param.zones(i).n_especes=val1;
      if (val1 > param.simu.n_especes_max) param.simu.n_especes_max = val1 ;

      for (j=1 ; j<= param.zones(i).n_especes ; j++) {
        read, f, string1, val1, val2, val3, val4, val5 ;
        if (param.zones(i).n_especes <= n_especes) {
          param.dust_pop(i,j).type = string1 ;
          param.dust_pop(i,j).n_components = val1 ; if (val1 > param.simu.n_components_max) param.simu.n_components_max = val1 ;
          param.dust_pop(i,j).mixing_rule = val2 ;
          param.dust_pop(i,j).porosity = val3;
          param.dust_pop(i,j).mass_fraction = val4;
          param.dust_pop(i,j).dhs_maxf = val5;
        }

        if (n_components < param.simu.n_components_max) return param

        for (k=1 ; k<=param.dust_pop(i,j).n_components ; k++) {
          read, f, string1, val1 ;
          if (param.dust_pop(i,j).n_components <= n_components) {
            param.dust_pop(i,j).file(k) = string1;  param.dust_pop(i,j).component_volume_fraction(k) = val1 ;
          }
        }
        read, f, val1; param.dust_pop(i,j).type_grain = val1;
        read, f, val1, val2, val3, val4 ; param.dust_pop(i,j).amin=val1;  param.dust_pop(i,j).amax=val2;  param.dust_pop(i,j).aexp=val3; param.dust_pop(i,j).n_grains=val4;
      }
    }


    if (param.dust_pop(1,1).type_grain == 2) {
      param.gridT.n_T_grain = sum(param.dust_pop.n_grains);
    } else {
      param.gridT.n_T_grain = 1;
    }
    rdline, f;
    rdline, f;


    // mol
    read, f, string1, string2, string3, val1; param.mol.lpop=string1 ; param.mol.lpop_precise=string2 ; param.mol.lLTE=string3 ; param.mol.profile_width = val1;

    read, f, val1 ; param.mol.vturb=val1 ;
    read, f, val1 ; param.simu.nmol=val1 ;

    if (param.simu.nmol == nmol)  {
      for (i=1 ; i <=  param.simu.nmol ; i++) {
        read, f, string1, val1; param.mol(i).molecule_file=string1 ; param.mol(i).level_max = val1;
        read, f, val1, val2 ;  param.mol(i).vmax=val1 ;  param.mol(i).n_speed = val2;
        read, f, string1, val1, string2 ; param.mol(i).cst_abundance=string1 ; param.mol(i).abundance=val1 ; param.mol(i).abundance_file = string2;
        read, f, string1, val1; param.mol(i).ray_tracing=string1 ; param.mol(i).n_lines_ray_tracing=val1;
        if (param.mol(i).n_lines_ray_tracing == 1) {
          read, f, val1 ; param.mol(i).transition_numbers(1) = val1;
        } else if (param.mol(i).n_lines_ray_tracing == 2) {
          read, f, val1, val2 ; param.mol(i).transition_numbers(1:2) = [val1,val2];
        } else if (param.mol(i).n_lines_ray_tracing == 3) {
          read, f, val1, val2, val3 ; param.mol(i).transition_numbers(1:3) = [val1,val2,val3];
        } else if (param.mol(i).n_lines_ray_tracing == 4) {
          read, f, val1, val2, val3, val4 ; param.mol(i).transition_numbers(1:4) = [val1,val2,val3,val4];
        } else if (param.mol(i).n_lines_ray_tracing == 5) {
          read, f, val1, val2, val3, val4, val5 ; param.mol(i).transition_numbers(1:5) = [val1,val2,val3,val4,val5];
        }
      }
    }
    rdline, f;
    rdline, f;


    // Star
    read, f, val1 ; param.simu.n_stars = val1;
    if (param.simu.n_stars == n_stars) {
      for (i=1 ; i <= param.simu.n_stars ; i++) {
        read, f, val1, val2, val3, val4, val5, val6, string1, val7; param.stars(i).T = val1; param.stars(i).R = val2 ;  param.stars(i).Mass = val3;
        param.stars(i).x = val4 ; param.stars(i).y = val5 ; param.stars(i).z = val6 ;
        param.stars(i).is_bb = string1;
        read, f, param.stars(i).spectre ;
        read, f, val1, val2 ;
        param.stars(i).fUV = val1 ;
        param.stars(i).slope_fUV = val2 ;
      }
    }

  } else if (abs(param.simu.version -2.20) < 1.0e-6) {
    rdline, f;

    //photon;
    read, f, val1 ; param.phot.n2 = val1;
    read, f, val1 ; param.phot.nlambda = val1;
    read, f, val1 ; param.phot.nimage = val1;

    rdline, f;
    rdline, f;

    //Wavelength
    read, f, val1, val2, val3 ;
    param.wave.nlamb1 = val1; param.wave.lambda_min = val2; param.wave.lambda_max = val3; lambda_min=val2 ; lambda_max=val3;
    read, f, string1, string2, string3 ; param.simu.ltemp = string1; param.simu.lsed = string2; param.simu.lsed_complete = string3;
    read, f, string1 ; param.wave.file = string1;
    param.simu.lem_disk_image = "T";
    read, f, string1, string2 ; param.simu.lsepar = string1; param.simu.lsepar_pola = string2;
    param.map.n_type=1;
    if (param.simu.lsepar=="T") param.map.n_type += 4;
    if (param.simu.lsepar_pola=="T") param.map.n_type += 3;

    rdline, f;
    rdline, f;

    //Grid
    read, f, val1; param.grid.geometry = val1;
    read, f, val1, val2, val3, val4 ;
    param.grid.nrad=val1; param.grid.nz=val2; param.grid.n_az=val3 ; param.grid.nrad_in=val4;

    rdline, f;
    rdline, f;

    //map;
    read, f, val1, val2, val3 ;
    param.map.nx = val1; param.map.ny=val2; param.map.map_size=val3; param.map.zoom = 1.0 ;
    param.map.size_neb = 0.5 * param.map.map_size ;
    read, f, val1, val2, val3 ;
    param.map.nthet = val1;   param.map.nphi = val2  ; param.map.inc = val3;
    param.map.delta = 1 ; param.map.angle_interest = 75.; param.map.lonly_bin_interest="F"
    read, f, val1, val2, val3, string1;
    param.map.RT_imin = val1 ; param.map.RT_imax = val2 ; param.map.RT_ntheta = val3 ; param.map.lRT_i_centered = string1 ;
    read, f, val1, val2, val3 ;
    param.map.RT_az_min = val1 ; param.map.RT_az_max = val2 ; param.map.RT_naz = val3 ;
    read, f, val1 ; param.map.dist = val1;
    read, f, val1 ; param.map.PA = val1;

    rdline, f;
    rdline, f;

    // Scattering method
    read, f, val1 ; param.simu.scatt_method = val1;
    read, f, val1 ; param.simu.aniso_method = val1;

    rdline, f;
    rdline, f;

    // Symetries
    read, f, string1; param.simu.lsym_ima = string1 ;
    read, f, string1; param.simu.lsym_centrale = string1 ;
    read, f, string1; param.simu.lsym_axiale = string1 ;

    rdline, f;
    rdline, f;

    // Dust global properties
    read, f, val1, val2, val3, com2; param.dust.strat_type = val1 ;  param.dust.strat = val2; param.dust.a_strat = val3;
    if (val1 > 0) {
      param.dust.lstrat = "T" ;
    } else {
      param.dust.lstrat = "F" ;
    }
    read, f, string1, val1; param.dust.lmigration = string1 ;
    read, f, string1, val1; param.simu.ldust_sublimation = string1 ; param.simu.correct_Rsub = val1 ;
    read, f, string1, val1; param.simu.lhydro_eq = string1 ;
    read, f, string1 , val1 ; param.simu.lchauff_int = string1 ; param.simu.viscosity = val1;

    read, f;
    read, f;

    //Nbre de zones
    read,f, val1 ; param.simu.n_zones = val1;

    // On arrete ici si on n'a pas le bon nombre de zones et on recharche la structure
    if (n_zones != param.simu.n_zones) return param

    read, f;
    read, f;

    //Disk
    for (i=1 ; i <=  param.simu.n_zones ; i++) {
      read, f, val1 ; param.zones(i).geometry = val1 ;
      read, f, val1, val2; param.zones(i).mass = val1; param.zones(i).gas_to_dust = val2;
      read, f, val1, val2, val3; param.zones(i).Ho = val1; param.zones(i).Ro = val2; param.zones(i).vert_exp = val3 ;
      read, f, val1, val2, val3, val4;  param.zones(i).rin = val1; param.zones(i).edge = val2; param.zones(i).rout = val3;  param.zones(i).Rc = val4;
      read, f, val1; param.zones(i).beta = val1;
      read, f, val1, val2; param.zones(i).dens = val1; param.zones(i).gamma_exp = val2;
    }
    rdline, f;
    rdline, f;

    // Cavity
    read, f, string1 ;  param.cavity.lcavity = string1;
    read, f, val1 ;  param.cavity.H0 = val1; param.cavity.R0 = val2 ;
    read, f, val1 ;  param.cavity.beta = val1;
    rdline, f;
    rdline, f;

    //Grain
    param.simu.n_especes_max = 0 ;
    param.simu.n_components_max = 0 ;

    for (i=1 ; i <=  param.simu.n_zones ; i++) {
      read, f, val1 ; param.zones(i).n_especes=val1;
      if (val1 > param.simu.n_especes_max) param.simu.n_especes_max = val1 ;

      for (j=1 ; j<= param.zones(i).n_especes ; j++) {
        read, f, string1, val1, val2, val3, val4, val5 ;
        if (param.zones(i).n_especes <= n_especes) {
          param.dust_pop(i,j).type = string1 ;
          param.dust_pop(i,j).n_components = val1 ; if (val1 > param.simu.n_components_max) param.simu.n_components_max = val1 ;
          param.dust_pop(i,j).mixing_rule = val2 ;
          param.dust_pop(i,j).porosity = val3;
          param.dust_pop(i,j).mass_fraction = val4;
          param.dust_pop(i,j).dhs_maxf = val5;
        }

        if (n_components < param.simu.n_components_max) return param

        for (k=1 ; k<=param.dust_pop(i,j).n_components ; k++) {
          read, f, string1, val1 ;
          if (param.dust_pop(i,j).n_components <= n_components) {
            param.dust_pop(i,j).file(k) = string1;  param.dust_pop(i,j).component_volume_fraction(k) = val1 ;
          }
        }
        read, f, val1; param.dust_pop(i,j).type_grain = val1;
        read, f, val1, val2, val3, val4 ; param.dust_pop(i,j).amin=val1;  param.dust_pop(i,j).amax=val2;  param.dust_pop(i,j).aexp=val3; param.dust_pop(i,j).n_grains=val4;
      }
    }


    if (param.dust_pop(1,1).type_grain == 2) {
      param.gridT.n_T_grain = sum(param.dust_pop.n_grains);
    } else {
      param.gridT.n_T_grain = 1;
    }
    rdline, f;
    rdline, f;


    // mol
    read, f, string1, string2, string3, val1; param.mol.lpop=string1 ; param.mol.lpop_precise=string2 ; param.mol.lLTE=string3 ; param.mol.profile_width = val1;

    read, f, val1 ; param.mol.vturb=val1 ;
    read, f, val1 ; param.simu.nmol=val1 ;

    if (param.simu.nmol == nmol)  {
      for (i=1 ; i <=  param.simu.nmol ; i++) {
        read, f, string1, val1; param.mol(i).molecule_file=string1 ; param.mol(i).level_max = val1;
        read, f, val1, val2 ;  param.mol(i).vmax=val1 ;  param.mol(i).n_speed = val2;
        read, f, string1, val1, string2 ; param.mol(i).cst_abundance=string1 ; param.mol(i).abundance=val1 ; param.mol(i).abundance_file = string2;
        read, f, string1, val1; param.mol(i).ray_tracing=string1 ; param.mol(i).n_lines_ray_tracing=val1;
        if (param.mol(i).n_lines_ray_tracing == 1) {
          read, f, val1 ; param.mol(i).transition_numbers(1) = val1;
        } else if (param.mol(i).n_lines_ray_tracing == 2) {
          read, f, val1, val2 ; param.mol(i).transition_numbers(1:2) = [val1,val2];
        } else if (param.mol(i).n_lines_ray_tracing == 3) {
          read, f, val1, val2, val3 ; param.mol(i).transition_numbers(1:3) = [val1,val2,val3];
        } else if (param.mol(i).n_lines_ray_tracing == 4) {
          read, f, val1, val2, val3, val4 ; param.mol(i).transition_numbers(1:4) = [val1,val2,val3,val4];
        } else if (param.mol(i).n_lines_ray_tracing == 5) {
          read, f, val1, val2, val3, val4, val5 ; param.mol(i).transition_numbers(1:5) = [val1,val2,val3,val4,val5];
        }
      }
    }
    rdline, f;
    rdline, f;


    // Star
    read, f, val1 ; param.simu.n_stars = val1;
    if (param.simu.n_stars == n_stars) {
      for (i=1 ; i <= param.simu.n_stars ; i++) {
        read, f, val1, val2, val3, val4, val5, val6, string1, val7; param.stars(i).T = val1; param.stars(i).R = val2 ;  param.stars(i).Mass = val3;
        param.stars(i).x = val4 ; param.stars(i).y = val5 ; param.stars(i).z = val6 ;
        param.stars(i).is_bb = string1;
        read, f, param.stars(i).spectre ;
        read, f, val1, val2 ;
        param.stars(i).fUV = val1 ;
        param.stars(i).slope_fUV = val2 ;
      }
    }

  } else if (abs(param.simu.version -2.19) < 1.0e-6) {
    rdline, f;

    //photon;
    read, f, val1 ; param.phot.n2 = val1;
    read, f, val1 ; param.phot.nlambda = val1;
    read, f, val1 ; param.phot.nimage = val1;

    rdline, f;
    rdline, f;

    //Wavelength
    read, f, val1, val2, val3 ;
    param.wave.nlamb1 = val1; param.wave.lambda_min = val2; param.wave.lambda_max = val3; lambda_min=val2 ; lambda_max=val3;
    read, f, string1, string2, string3 ; param.simu.ltemp = string1; param.simu.lsed = string2; param.simu.lsed_complete = string3;
    read, f, string1 ; param.wave.file = string1;
    param.simu.lem_disk_image = "T";
    read, f, string1, string2 ; param.simu.lsepar = string1; param.simu.lsepar_pola = string2;
    param.map.n_type=1;
    if (param.simu.lsepar=="T") param.map.n_type += 4;
    if (param.simu.lsepar_pola=="T") param.map.n_type += 3;

    rdline, f;
    rdline, f;

    //Grid
    read, f, val1; param.grid.geometry = val1;
    read, f, val1, val2, val3, val4 ;
    param.grid.nrad=val1; param.grid.nz=val2; param.grid.n_az=val3 ; param.grid.nrad_in=val4;

    rdline, f;
    rdline, f;

    //map;
    read, f, val1, val2, val3 ;
    param.map.nx = val1; param.map.ny=val2; param.map.map_size=val3; param.map.zoom = 1.0 ;
    param.map.size_neb = 0.5 * param.map.map_size ;
    read, f, val1, val2, val3 ;
    param.map.nthet = val1;   param.map.nphi = val2  ; param.map.inc = val3;
    param.map.delta = 1 ; param.map.angle_interest = 75.; param.map.lonly_bin_interest="F"
    read, f, val1, val2, val3, string1;
    param.map.RT_imin = val1 ; param.map.RT_imax = val2 ; param.map.RT_ntheta = val3 ; param.map.lRT_i_centered = string1 ;
    read, f, val1 ; param.map.dist = val1;
    read, f, val1 ; param.map.PA = val1;

    rdline, f;
    rdline, f;

    // Scattering method
    read, f, val1 ; param.simu.scatt_method = val1;
    read, f, val1 ; param.simu.aniso_method = val1;

    rdline, f;
    rdline, f;

    // Symetries
    read, f, string1; param.simu.lsym_ima = string1 ;
    read, f, string1; param.simu.lsym_centrale = string1 ;
    read, f, string1; param.simu.lsym_axiale = string1 ;

    rdline, f;
    rdline, f;

    // Dust global properties
    read, f, val1, val2, val3, com2; param.dust.strat_type = val1 ;  param.dust.strat = val2; param.dust.a_strat = val3;
    if (val1 > 0) {
      param.dust.lstrat = "T" ;
    } else {
      param.dust.lstrat = "F" ;
    }
    read, f, string1, val1; param.dust.lmigration = string1 ;
    read, f, string1, val1; param.simu.ldust_sublimation = string1 ; param.simu.correct_Rsub = val1 ;
    read, f, string1, val1; param.simu.lhydro_eq = string1 ;
    read, f, string1 , val1 ; param.simu.lchauff_int = string1 ; param.simu.viscosity = val1;

    read, f;
    read, f;

    //Nbre de zones
    read,f, val1 ; param.simu.n_zones = val1;

    // On arrete ici si on n'a pas le bon nombre de zones et on recharche la structure
    if (n_zones != param.simu.n_zones) return param

    read, f;
    read, f;

    //Disk
    for (i=1 ; i <=  param.simu.n_zones ; i++) {
      read, f, val1 ; param.zones(i).geometry = val1 ;
      read, f, val1, val2; param.zones(i).mass = val1; param.zones(i).gas_to_dust = val2;
      read, f, val1, val2, val3; param.zones(i).Ho = val1; param.zones(i).Ro = val2; param.zones(i).vert_exp = val3 ;
      read, f, val1, val2, val3, val4;  param.zones(i).rin = val1; param.zones(i).edge = val2; param.zones(i).rout = val3;  param.zones(i).Rc = val4;
      read, f, val1; param.zones(i).beta = val1;
      read, f, val1, val2; param.zones(i).dens = val1; param.zones(i).gamma_exp = val2;
    }
    rdline, f;
    rdline, f;

    // Cavity
    read, f, string1 ;  param.cavity.lcavity = string1;
    read, f, val1 ;  param.cavity.H0 = val1; param.cavity.R0 = val2 ;
    read, f, val1 ;  param.cavity.beta = val1;
    rdline, f;
    rdline, f;

    //Grain
    param.simu.n_especes_max = 0 ;
    param.simu.n_components_max = 0 ;

    for (i=1 ; i <=  param.simu.n_zones ; i++) {
      read, f, val1 ; param.zones(i).n_especes=val1;
      if (val1 > param.simu.n_especes_max) param.simu.n_especes_max = val1 ;

      for (j=1 ; j<= param.zones(i).n_especes ; j++) {
        read, f, string1, val1, val2, val3, val4, val5 ;
        if (param.zones(i).n_especes <= n_especes) {
          param.dust_pop(i,j).type = string1 ;
          param.dust_pop(i,j).n_components = val1 ; if (val1 > param.simu.n_components_max) param.simu.n_components_max = val1 ;
          param.dust_pop(i,j).mixing_rule = val2 ;
          param.dust_pop(i,j).porosity = val3;
          param.dust_pop(i,j).mass_fraction = val4;
          param.dust_pop(i,j).dhs_maxf = val5;
        }

        if (n_components != param.simu.n_components_max) return param

        for (k=1 ; k<=param.dust_pop(i,j).n_components ; k++) {
          read, f, string1, val1 ;
          if (param.dust_pop(i,j).n_components <= n_components) {
            param.dust_pop(i,j).file(k) = string1;  param.dust_pop(i,j).component_volume_fraction(k) = val1 ;
          }
        }
        read, f, val1; param.dust_pop(i,j).type_grain = val1;
        read, f, val1, val2, val3, val4 ; param.dust_pop(i,j).amin=val1;  param.dust_pop(i,j).amax=val2;  param.dust_pop(i,j).aexp=val3; param.dust_pop(i,j).n_grains=val4;
      }
    }


    if (param.dust_pop(1,1).type_grain == 2) {
      param.gridT.n_T_grain = sum(param.dust_pop.n_grains);
    } else {
      param.gridT.n_T_grain = 1;
    }
    rdline, f;
    rdline, f;


    // mol
    read, f, string1, string2, string3, val1; param.mol.lpop=string1 ; param.mol.lpop_precise=string2 ; param.mol.lLTE=string3 ; param.mol.profile_width = val1;

    read, f, val1 ; param.mol.vturb=val1 ;
    read, f, val1 ; param.simu.nmol=val1 ;

    if (param.simu.nmol == nmol)  {
      for (i=1 ; i <=  param.simu.nmol ; i++) {
        read, f, string1, val1; param.mol(i).molecule_file=string1 ; param.mol(i).level_max = val1;
        read, f, val1, val2 ;  param.mol(i).vmax=val1 ;  param.mol(i).n_speed = val2;
        read, f, string1, val1, string2 ; param.mol(i).cst_abundance=string1 ; param.mol(i).abundance=val1 ; param.mol(i).abundance_file = string2;
        read, f, string1, val1; param.mol(i).ray_tracing=string1 ; param.mol(i).n_lines_ray_tracing=val1;
        if (param.mol(i).n_lines_ray_tracing == 1) {
          read, f, val1 ; param.mol(i).transition_numbers(1) = val1;
        } else if (param.mol(i).n_lines_ray_tracing == 2) {
          read, f, val1, val2 ; param.mol(i).transition_numbers(1:2) = [val1,val2];
        } else if (param.mol(i).n_lines_ray_tracing == 3) {
          read, f, val1, val2, val3 ; param.mol(i).transition_numbers(1:3) = [val1,val2,val3];
        } else if (param.mol(i).n_lines_ray_tracing == 4) {
          read, f, val1, val2, val3, val4 ; param.mol(i).transition_numbers(1:4) = [val1,val2,val3,val4];
        } else if (param.mol(i).n_lines_ray_tracing == 5) {
          read, f, val1, val2, val3, val4, val5 ; param.mol(i).transition_numbers(1:5) = [val1,val2,val3,val4,val5];
        }
      }
    }
    rdline, f;
    rdline, f;


    // Star
    read, f, val1 ; param.simu.n_stars = val1;
    if (param.simu.n_stars == n_stars) {
      for (i=1 ; i <= param.simu.n_stars ; i++) {
        read, f, val1, val2, val3, val4, val5, val6, string1, val7; param.stars(i).T = val1; param.stars(i).R = val2 ;  param.stars(i).Mass = val3;
        param.stars(i).x = val4 ; param.stars(i).y = val5 ; param.stars(i).z = val6 ;
        param.stars(i).is_bb = string1;
        read, f, param.stars(i).spectre ;
        read, f, val1, val2 ;
        param.stars(i).fUV = val1 ;
        param.stars(i).slope_fUV = val2 ;
      }
    }


  } else if (abs(param.simu.version -2.18) < 1.0e-6) {
    rdline, f;

    //photon;
    read, f, val1 ; param.phot.n2 = val1;
    read, f, val1 ; param.phot.nlambda = val1;
    read, f, val1 ; param.phot.nimage = val1;

    rdline, f;
    rdline, f;

    //Wavelength
    read, f, val1, val2, val3 ;
    param.wave.nlamb1 = val1; param.wave.lambda_min = val2; param.wave.lambda_max = val3; lambda_min=val2 ; lambda_max=val3;
    read, f, string1, string2, string3 ; param.simu.ltemp = string1; param.simu.lsed = string2; param.simu.lsed_complete = string3;
    read, f, string1 ; param.wave.file = string1;
    param.simu.lem_disk_image = "T";
    read, f, string1, string2 ; param.simu.lsepar = string1; param.simu.lsepar_pola = string2;
    param.map.n_type=1;
    if (param.simu.lsepar=="T") param.map.n_type += 4;
    if (param.simu.lsepar_pola=="T") param.map.n_type += 3;

    rdline, f;
    rdline, f;

    //Grid
    read, f, val1; param.grid.geometry = val1;
    read, f, val1, val2, val3, val4 ;
    param.grid.nrad=val1; param.grid.nz=val2; param.grid.n_az=val3 ; param.grid.nrad_in=val4;

    rdline, f;
    rdline, f;

    //map;
    read, f, val1, val2, val3 ;
    param.map.nx = val1; param.map.ny=val2; param.map.map_size=val3; param.map.zoom = 1.0 ;
    param.map.size_neb = 0.5 * param.map.map_size ;
    read, f, val1, val2, val3 ;
    param.map.nthet = val1;   param.map.nphi = val2  ; param.map.inc = val3;
    param.map.delta = 1 ; param.map.angle_interest = 75.; param.map.lonly_bin_interest="F"
    read, f, val1, val2, val3, string1;
    param.map.RT_imin = val1 ; param.map.RT_imax = val2 ; param.map.RT_ntheta = val3 ; param.map.lRT_i_centered = string1 ;
    read, f, val1 ; param.map.dist = val1;
    read, f, val1 ; param.map.PA = val1;

    rdline, f;
    rdline, f;

    // Scattering method
    read, f, val1 ; param.simu.scatt_method = val1;
    read, f, val1 ; param.simu.aniso_method = val1;

    rdline, f;
    rdline, f;

    // Symetries
    read, f, string1; param.simu.lsym_ima = string1 ;
    read, f, string1; param.simu.lsym_centrale = string1 ;
    read, f, string1; param.simu.lsym_axiale = string1 ;

    rdline, f;
    rdline, f;

    // Dust global properties
    read, f, val1, val2, val3, com2; param.dust.strat_type = val1 ;  param.dust.strat = val2; param.dust.a_strat = val3;
    if (val1 > 0) {
      param.dust.lstrat = "T" ;
    } else {
      param.dust.lstrat = "F" ;
    }
    read, f, string1, val1; param.dust.lmigration = string1 ;
    read, f, string1, val1; param.simu.ldust_sublimation = string1 ; param.simu.correct_Rsub = val1 ;
    read, f, string1, val1; param.simu.lhydro_eq = string1 ;
    read, f, string1 , val1 ; param.simu.lchauff_int = string1 ; param.simu.viscosity = val1;

    read, f;
    read, f;

    //Nbre de zones
    read,f, val1 ; param.simu.n_zones = val1;

    // On arrete ici si on n'a pas le bon nombre de zones et on recharche la structure
    if (n_zones != param.simu.n_zones) return param

    read, f;
    read, f;

    //Disk
    for (i=1 ; i <=  param.simu.n_zones ; i++) {
      read, f, val1 ; param.zones(i).geometry = val1 ;
      read, f, val1, val2; param.zones(i).mass = val1; param.zones(i).gas_to_dust = val2;
      read, f, val1, val2; param.zones(i).Ho = val1; param.zones(i).Ro = val2;
      read, f, val1, val2, val3, val4;  param.zones(i).rin = val1; param.zones(i).edge = val2; param.zones(i).rout = val3;  param.zones(i).Rc = val4;
      read, f, val1; param.zones(i).beta = val1;
      read, f, val1, val2; param.zones(i).dens = val1; param.zones(i).gamma_exp = val2;
    }
    rdline, f;
    rdline, f;

    // Cavity
    read, f, string1 ;  param.cavity.lcavity = string1;
    read, f, val1 ;  param.cavity.H0 = val1; param.cavity.R0 = val2 ;
    read, f, val1 ;  param.cavity.beta = val1;
    rdline, f;
    rdline, f;

    //Grain
    param.simu.n_especes_max = 0 ;
    param.simu.n_components_max = 0 ;

    for (i=1 ; i <=  param.simu.n_zones ; i++) {
      read, f, val1 ; param.zones(i).n_especes=val1;
      if (val1 > param.simu.n_especes_max) param.simu.n_especes_max = val1 ;

      for (j=1 ; j<= param.zones(i).n_especes ; j++) {
        read, f, string1, val1, val2, val3, val4, val5 ;
        if (param.zones(i).n_especes <= n_especes) {
          param.dust_pop(i,j).type = string1 ;
          param.dust_pop(i,j).n_components = val1 ;  if (val1 > param.simu.n_components_max) param.simu.n_components_max = val1 ;
          param.dust_pop(i,j).mixing_rule = val2 ;
          param.dust_pop(i,j).porosity = val3;
          param.dust_pop(i,j).mass_fraction = val4;
          param.dust_pop(i,j).dhs_maxf = val5;
        }

        for (k=1 ; k<=param.dust_pop(i,j).n_components ; k++) {
          read, f, string1, val1 ;
          if (param.dust_pop(i,j).n_components <= n_components) {
            param.dust_pop(i,j).file(k) = string1;  param.dust_pop(i,j).component_volume_fraction(k) = val1 ;
          }
        }
        read, f, val1; param.dust_pop(i,j).type_grain = val1;
        read, f, val1, val2, val3, val4 ; param.dust_pop(i,j).amin=val1;  param.dust_pop(i,j).amax=val2;  param.dust_pop(i,j).aexp=val3; param.dust_pop(i,j).n_grains=val4;
      }
    }


    if (param.dust_pop(1,1).type_grain == 2) {
      param.gridT.n_T_grain = sum(param.dust_pop.n_grains);
    } else {
      param.gridT.n_T_grain = 1;
    }
    rdline, f;
    rdline, f;


    // mol
    read, f, string1, string2, string3, val1; param.mol.lpop=string1 ; param.mol.lpop_precise=string2 ; param.mol.lLTE=string3 ; param.mol.profile_width = val1;

    read, f, val1 ; param.mol.vturb=val1 ;
    read, f, val1 ; param.simu.nmol=val1 ;

    if (param.simu.nmol == nmol)  {
      for (i=1 ; i <=  param.simu.nmol ; i++) {
        read, f, string1, val1; param.mol(i).molecule_file=string1 ; param.mol(i).level_max = val1;
        read, f, val1, val2 ;  param.mol(i).vmax=val1 ;  param.mol(i).n_speed = val2;
        read, f, string1, val1, string2 ; param.mol(i).cst_abundance=string1 ; param.mol(i).abundance=val1 ; param.mol(i).abundance_file = string2;
        read, f, string1, val1; param.mol(i).ray_tracing=string1 ; param.mol(i).n_lines_ray_tracing=val1;
        if (param.mol(i).n_lines_ray_tracing == 1) {
          read, f, val1 ; param.mol(i).transition_numbers(1) = val1;
        } else if (param.mol(i).n_lines_ray_tracing == 2) {
          read, f, val1, val2 ; param.mol(i).transition_numbers(1:2) = [val1,val2];
        } else if (param.mol(i).n_lines_ray_tracing == 3) {
          read, f, val1, val2, val3 ; param.mol(i).transition_numbers(1:3) = [val1,val2,val3];
        } else if (param.mol(i).n_lines_ray_tracing == 4) {
          read, f, val1, val2, val3, val4 ; param.mol(i).transition_numbers(1:4) = [val1,val2,val3,val4];
        } else if (param.mol(i).n_lines_ray_tracing == 5) {
          read, f, val1, val2, val3, val4, val5 ; param.mol(i).transition_numbers(1:5) = [val1,val2,val3,val4,val5];
        }
      }
    }
    rdline, f;
    rdline, f;


    // Star
    read, f, val1 ; param.simu.n_stars = val1;
    if (param.simu.n_stars == n_stars) {
      for (i=1 ; i <= param.simu.n_stars ; i++) {
        read, f, val1, val2, val3, val4, val5, val6, string1, val7; param.stars(i).T = val1; param.stars(i).R = val2 ;  param.stars(i).Mass = val3;
        param.stars(i).x = val4 ; param.stars(i).y = val5 ; param.stars(i).z = val6 ;
        param.stars(i).is_bb = string1;
        read, f, param.stars(i).spectre ;
        read, f, val1, val2 ;
        param.stars(i).fUV = val1 ;
        param.stars(i).slope_fUV = val2 ;
      }
    }


  } else if (abs(param.simu.version -2.17) < 1.0e-6) {
    rdline, f;

    //photon;
    read, f, val1 ; param.phot.n2 = val1;
    read, f, val1 ; param.phot.nlambda = val1;
    read, f, val1 ; param.phot.nimage = val1;

    rdline, f;
    rdline, f;

    //Wavelength
    read, f, val1, val2, val3 ;
    param.wave.nlamb1 = val1; param.wave.lambda_min = val2; param.wave.lambda_max = val3; lambda_min=val2 ; lambda_max=val3;
    read, f, string1, string2, string3 ; param.simu.ltemp = string1; param.simu.lsed = string2; param.simu.lsed_complete = string3;
    read, f, string1 ; param.wave.file = string1;
    param.simu.lem_disk_image = "T";
    read, f, string1, string2 ; param.simu.lsepar = string1; param.simu.lsepar_pola = string2;
    param.map.n_type=1;
    if (param.simu.lsepar=="T") param.map.n_type += 4;
    if (param.simu.lsepar_pola=="T") param.map.n_type += 3;

    rdline, f;
    rdline, f;

    //Grid
    read, f, val1; param.grid.geometry = val1;
    read, f, val1, val2, val3, val4 ;
    param.grid.nrad=val1; param.grid.nz=val2; param.grid.n_az=val3 ; param.grid.nrad_in=val4;

    rdline, f;
    rdline, f;

    //map;
    read, f, val1, val2, val3, val4 ;
    param.map.nx = val1; param.map.ny=val2; param.map.map_size=val3; param.map.zoom = val4 ;
    param.map.size_neb = 0.5 * param.map.map_size ;
    read, f, val1, val2, val3 ;
    param.map.nthet = val1;   param.map.nphi = val2  ; param.map.inc = val3;
    param.map.delta = 1 ; param.map.angle_interest = 75.; param.map.lonly_bin_interest="F"
    read, f, val1, val2, val3, string1;
    param.map.RT_imin = val1 ; param.map.RT_imax = val2 ; param.map.RT_ntheta = val3 ; param.map.lRT_i_centered = string1 ;
    read, f, val1 ; param.map.dist = val1;
    read, f, val1 ; param.map.PA = val1;

    rdline, f;
    rdline, f;

    // Scattering method
    read, f, val1 ; param.simu.scatt_method = val1;
    read, f, val1 ; param.simu.aniso_method = val1;

    rdline, f;
    rdline, f;

    // Symetries
    read, f, string1; param.simu.lsym_ima = string1 ;
    read, f, string1; param.simu.lsym_centrale = string1 ;
    read, f, string1; param.simu.lsym_axiale = string1 ;

    rdline, f;
    rdline, f;

    // Dust global properties
    read, f, val1, val2, val3, com2; param.dust.strat_type = val1 ;  param.dust.strat = val2; param.dust.a_strat = val3;
    if (val1 > 0) {
      param.dust.lstrat = "T" ;
    } else {
      param.dust.lstrat = "F" ;
    }
    read, f, string1, val1; param.simu.ldust_sublimation = string1 ; param.simu.correct_Rsub = val1 ;
    read, f, string1 , val1 ; param.simu.lchauff_int = string1 ; param.simu.viscosity = val1;

    read, f;
    read, f;

    //Nbre de zones
    read,f, val1 ; param.simu.n_zones = val1;

    // On arrete ici si on n'a pas le bon nombre de zones et on recharche la structure
    if (n_zones != param.simu.n_zones) return param

    read, f;
    read, f;

    //Disk
    for (i=1 ; i <=  param.simu.n_zones ; i++) {
      read, f, val1 ; param.zones(i).geometry = val1 ;
      read, f, val1, val2; param.zones(i).mass = val1; param.zones(i).gas_to_dust = val2;
      read, f, val1, val2; param.zones(i).Ho = val1; param.zones(i).Ro = val2;
      read, f, val1, val2, val3, val4;  param.zones(i).rin = val1; param.zones(i).rout = val2; param.zones(i).edge = val3;
      read, f, val1; param.zones(i).beta = val1;
      read, f, val1, val2; param.zones(i).dens = val1; param.zones(i).gamma_exp = val2;
    }
    rdline, f;
    rdline, f;

    // Cavity
    read, f, string1 ;  param.cavity.lcavity = string1;
    read, f, val1 ;  param.cavity.H0 = val1; param.cavity.R0 = val2 ;
    read, f, val1 ;  param.cavity.beta = val1;
    rdline, f;
    rdline, f;

    //Grain
    param.simu.n_especes_max = 0 ;
    param.simu.n_components_max = 0 ;

    for (i=1 ; i <=  param.simu.n_zones ; i++) {
      read, f, val1 ; param.zones(i).n_especes=val1;
      if (val1 > param.simu.n_especes_max) param.simu.n_especes_max = val1 ;

      for (j=1 ; j<= param.zones(i).n_especes ; j++) {
        read, f, string1, val1, val2, val3, val4, val5 ;
        if (param.zones(i).n_especes <= n_especes) {
          param.dust_pop(i,j).type = string1 ;
          param.dust_pop(i,j).n_components = val1 ;  if (val1 > param.simu.n_components_max) param.simu.n_components_max = val1 ;
          param.dust_pop(i,j).mixing_rule = val2 ;
          param.dust_pop(i,j).porosity = val3;
          param.dust_pop(i,j).mass_fraction = val4;
          param.dust_pop(i,j).dhs_maxf = val5;
        }

        for (k=1 ; k<=param.dust_pop(i,j).n_components ; k++) {
          read, f, string1, val1 ;
          if (param.dust_pop(i,j).n_components <= n_components) {
            param.dust_pop(i,j).file(k) = string1;  param.dust_pop(i,j).component_volume_fraction(k) = val1 ;
          }
        }
        read, f, val1; param.dust_pop(i,j).type_grain = val1;
        read, f, val1, val2, val3, val4 ; param.dust_pop(i,j).amin=val1;  param.dust_pop(i,j).amax=val2;  param.dust_pop(i,j).aexp=val3; param.dust_pop(i,j).n_grains=val4;
      }
    }


    if (param.dust_pop(1,1).type_grain == 2) {
      param.gridT.n_T_grain = sum(param.dust_pop.n_grains);
    } else {
      param.gridT.n_T_grain = 1;
    }
    rdline, f;
    rdline, f;


    // mol
    read, f, string1, string2, string3, val1; param.mol.lpop=string1 ; param.mol.lpop_precise=string2 ; param.mol.lLTE=string3 ; param.mol.profile_width = val1;

    read, f, val1 ; param.mol.vturb=val1 ;
    read, f, val1 ; param.simu.nmol=val1 ;

    if (param.simu.nmol == nmol)  {
      for (i=1 ; i <=  param.simu.nmol ; i++) {
        read, f, string1, val1; param.mol(i).molecule_file=string1 ; param.mol(i).level_max = val1;
        read, f, val1, val2 ;  param.mol(i).vmax=val1 ;  param.mol(i).n_speed = val2;
        read, f, string1, val1, string2 ; param.mol(i).cst_abundance=string1 ; param.mol(i).abundance=val1 ; param.mol(i).abundance_file = string2;
        read, f, string1, val1; param.mol(i).ray_tracing=string1 ; param.mol(i).n_lines_ray_tracing=val1;
        if (param.mol(i).n_lines_ray_tracing == 1) {
          read, f, val1 ; param.mol(i).transition_numbers(1) = val1;
        } else if (param.mol(i).n_lines_ray_tracing == 2) {
          read, f, val1, val2 ; param.mol(i).transition_numbers(1:2) = [val1,val2];
        } else if (param.mol(i).n_lines_ray_tracing == 3) {
          read, f, val1, val2, val3 ; param.mol(i).transition_numbers(1:3) = [val1,val2,val3];
        } else if (param.mol(i).n_lines_ray_tracing == 4) {
          read, f, val1, val2, val3, val4 ; param.mol(i).transition_numbers(1:4) = [val1,val2,val3,val4];
        } else if (param.mol(i).n_lines_ray_tracing == 5) {
          read, f, val1, val2, val3, val4, val5 ; param.mol(i).transition_numbers(1:5) = [val1,val2,val3,val4,val5];
        }
      }
    }
    rdline, f;
    rdline, f;


    // Star
    read, f, val1 ; param.simu.n_stars = val1;
    if (param.simu.n_stars == n_stars) {
      for (i=1 ; i <= param.simu.n_stars ; i++) {
        read, f, val1, val2, val3, val4, val5, val6, string1, val7; param.stars(i).T = val1; param.stars(i).R = val2 ;  param.stars(i).Mass = val3;
        param.stars(i).x = val4 ; param.stars(i).y = val5 ; param.stars(i).z = val6 ;
        param.stars(i).is_bb = string1;
        read, f, param.stars(i).spectre ;
        read, f, val1, val2 ;
        param.stars(i).fUV = val1 ;
        param.stars(i).slope_fUV = val2 ;
      }
    }


  } else if (abs(param.simu.version -2.16) < 1.0e-6) {
    rdline, f;

    //photon;
    read, f, val1; param.phot.n2 = val1;
    read, f, val1; param.phot.nlambda = val1;
    read, f, val1; param.phot.nimage = val1;

    rdline, f;
    rdline, f;

    //Wavelength
    read, f, val1, val2, val3 ;
    param.wave.nlamb1 = val1; param.wave.lambda_min = val2; param.wave.lambda_max = val3; lambda_min=val2 ; lambda_max=val3;
    read, f, string1, string2, string3 ; param.simu.ltemp = string1; param.simu.lsed = string2; param.simu.lsed_complete = string3;
    read, f, string1, com2; param.wave.file = string1;
    param.simu.lem_disk_image = "T";
    read, f, string1, string2 ; param.simu.lsepar = string1; param.simu.lsepar_pola = string2;
    param.map.n_type=1;
    if (param.simu.lsepar=="T") param.map.n_type += 4;
    if (param.simu.lsepar_pola=="T") param.map.n_type += 3;

    rdline, f;
    rdline, f;

    //Grid
    read, f, val1; param.grid.geometry = val1;
    read, f, val1, val2, val3, val4;
    param.grid.nrad=val1; param.grid.nz=val2; param.grid.n_az=val3 ; param.grid.nrad_in=val4;

    rdline, f;
    rdline, f;

    //map;
    read, f, val1, val2, val3, val4;
    param.map.nx = val1; param.map.ny=val2; param.map.map_size=val3; param.map.zoom = val4 ;
    param.map.size_neb = 0.5 * param.map.map_size ;
    read, f, val1, val2, val3;
    param.map.nthet = val1;   param.map.nphi = val2  ; param.map.inc = val3;
    param.map.delta = 1 ; param.map.angle_interest = 75.; param.map.lonly_bin_interest="F"
    read, f, val1, val2, val3, string1;
    param.map.RT_imin = val1 ; param.map.RT_imax = val2 ; param.map.RT_ntheta = val3 ; param.map.lRT_i_centered = string1 ;
    read, f, val1 ; param.map.dist = val1;
    read, f, val1 ; param.map.PA = val1;

    rdline, f;
    rdline, f;

    // Scattering method
    read, f, val1 ; param.simu.scatt_method = val1;
    read, f, val1 ; param.simu.aniso_method = val1;

    rdline, f;
    rdline, f;

    // Symetries
    read, f, string1; param.simu.lsym_ima = string1 ;
    read, f, string1; param.simu.lsym_centrale = string1 ;
    read, f, string1; param.simu.lsym_axiale = string1 ;

    rdline, f;
    rdline, f;

    // Dust global properties
    read, f, string1, val1, val2, com2; param.dust.lstrat = string1 ;  param.dust.strat = val1; param.dust.a_strat = val2;
    if (param.dust.lstrat=="T") {
      param.dust.strat_type=1 ;
    } else {
      param.dust.strat_type=0 ;
    }
    read, f, string1, val1; param.simu.ldust_sublimation = string1 ; param.simu.correct_Rsub = val1 ;
    read, f, string1 , val1 ; param.simu.lchauff_int = string1 ; param.simu.viscosity = val1;

    read, f;
    read, f;

    //Nbre de zones
    read,f, val1 ; param.simu.n_zones = val1;

    // On arrete ici si on n'a pas le bon nombre de zones et on recharche la structure
    if (n_zones != param.simu.n_zones) return param

    read, f;
    read, f;

    //Disk
    for (i=1 ; i <=  param.simu.n_zones ; i++) {
      read, f, val1 ; param.zones(i).geometry = val1 ;
      read, f, val1, val2; param.zones(i).mass = val1; param.zones(i).gas_to_dust = val2;
      read, f, val1, val2; param.zones(i).Ho = val1; param.zones(i).Ro = val2;
      read, f, val1, val2, val3, val4;  param.zones(i).rin = val1; param.zones(i).rout = val2; param.zones(i).edge = val3;
      read, f, val1; param.zones(i).beta = val1;
      read, f, val1; param.zones(i).dens = val1;
    }
    rdline, f;
    rdline, f;

    // Cavity
    read, f, string1 ;  param.cavity.lcavity = string1;
    read, f, val1 ;  param.cavity.H0 = val1; param.cavity.R0 = val2 ;
    read, f, val1 ;  param.cavity.beta = val1;
    rdline, f;
    rdline, f;

    //Grain
    param.simu.n_especes_max = 0 ;
    param.simu.n_components_max = 0 ;

    for (i=1 ; i <=  param.simu.n_zones ; i++) {
      read, f, val1; param.zones(i).n_especes=val1;
      if (val1 > param.simu.n_especes_max) param.simu.n_especes_max = val1 ;

      for (j=1 ; j<= param.zones(i).n_especes ; j++) {
        read, f, string1, val1, val2, val3, val4 ;
        if (param.zones(i).n_especes <= n_especes) {
          param.dust_pop(i,j).type = string1 ;
          param.dust_pop(i,j).n_components = val1 ;  if (val1 > param.simu.n_components_max) param.simu.n_components_max = val1 ;
          param.dust_pop(i,j).mixing_rule = val2 ;
          param.dust_pop(i,j).porosity = val3;
          param.dust_pop(i,j).mass_fraction = val4;
        }

        for (k=1 ; k<=param.dust_pop(i,j).n_components ; k++) {
          read, f, string1, val1 ;
          if (param.dust_pop(i,j).n_components <= n_components) {
            param.dust_pop(i,j).file(k) = string1;  param.dust_pop(i,j).component_volume_fraction(k) = val1 ;
          }
        }
        read, f, val1; param.dust_pop(i,j).type_grain = val1;
        read, f, val1, val2, val3, val4 ; param.dust_pop(i,j).amin=val1;  param.dust_pop(i,j).amax=val2;  param.dust_pop(i,j).aexp=val3; param.dust_pop(i,j).n_grains=val4;
      }
    }


    if (param.dust_pop(1,1).type_grain == 2) {
      param.gridT.n_T_grain = sum(param.dust_pop.n_grains);
    } else {
      param.gridT.n_T_grain = 1;
    }
    rdline, f;
    rdline, f;


    // mol
    read, f, string1, string2, string3, val1; param.mol.lpop=string1 ; param.mol.lpop_precise=string2 ; param.mol.lLTE=string3 ; param.mol.profile_width = val1;

    read, f, val1 ; param.mol.vturb=val1 ;
    read, f, val1 ; param.simu.nmol=val1 ;

    if (param.simu.nmol == nmol)  {
      for (i=1 ; i <=  param.simu.nmol ; i++) {
        read, f, string1, val1; param.mol(i).molecule_file=string1 ; param.mol(i).level_max = val1;
        read, f, val1, val2 ;  param.mol(i).vmax=val1 ;  param.mol(i).n_speed = val2;
        read, f, string1, val1, string2 ; param.mol(i).cst_abundance=string1 ; param.mol(i).abundance=val1 ; param.mol(i).abundance_file = string2;
        read, f, string1, val1; param.mol(i).ray_tracing=string1 ; param.mol(i).n_lines_ray_tracing=val1;
        if (param.mol(i).n_lines_ray_tracing == 1) {
          read, f, val1 ; param.mol(i).transition_numbers(1) = val1;
        } else if (param.mol(i).n_lines_ray_tracing == 2) {
          read, f, val1, val2 ; param.mol(i).transition_numbers(1:2) = [val1,val2];
        } else if (param.mol(i).n_lines_ray_tracing == 3) {
          read, f, val1, val2, val3 ; param.mol(i).transition_numbers(1:3) = [val1,val2,val3];
        } else if (param.mol(i).n_lines_ray_tracing == 4) {
          read, f, val1, val2, val3, val4 ; param.mol(i).transition_numbers(1:4) = [val1,val2,val3,val4];
        } else if (param.mol(i).n_lines_ray_tracing == 5) {
          read, f, val1, val2, val3, val4, val5 ; param.mol(i).transition_numbers(1:5) = [val1,val2,val3,val4,val5];
        }
      }
    }
    rdline, f;
    rdline, f;


    // Star
    read, f, val1; param.simu.n_stars = val1;
    if (param.simu.n_stars == n_stars) {
      for (i=1 ; i <= param.simu.n_stars ; i++) {
        read, f, val1, val2, val3, val4, val5, val6, string1, val7; param.stars(i).T = val1; param.stars(i).R = val2 ;  param.stars(i).Mass = val3;
        param.stars(i).x = val4 ; param.stars(i).y = val5 ; param.stars(i).z = val6 ;
        param.stars(i).is_bb = string1;
        read, f, param.stars(i).spectre ;
        read, f, val1, val2 ;
        param.stars(i).fUV = val1 ;
        param.stars(i).slope_fUV = val2 ;
      }
    }

  } else if (abs(param.simu.version -2.15) < 1.0e-6) {
    rdline, f;

    //photon;
    read, f, val1; param.phot.n2 = val1;
    read, f, val1; param.phot.nlambda = val1;
    read, f, val1; param.phot.nimage = val1;

    rdline, f;
    rdline, f;

    //Wavelength
    read, f, val1, val2, val3;
    param.wave.nlamb1 = val1; param.wave.lambda_min = val2; param.wave.lambda_max = val3; lambda_min=val2 ; lambda_max=val3;
    read, f, string1, string2, string3 ; param.simu.ltemp = string1; param.simu.lsed = string2; param.simu.lsed_complete = string3;
    read, f, string1; param.wave.file = string1;
    param.simu.lem_disk_image = "T";
    read, f, string1, string2 ; param.simu.lsepar = string1; param.simu.lsepar_pola = string2;
    param.map.n_type=1;
    if (param.simu.lsepar=="T") param.map.n_type += 4;
    if (param.simu.lsepar_pola=="T") param.map.n_type += 3;

    rdline, f;
    rdline, f;

    //Grid
    read, f, val1; param.grid.geometry = val1;
    read, f, val1, val2, val3, val4;
    param.grid.nrad=val1; param.grid.nz=val2; param.grid.n_az=val3 ; param.grid.nrad_in=val4;

    rdline, f;
    rdline, f;

    //map;
    read, f, val1, val2, val3, val4;
    param.map.nx = val1; param.map.ny=val2; param.map.map_size=val3; param.map.zoom = val4 ;
    param.map.size_neb = 0.5 * param.map.map_size ;
    read, f, val1, val2, val3;
    param.map.nthet = val1;   param.map.nphi = val2  ; param.map.inc = val3;
    param.map.delta = 1 ; param.map.angle_interest = 75.; param.map.lonly_bin_interest="F"
    read, f, val1, val2, val3, string1;
    param.map.RT_imin = val1 ; param.map.RT_imax = val2 ; param.map.RT_ntheta = val3 ; param.map.lRT_i_centered = string1 ;
    read, f, val1 ; param.map.dist = val1;
    read, f, val1 ; param.map.PA = val1;

    rdline, f;
    rdline, f;

    // Scattering method
    read, f, val1 ; param.simu.scatt_method = val1;
    read, f, val1 ; param.simu.aniso_method = val1;

    rdline, f;
    rdline, f;

    // Symetries
    read, f, string1; param.simu.lsym_ima = string1 ;
    read, f, string1; param.simu.lsym_centrale = string1 ;
    read, f, string1; param.simu.lsym_axiale = string1 ;

    rdline, f;
    rdline, f;

    // Dust global properties
    read, f, string1, val1, val2, com2; param.dust.lstrat = string1 ;  param.dust.strat = val1; param.dust.a_strat = val2;
    if (param.dust.lstrat=="T") {
      param.dust.strat_type=1 ;
    } else {
      param.dust.strat_type=0 ;
    }
    read, f, string1, val1; param.simu.ldust_sublimation = string1 ; param.simu.correct_Rsub = val1 ;
    read, f, string1 , val1 ; param.simu.lchauff_int = string1 ; param.simu.viscosity = val1;

    read, f;
    read, f;

    //Nbre de zones
    read,f, val1 ; param.simu.n_zones = val1;

    // On arrete ici si on n'a pas le bon nombre de zones et on recharche la structure
    if (n_zones != param.simu.n_zones) return param

    read, f;
    read, f;

    //Disk
    for (i=1 ; i <=  param.simu.n_zones ; i++) {
      read, f, val1 ; param.zones(i).geometry = val1 ;
      read, f, val1, val2; param.zones(i).mass = val1; param.zones(i).gas_to_dust = val2;
      read, f, val1, val2; param.zones(i).Ho = val1; param.zones(i).Ro = val2;
      read, f, val1, val2, val3, val4;  param.zones(i).rin = val1; param.zones(i).rout = val2; param.zones(i).edge = val3;
      read, f, val1; param.zones(i).beta = val1;
      read, f, val1; param.zones(i).dens = val1;
    }
    rdline, f;
    rdline, f;

    // Cavity
    read, f, string1 ;  param.cavity.lcavity = string1;
    read, f, val1 ;  param.cavity.H0 = val1; param.cavity.R0 = val2 ;
    read, f, val1 ;  param.cavity.beta = val1;
    rdline, f;
    rdline, f;

    //Grain
    param.simu.n_especes_max = 0 ;
    param.simu.n_components_max = 0 ;

    for (i=1 ; i <=  param.simu.n_zones ; i++) {
      read, f, val1; param.zones(i).n_especes=val1;
      if (val1 > param.simu.n_especes_max) param.simu.n_especes_max = val1 ;

      for (j=1 ; j<= param.zones(i).n_especes ; j++) {
        read, f, val1, val2, val3, val4 ;
        if (param.zones(i).n_especes <= n_especes) {
          param.dust_pop(i,j).type = "Mie" ;
          param.dust_pop(i,j).n_components = val1 ;  if (val1 > param.simu.n_components_max) param.simu.n_components_max = val1 ;
          param.dust_pop(i,j).mixing_rule = val2 ;
          param.dust_pop(i,j).porosity = val3;
          param.dust_pop(i,j).mass_fraction = val4;
        }

        for (k=1 ; k<=param.dust_pop(i,j).n_components ; k++) {
          read, f, string1, val1 ;
          if (param.dust_pop(i,j).n_components <= n_components) {
            param.dust_pop(i,j).file(k) = string1;  param.dust_pop(i,j).component_volume_fraction(k) = val1 ;
          }
        }
        read, f, val1; param.dust_pop(i,j).type_grain = val1;
        read, f, val1, val2, val3, val4 ; param.dust_pop(i,j).amin=val1;  param.dust_pop(i,j).amax=val2;  param.dust_pop(i,j).aexp=val3; param.dust_pop(i,j).n_grains=val4;
      }
    }


    if (param.dust_pop(1,1).type_grain == 2) {
      param.gridT.n_T_grain = sum(param.dust_pop.n_grains);
    } else {
      param.gridT.n_T_grain = 1;
    }
    rdline, f;
    rdline, f;


    // mol
    read, f, string1, string2, string3, val1; param.mol.lpop=string1 ; param.mol.lpop_precise=string2 ; param.mol.lLTE=string3 ; param.mol.profile_width = val1;

    read, f, val1 ; param.mol.vturb=val1 ;
    read, f, val1 ; param.simu.nmol=val1 ;

    if (param.simu.nmol == nmol)  {
      for (i=1 ; i <=  param.simu.nmol ; i++) {
        read, f, string1, val1; param.mol(i).molecule_file=string1 ; param.mol(i).level_max = val1;
        read, f, val1, val2 ;  param.mol(i).vmax=val1 ;  param.mol(i).n_speed = val2;
        read, f, string1, val1, string2 ; param.mol(i).cst_abundance=string1 ; param.mol(i).abundance=val1 ; param.mol(i).abundance_file = string2;
        read, f, string1, val1; param.mol(i).ray_tracing=string1 ; param.mol(i).n_lines_ray_tracing=val1;
        if (param.mol(i).n_lines_ray_tracing == 1) {
          read, f, val1 ; param.mol(i).transition_numbers(1) = val1;
        } else if (param.mol(i).n_lines_ray_tracing == 2) {
          read, f, val1, val2 ; param.mol(i).transition_numbers(1:2) = [val1,val2];
        } else if (param.mol(i).n_lines_ray_tracing == 3) {
          read, f, val1, val2, val3 ; param.mol(i).transition_numbers(1:3) = [val1,val2,val3];
        } else if (param.mol(i).n_lines_ray_tracing == 4) {
          read, f, val1, val2, val3, val4 ; param.mol(i).transition_numbers(1:4) = [val1,val2,val3,val4];
        } else if (param.mol(i).n_lines_ray_tracing == 5) {
          read, f, val1, val2, val3, val4, val5 ; param.mol(i).transition_numbers(1:5) = [val1,val2,val3,val4,val5];
        }
      }
    }
    rdline, f;
    rdline, f;


    // Star
    read, f, val1; param.simu.n_stars = val1;
    if (param.simu.n_stars == n_stars) {
      for (i=1 ; i <= param.simu.n_stars ; i++) {
        read, f, val1, val2, val3, val4, val5, val6, string1, val7; param.stars(i).T = val1; param.stars(i).R = val2 ;  param.stars(i).Mass = val3;
        param.stars(i).x = val4 ; param.stars(i).y = val5 ; param.stars(i).z = val6 ;
        param.stars(i).is_bb = string1;
        read, f, param.stars(i).spectre ;
        read, f, val1, val2 ;
        param.stars(i).fUV = val1 ;
        param.stars(i).slope_fUV = val2 ;
      }
    }


  } else if (abs(param.simu.version -2.13) < 1.0e-6) {
    rdline, f;

    //photon;
    read, f, val1, com1; param.phot.n1 = val1;
    read, f, val1, com1; param.phot.n2 = val1;
    read, f, val1, com1; param.phot.nlambda = val1;
    read, f, val1, com1; param.phot.nimage = val1;

    rdline, f;
    rdline, f;

    //Wavelength
    read, f, val1, val2, val3, com1;
    param.wave.nlamb1 = val1; param.wave.lambda_min = val2; param.wave.lambda_max = val3; lambda_min=val2 ; lambda_max=val3;
    read, f, string1, string2, string3 ; param.simu.ltemp = string1; param.simu.lsed = string2; param.simu.lsed_complete = string3;
    read, f, com1, com2; param.wave.file = com1;
    read, f, string1 ; param.simu.lem_disk_image = string1;
    read, f, string1, string2 ; param.simu.lsepar = string1; param.simu.lsepar_pola = string2;
    param.map.n_type=1;
    if (param.simu.lsepar=="T") param.map.n_type += 4;
    if (param.simu.lsepar_pola=="T") param.map.n_type += 3;

    rdline, f;
    rdline, f;

    //Grid
    read, f, val1; param.grid.geometry = val1;
    read, f, val1, val2, val3, val4, com1;
    param.grid.nrad=val1; param.grid.nz=val2; param.grid.n_az=val3 ; param.grid.nrad_in=val4;

    rdline, f;
    rdline, f;

    //map;
    read, f, val1, val2, val3, val4, val5, com1;
    param.map.nthet = val1;   param.map.nphi = val2
    param.map.nx = val3; param.map.ny=val4; param.map.zoom = val5;
    read, f, val1, val2, val3, string1;
    param.map.inc = val1; param.map.delta = val2; param.map.angle_interest =val3;param.map.lonly_bin_interest = string1;
    if (param.map.lonly_bin_interest=="T") {
      capt_inf = max(1,param.map.inc-param.map.delta);
      capt_sup = min(param.map.nthet,param.map.inc+param.map.delta);
      param.map.nthet = capt_sup-capt_inf+1;
    }
    read, f, val1, val2, val3, string1;
    param.map.RT_imin = val1 ; param.map.RT_imax = val2 ; param.map.RT_ntheta = val3 ; param.map.lRT_i_centered = string1 ;
    read, f, val1 ; param.map.dist = val1;
    read, f, val1 ; param.map.PA = val1;

    rdline, f;
    rdline, f;

    // Scattering method
    read, f, val1 ; param.simu.scatt_method = val1;
    read, f, val1 ; param.simu.aniso_method = val1;

    rdline, f;
    rdline, f;

    // Symetries
    read, f, string1; param.simu.lsym_ima = string1 ;
    read, f, string1; param.simu.lsym_centrale = string1 ;
    read, f, string1; param.simu.lsym_axiale = string1 ;

    rdline, f;
    rdline, f;

    // Dust global properties
    read, f, val1 ; param.dust.gas_dust = val1;
    read, f, string1, val1, val2, com2; param.dust.lstrat = string1 ;  param.dust.strat = val1; param.dust.a_strat = val2;
    if (param.dust.lstrat=="T") {
      param.dust.strat_type=1 ;
    } else {
      param.dust.strat_type=0 ;
    }
    read, f, string1, val1; param.simu.ldust_sublimation = string1 ; param.simu.correct_Rsub = val1 ;
    read, f, string1 , val1 ; param.simu.lchauff_int = string1 ; param.simu.viscosity = val1;
    read, f, val1, val2, val3; param.gridT.T_min=val1; param.gridT.T_max=val2; param.gridT.n_T=val3;

    read, f;
    read, f;

    //Nbre de zones
    read,f, val1 ; param.simu.n_zones = val1;

    // On arrete ici si on n'a pas le bon nombre de zones et on recharge la structure
    if (n_zones != param.simu.n_zones) return param

    read, f;
    read, f;

    //Disk
    for (i=1 ; i <=  param.simu.n_zones ; i++) {
      read, f, val1 ; param.zones(i).geometry = val1 ;
      read, f, val1; param.zones(i).mass = val1;
      read, f, val1, val2; param.zones(i).Ho = val1; param.zones(i).Ro = val2;
      read, f, val1, val2, val3, val4;  param.zones(i).rin = val1; param.zones(i).rout = val2; param.zones(i).size_neb = val3 ; param.zones(i).edge = val4;
      if (param.zones(i).size_neb > param.map.size_neb) param.map.size_neb = param.zones(i).size_neb ;
      read, f, val1; param.zones(i).beta = val1;
      read, f, val1; param.zones(i).dens = val1;
    }
    rdline, f;
    rdline, f;

    // Cavity
    read, f, string1 ;  param.cavity.lcavity = string1;
    read, f, val1 ;  param.cavity.H0 = val1; param.cavity.R0 = val2 ;
    read, f, val1 ;  param.cavity.beta = val1;
    rdline, f;
    rdline, f;

    //Grain
    param.simu.n_especes_max = 0 ;
    param.simu.n_components_max = 0 ;

    for (i=1 ; i <=  param.simu.n_zones ; i++) {
      read, f, val1, com1; param.zones(i).n_especes=val1;
      if (val1 > param.simu.n_especes_max) param.simu.n_especes_max = val1 ;

      for (j=1 ; j<= param.zones(i).n_especes ; j++) {
        read, f, val1, val2, val3, val4 ;
        if (param.zones(i).n_especes <= n_especes) {
          param.dust_pop(i,j).type = "Mie" ;
          param.dust_pop(i,j).n_components = val1 ;
          param.dust_pop(i,j).mixing_rule = val2 ;
          param.dust_pop(i,j).porosity = val3;
          param.dust_pop(i,j).mass_fraction = val4;
        }
        for (k=1 ; k<=param.dust_pop(i,j).n_components ; k++) {
          read, f, string1, val1 ;
          if (val1 > param.simu.n_components_max) param.simu.n_components_max = val1 ;
          if (param.dust_pop(i,j).n_components <= n_components) {
            param.dust_pop(i,j).file(k) = string1;  param.dust_pop(i,j).component_volume_fraction(k) = val1 ;
          }
        }
        read, f, val1; param.dust_pop(i,j).type_grain = val1;
        read, f, val1, val2, val3, val4 ; param.dust_pop(i,j).amin=val1;  param.dust_pop(i,j).amax=val2;  param.dust_pop(i,j).aexp=val3; param.dust_pop(i,j).n_grains=val4;
      }
    }


    if (param.dust_pop(1,1).type_grain == 2) {
      param.gridT.n_T_grain = sum(param.dust_pop.n_grains);
    } else {
      param.gridT.n_T_grain = 1;
    }
    rdline, f;
    rdline, f;


    // mol
    read, f, string1, string2, string3, val1; param.mol.lpop=string1 ; param.mol.lpop_precise=string2 ; param.mol.lLTE=string3 ; param.mol.profile_width = val1;

    read, f, val1 ; param.mol.vturb=val1 ;
    read, f, val1 ; param.simu.nmol=val1 ;

    if (param.simu.nmol == nmol)  {
      for (i=1 ; i <=  param.simu.nmol ; i++) {
        read, f, string1, val1; param.mol(i).molecule_file=string1 ; param.mol(i).level_max = val1;
        read, f, val1, val2 ;  param.mol(i).vmax=val1 ;  param.mol(i).n_speed = val2;
        read, f, string1, val1, string2 ; param.mol(i).cst_abundance=string1 ; param.mol(i).abundance=val1 ; param.mol(i).abundance_file = string2;
        read, f, string1, val1; param.mol(i).ray_tracing=string1 ; param.mol(i).n_lines_ray_tracing=val1;
        if (param.mol(i).n_lines_ray_tracing == 1) {
          read, f, val1 ; param.mol(i).transition_numbers(1) = val1;
        } else if (param.mol(i).n_lines_ray_tracing == 2) {
          read, f, val1, val2 ; param.mol(i).transition_numbers(1:2) = [val1,val2];
        } else if (param.mol(i).n_lines_ray_tracing == 3) {
          read, f, val1, val2, val3 ; param.mol(i).transition_numbers(1:3) = [val1,val2,val3];
        } else if (param.mol(i).n_lines_ray_tracing == 4) {
          read, f, val1, val2, val3, val4 ; param.mol(i).transition_numbers(1:4) = [val1,val2,val3,val4];
        } else if (param.mol(i).n_lines_ray_tracing == 5) {
          read, f, val1, val2, val3, val4, val5 ; param.mol(i).transition_numbers(1:5) = [val1,val2,val3,val4,val5];
        }
      }
    }
    rdline, f;
    rdline, f;


    // Star
    read, f, val1, com1; param.simu.n_stars = val1;
    if (param.simu.n_stars == n_stars) {
      for (i=1 ; i <= param.simu.n_stars ; i++) {
        read, f, val1, val2, val3, val4, val5, val6, string1, val7; param.stars(i).T = val1; param.stars(i).R = val2 ;  param.stars(i).Mass = val3;
        param.stars(i).x = val4 ; param.stars(i).y = val5 ; param.stars(i).z = val6 ;
        param.stars(i).is_bb = string1;
        read, f, param.stars(i).spectre ;
        read, f, val1, val2 ;
        param.stars(i).fUV = val1 ;
        param.stars(i).slope_fUV = val2 ;
      }
    }
  } else if (abs(param.simu.version -2.12) < 1.0e-6) {
    rdline, f;

    //photon;
    read, f, val1, com1; param.phot.n1 = val1;
    read, f, val1, com1; param.phot.n2 = val1;
    read, f, val1, com1; param.phot.nlambda = val1;
    read, f, val1, com1; param.phot.nimage = val1;

    rdline, f;
    rdline, f;

    //Wavelength
    read, f, val1, val2, val3, com1;
    param.wave.nlamb1 = val1; param.wave.lambda_min = val2; param.wave.lambda_max = val3; lambda_min=val2 ; lambda_max=val3;
    read, f, string1, string2, string3 ; param.simu.ltemp = string1; param.simu.lsed = string2; param.simu.lsed_complete = string3;
    read, f, com1, com2; param.wave.file = com1;
    read, f, string1 ; param.simu.lem_disk_image = string1;
    read, f, string1, string2 ; param.simu.lsepar = string1; param.simu.lsepar_pola = string2;
    param.map.n_type=1;
    if (param.simu.lsepar=="T") param.map.n_type += 4;
    if (param.simu.lsepar_pola=="T") param.map.n_type += 3;

    rdline, f;
    rdline, f;

    //Grid
    read, f, val1; param.grid.geometry = val1;
    read, f, val1, val2, val3, val4, com1;
    param.grid.nrad=val1; param.grid.nz=val2; param.grid.n_az=val3 ; param.grid.nrad_in=val4;

    rdline, f;
    rdline, f;

    //map;
    read, f, val1, val2, val3, val4, val5, com1;
    param.map.nthet = val1;   param.map.nphi = val2
    param.map.nx = val3; param.map.ny=val4; param.map.zoom = val5;
    read, f, val1, val2, val3, string1;
    param.map.inc = val1; param.map.delta = val2; param.map.angle_interest =val3;param.map.lonly_bin_interest = string1;
    if (param.map.lonly_bin_interest=="T") {
      capt_inf = max(1,param.map.inc-param.map.delta);
      capt_sup = min(param.map.nthet,param.map.inc+param.map.delta);
      param.map.nthet = capt_sup-capt_inf+1;
    }
    read, f, val1, val2, val3, string1;
    param.map.RT_imin = val1 ; param.map.RT_imax = val2 ; param.map.RT_ntheta = val3 ; param.map.lRT_i_centered = string1 ;
    read, f, val1 ; param.map.dist = val1;
    read, f, val1 ; param.map.PA = val1;

    rdline, f;
    rdline, f;

    // Scattering method
    read, f, val1 ; param.simu.scatt_method = val1;
    read, f, val1 ; param.simu.aniso_method = val1;

    rdline, f;
    rdline, f;

    // Symetries
    read, f, string1; param.simu.lsym_ima = string1 ;
    read, f, string1; param.simu.lsym_centrale = string1 ;
    read, f, string1; param.simu.lsym_axiale = string1 ;

    rdline, f;
    rdline, f;

    // Dust global properties
    read, f, val1 ; param.dust.gas_dust = val1;
    read, f, string1, val1, val2, com2; param.dust.lstrat = string1 ;  param.dust.strat = val1; param.dust.a_strat = val2;
    if (param.dust.lstrat=="T") {
      param.dust.strat_type=1 ;
    } else {
      param.dust.strat_type=0 ;
    }
    read, f, string1, val1; param.simu.ldust_sublimation = string1 ; param.simu.correct_Rsub = val1 ;
    read, f, string1 , val1 ; param.simu.lchauff_int = string1 ; param.simu.viscosity = val1;
    read, f, val1, val2, val3; param.gridT.T_min=val1; param.gridT.T_max=val2; param.gridT.n_T=val3;

    read, f;
    read, f;

    //Nbre de zones
    read,f, val1 ; param.simu.n_zones = val1;

    read, f;
    read, f;

    if (n_zones!=param.simu.n_zones) { // on veut juste lire n_especes
      for (i=1 ; i <=  param.simu.n_zones ; i++) {
        rdline, f;
        rdline, f;
        rdline, f;
        rdline, f;
        rdline, f;
        rdline, f;

        rdline, f;
        rdline, f;

        rdline, f;
        rdline, f;
        rdline, f;

        rdline, f;
        rdline, f;

        read, f, val1, com1;
        for (j=1 ; j<= val1 ; j++) {
          rdline, f;
          rdline, f;
          rdline, f;
        }

        rdline, f;
        rdline, f;
      }
    } else {
      //Disk
      for (i=1 ; i <=  param.simu.n_zones ; i++) {
        read, f, val1 ; param.zones(i).geometry = val1 ;
        read, f, val1; param.zones(i).mass = val1;
        read, f, val1, val2; param.zones(i).Ho = val1; param.zones(i).Ro = val2;
        read, f, val1, val2, val3, val4;  param.zones(i).rin = val1; param.zones(i).rout = val2; param.zones(i).size_neb = val3 ; param.zones(i).edge = val4;
        if (param.zones(i).size_neb > param.map.size_neb) {
          param.map.size_neb = param.zones(i).size_neb ;
          param.map.map_size = 2*param.map.size_neb ;
        }
        read, f, val1; param.zones(i).beta = val1;
        read, f, val1; param.zones(i).dens = val1;
      }
      rdline, f;
      rdline, f;

      // Cavity
      read, f, string1 ;  param.cavity.lcavity = string1;
      read, f, val1 ;  param.cavity.H0 = val1; param.cavity.R0 = val2 ;
      read, f, val1 ;  param.cavity.beta = val1;
      rdline, f;
      rdline, f;

      //Grain
      for (i=1 ; i <=  param.simu.n_zones ; i++) {
        read, f, val1, com1; param.zones(i).n_especes=val1;
        if (param.zones(i).n_especes > n_especes) { // on ne lit pas les valeurs
          "**************************************************" ;
          "ERROR : reading parameter files : pb with n_espece" ;
          write, param.zones(i).n_especes, n_especes ;
          "**************************************************" ;

          for (j=1 ; j<= param.zones(i).n_especes ; j++) {
            read, f, string1, val1, val2;
            read, f, val1;
            read, f, val1, val2, val3, val4 ;
          }
        } else { // On lit les valeurs
          for (j=1 ; j<= param.zones(i).n_especes ; j++) {
            read, f, string1, val1, val2; param.dust_pop(i,j).file = string1;  param.dust_pop(i,j).porosity = val1; param.dust_pop(i,j).mass_fraction = val2;
            read, f, val1; param.dust_pop(i,j).type_grain = val1;
            read, f, val1, val2, val3, val4 ; param.dust_pop(i,j).amin=val1;  param.dust_pop(i,j).amax=val2;  param.dust_pop(i,j).aexp=val3; param.dust_pop(i,j).n_grains=val4;
          }
        } // test n_especes
      }


      if (param.dust_pop(1,1).type_grain == 2) {
        param.gridT.n_T_grain = sum(param.dust_pop.n_grains);
      } else {
        param.gridT.n_T_grain = 1;
      }
      rdline, f;
      rdline, f;

    }

    // mol
    read, f, string1, string2, string3, val1; param.mol.lpop=string1 ; param.mol.lpop_precise=string2 ; param.mol.lLTE=string3 ; param.mol.profile_width = val1;

    read, f, val1 ; param.mol.vturb=val1 ;
    read, f, val1 ; param.simu.nmol=val1 ;

    if (param.simu.nmol == nmol)  {
      for (i=1 ; i <=  param.simu.nmol ; i++) {
        read, f, string1, val1; param.mol(i).molecule_file=string1 ; param.mol(i).level_max = val1;
        read, f, val1, val2 ;  param.mol(i).vmax=val1 ;  param.mol(i).n_speed = val2;
        read, f, string1, val1, string2 ; param.mol(i).cst_abundance=string1 ; param.mol(i).abundance=val1 ; param.mol(i).abundance_file = string2;
        read, f, string1, val1; param.mol(i).ray_tracing=string1 ; param.mol(i).n_lines_ray_tracing=val1;
        if (param.mol(i).n_lines_ray_tracing == 1) {
          read, f, val1 ; param.mol(i).transition_numbers(1) = val1;
        } else if (param.mol(i).n_lines_ray_tracing == 2) {
          read, f, val1, val2 ; param.mol(i).transition_numbers(1:2) = [val1,val2];
        } else if (param.mol(i).n_lines_ray_tracing == 3) {
          read, f, val1, val2, val3 ; param.mol(i).transition_numbers(1:3) = [val1,val2,val3];
        } else if (param.mol(i).n_lines_ray_tracing == 4) {
          read, f, val1, val2, val3, val4 ; param.mol(i).transition_numbers(1:4) = [val1,val2,val3,val4];
        } else if (param.mol(i).n_lines_ray_tracing == 5) {
          read, f, val1, val2, val3, val4, val5 ; param.mol(i).transition_numbers(1:5) = [val1,val2,val3,val4,val5];
        }
      }
    }
    rdline, f;
    rdline, f;


    // Star
    read, f, val1, com1; param.simu.n_stars = val1;
    if (param.simu.n_stars == n_stars) {
      for (i=1 ; i <= param.simu.n_stars ; i++) {
        read, f, val1, val2, val3, val4, val5, val6, string1, val7; param.stars(i).T = val1; param.stars(i).R = val2 ;  param.stars(i).Mass = val3;
        param.stars(i).x = val4 ; param.stars(i).y = val5 ; param.stars(i).z = val6 ;
        param.stars(i).is_bb = string1;
        read, f, param.stars(i).spectre ;
        read, f, val1, val2 ;
        param.stars(i).fUV = val1 ;
        param.stars(i).slope_fUV = val2 ;
      }
    }
  } else  if (abs(param.simu.version -2.11) < 1.0e-6) {
    rdline, f;

    //photon;
    read, f, val1, com1; param.phot.n1 = val1;
    read, f, val1, com1; param.phot.n2 = val1;
    read, f, val1, com1; param.phot.nlambda = val1;
    read, f, val1, com1; param.phot.nimage = val1;

    rdline, f;
    rdline, f;

    //Wavelength
    read, f, val1, val2, val3, com1;
    param.wave.nlamb1 = val1; param.wave.lambda_min = val2; param.wave.lambda_max = val3; lambda_min=val2 ; lambda_max=val3;
    read, f, string1, string2, string3 ; param.simu.ltemp = string1; param.simu.lsed = string2; param.simu.lsed_complete = string3;
    read, f, com1, com2; param.wave.file = com1;
    read, f, string1 ; param.simu.lem_disk_image = string1;
    read, f, string1, string2 ; param.simu.lsepar = string1; param.simu.lsepar_pola = string2;
    param.map.n_type=1;
    if (param.simu.lsepar=="T") param.map.n_type += 4;
    if (param.simu.lsepar_pola=="T") param.map.n_type += 3;

    rdline, f;
    rdline, f;

    //Grid
    read, f, val1; param.grid.geometry = val1;
    read, f, val1, val2, val3, val4, com1;
    param.grid.nrad=val1; param.grid.nz=val2; param.grid.n_az=val3 ; param.grid.nrad_in=val4;

    rdline, f;
    rdline, f;

    //map;
    read, f, val1, val2, val3, val4, val5, com1;
    param.map.nthet = val1;   param.map.nphi = val2
    param.map.nx = val3; param.map.ny=val4; param.map.zoom = val5;
    read, f, val1, val2, val3, string1;
    param.map.inc = val1; param.map.delta = val2; param.map.angle_interest =val3;param.map.lonly_bin_interest = string1;
    if (param.map.lonly_bin_interest=="T") {
      capt_inf = max(1,param.map.inc-param.map.delta);
      capt_sup = min(param.map.nthet,param.map.inc+param.map.delta);
      param.map.nthet = capt_sup-capt_inf+1;
    }
    read, f, val1, val2, val3, string1;
    param.map.RT_imin = val1 ; param.map.RT_imax = val2 ; param.map.RT_ntheta = val3 ; param.map.lRT_i_centered = string1 ;
    read, f, val1 ; param.map.dist = val1;
    read, f, val1 ; param.map.PA = val1;

    rdline, f;
    rdline, f;

    // Scattering method
    read, f, val1 ; param.simu.scatt_method = val1;
    read, f, val1 ; param.simu.aniso_method = val1;

    rdline, f;
    rdline, f;

    // Symetries
    read, f, string1; param.simu.lsym_ima = string1 ;
    read, f, string1; param.simu.lsym_centrale = string1 ;
    read, f, string1; param.simu.lsym_axiale = string1 ;

    rdline, f;
    rdline, f;

    // Dust global properties
    read, f, val1 ; param.dust.gas_dust = val1;
    read, f, string1, val1, val2, com2; param.dust.lstrat = string1 ;  param.dust.strat = val1; param.dust.a_strat = val2;
    if (param.dust.lstrat=="T") {
      param.dust.strat_type=1 ;
    } else {
      param.dust.strat_type=0 ;
    }
    read, f, string1, val1; param.simu.ldust_sublimation = string1 ; param.simu.correct_Rsub = val1 ;
    read, f, string1 , val1 ; param.simu.lchauff_int = string1 ; param.simu.viscosity = val1;
    read, f, val1, val2, val3; param.gridT.T_min=val1; param.gridT.T_max=val2; param.gridT.n_T=val3;

    read, f;
    read, f;

    //Nbre de zones
    read,f, val1 ; param.simu.n_zones = val1;

    read, f;
    read, f;

    if (n_zones!=param.simu.n_zones) { // on veut juste lire n_especes
      for (i=1 ; i <=  param.simu.n_zones ; i++) {
        rdline, f;
        rdline, f;
        rdline, f;
        rdline, f;
        rdline, f;
        rdline, f;

        rdline, f;
        rdline, f;

        rdline, f;
        rdline, f;
        rdline, f;

        rdline, f;
        rdline, f;

        read, f, val1, com1;
        for (j=1 ; j<= val1 ; j++) {
          rdline, f;
          rdline, f;
          rdline, f;
        }

        rdline, f;
        rdline, f;
      }
    } else {
      //Disk
      for (i=1 ; i <=  param.simu.n_zones ; i++) {
        read, f, val1 ; param.zones(i).geometry = val1 ;
        read, f, val1; param.zones(i).mass = val1;
        read, f, val1, val2; param.zones(i).Ho = val1; param.zones(i).Ro = val2;
        read, f, val1, val2, val3, val4;  param.zones(i).rin = val1; param.zones(i).rout = val2; param.zones(i).size_neb = val3 ; param.zones(i).edge = val4;
        if (param.zones(i).size_neb > param.map.size_neb) param.map.size_neb = param.zones(i).size_neb ;
        read, f, val1; param.zones(i).beta = val1;
        read, f, val1; param.zones(i).dens = val1;
      }
      rdline, f;
      rdline, f;

      // Cavity
      read, f, string1 ;  param.cavity.lcavity = string1;
      read, f, val1 ;  param.cavity.H0 = val1; param.cavity.R0 = val2 ;
      read, f, val1 ;  param.cavity.beta = val1;
      rdline, f;
      rdline, f;

      //Grain
      for (i=1 ; i <=  param.simu.n_zones ; i++) {
        read, f, val1, com1; param.zones(i).n_especes=val1;
        if (param.zones(i).n_especes > n_especes) { // on ne lit pas les valeurs
          for (j=1 ; j<= param.zones(i).n_especes ; j++) {
            read, f, string1, val1, val2;
            read, f, val1;
            read, f, val1, val2, val3, val4 ;
          }
        } else { // On lit les valeurs
          for (j=1 ; j<= param.zones(i).n_especes ; j++) {
            read, f, string1, val1, val2; param.dust_pop(i,j).file = string1;  param.dust_pop(i,j).porosity = val1; param.dust_pop(i,j).mass_fraction = val2;
            read, f, val1; param.dust_pop(i,j).type_grain = val1;
            read, f, val1, val2, val3, val4 ; param.dust_pop(i,j).amin=val1;  param.dust_pop(i,j).amax=val2;  param.dust_pop(i,j).aexp=val3; param.dust_pop(i,j).n_grains=val4;
          }
        } // test n_especes
      }


      if (param.dust_pop(1,1).type_grain == 2) {
        param.gridT.n_T_grain = sum(param.dust_pop.n_grains);
      } else {
        param.gridT.n_T_grain = 1;
      }
      rdline, f;
      rdline, f;

    }

    // mol
    read, f, string1, string2, string3, val1; param.mol.lpop=string1 ; param.mol.lpop_precise=string2 ; param.mol.lLTE=string3 ; param.mol.profile_width = val1;

    read, f, val1 ; param.mol.vturb=val1 ;
    read, f, val1 ; param.simu.nmol=val1 ;

    if (param.simu.nmol == nmol)  {
      for (i=1 ; i <=  param.simu.nmol ; i++) {
        read, f, string1, val1; param.mol(i).molecule_file=string1 ; param.mol(i).level_max = val1;
        read, f, val1, val2 ;  param.mol(i).vmax=val1 ;  param.mol(i).n_speed = val2;
        read, f, string1, val1, string2 ; param.mol(i).cst_abundance=string1 ; param.mol(i).abundance=val1 ; param.mol(i).abundance_file = string2;
        read, f, string1, val1; param.mol(i).ray_tracing=string1 ; param.mol(i).n_lines_ray_tracing=val1;
        if (param.mol(i).n_lines_ray_tracing == 1) {
          read, f, val1 ; param.mol(i).transition_numbers(1) = val1;
        } else if (param.mol(i).n_lines_ray_tracing == 2) {
          read, f, val1, val2 ; param.mol(i).transition_numbers(1:2) = [val1,val2];
        } else if (param.mol(i).n_lines_ray_tracing == 3) {
          read, f, val1, val2, val3 ; param.mol(i).transition_numbers(1:3) = [val1,val2,val3];
        } else if (param.mol(i).n_lines_ray_tracing == 4) {
          read, f, val1, val2, val3, val4 ; param.mol(i).transition_numbers(1:4) = [val1,val2,val3,val4];
        } else if (param.mol(i).n_lines_ray_tracing == 5) {
          read, f, val1, val2, val3, val4, val5 ; param.mol(i).transition_numbers(1:5) = [val1,val2,val3,val4,val5];
        }
      }
    }
    rdline, f;
    rdline, f;


    // Star
    read, f, val1, com1; param.simu.n_stars = val1;
    if (param.simu.n_stars == n_stars) {
      for (i=1 ; i <= param.simu.n_stars ; i++) {
        read, f, val1, val2, val3, val4, val5, val6, string1, val7; param.stars(i).T = val1; param.stars(i).R = val2 ;  param.stars(i).Mass = val3;
        param.stars(i).x = val4 ; param.stars(i).y = val5 ; param.stars(i).z = val6 ;
        param.stars(i).is_bb = string1;
        read, f, param.stars(i).spectre ;
        param.stars(i).fUV = val7 ;
      }
    }
  } else if (abs(param.simu.version -2.10) < 1.0e-6) {
    rdline, f;

    //photon;
    read, f, val1, com1; param.phot.n1 = val1;
    read, f, val1, com1; param.phot.n2 = val1;
    read, f, val1, com1; param.phot.nlambda = val1;
    read, f, val1, com1; param.phot.nimage = val1;
    read, f, string1, val1; param.simu.lcheckpoint = string1 ; param.simu.checkpoint_time = val1;

    rdline, f;
    rdline, f;

    //Wavelength
    read, f, val1, val2, val3, com1;
    param.wave.nlamb1 = val1; param.wave.lambda_min = val2; param.wave.lambda_max = val3; lambda_min=val2 ; lambda_max=val3;
    read, f, string1, string2, string3 ; param.simu.ltemp = string1; param.simu.lsed = string2; param.simu.lsed_complete = string3;
    read, f, com1, com2; param.wave.file = com1;
    read, f, string1 ; param.simu.lem_disk_image = string1;
    read, f, string1, string2 ; param.simu.lsepar = string1; param.simu.lsepar_pola = string2;
    param.map.n_type=1;
    if (param.simu.lsepar=="T") param.map.n_type += 4;
    if (param.simu.lsepar_pola=="T") param.map.n_type += 3;
    read, f, val1, val2 ; param.simu.tau_seuil = val1 ; param.simu.wl_seuil = val2;

    rdline, f;
    rdline, f;

    //Grid
    read, f, val1; param.grid.geometry = val1;
    read, f, val1, val2, val3, val4, com1;
    param.grid.nrad=val1; param.grid.nz=val2; param.grid.n_az=val3 ; param.grid.nrad_in=val4;

    rdline, f;
    rdline, f;

    //map;
    read, f, val1, val2, val3, val4, val5, com1;
    param.map.nthet = val1;   param.map.nphi = val2
    param.map.nx = val3; param.map.ny=val4; param.map.zoom = val5;
    read, f, val1, val2, val3, string1;
    param.map.inc = val1; param.map.delta = val2; param.map.angle_interest =val3;param.map.lonly_bin_interest = string1;
    if (param.map.lonly_bin_interest=="T") {
      capt_inf = max(1,param.map.inc-param.map.delta);
      capt_sup = min(param.map.nthet,param.map.inc+param.map.delta);
      param.map.nthet = capt_sup-capt_inf+1;
    }
    read, f, val1, val2, val3, string1;
    param.map.RT_imin = val1 ; param.map.RT_imax = val2 ; param.map.RT_ntheta = val3 ; param.map.lRT_i_centered = string1 ;
    read, f, val1 ; param.map.dist = val1;
    read, f, val1 ; param.map.PA = val1;

    rdline, f;
    rdline, f;

    // Scattering method
    read, f, val1 ; param.simu.scatt_method = val1;
    read, f, val1 ; param.simu.aniso_method = val1;

    rdline, f;
    rdline, f;

    // Symetries
    read, f, string1; param.simu.lsym_ima = string1 ;
    read, f, string1; param.simu.lsym_centrale = string1 ;
    read, f, string1; param.simu.lsym_axiale = string1 ;

    rdline, f;
    rdline, f;

    // Dust global properties
    read, f, val1 ; param.dust.gas_dust = val1;
    read, f, string1, val1, val2, com2; param.dust.lstrat = string1 ;  param.dust.strat = val1; param.dust.a_strat = val2;
    if (param.dust.lstrat=="T") {
      param.dust.strat_type=1 ;
    } else {
      param.dust.strat_type=0 ;
    }
    read, f, string1; param.simu.ldust_sublimation = string1 ;
    read, f, string1 , val1 ; param.simu.lchauff_int = string1 ; param.simu.viscosity = val1;
    read, f, val1, val2, val3; param.gridT.T_min=val1; param.gridT.T_max=val2; param.gridT.n_T=val3;

    read, f;
    read, f;

    //Nbre de zones
    read,f, val1 ; param.simu.n_zones = val1;

    read, f;
    read, f;

    if (n_zones!=param.simu.n_zones) { // on veut juste lire n_especes
      for (i=1 ; i <=  param.simu.n_zones ; i++) {
        rdline, f;
        rdline, f;
        rdline, f;
        rdline, f;
        rdline, f;
        rdline, f;

        rdline, f;
        rdline, f;

        rdline, f;
        rdline, f;
        rdline, f;

        rdline, f;
        rdline, f;

        read, f, val1, com1;
        for (j=1 ; j<= val1 ; j++) {
          rdline, f;
          rdline, f;
          rdline, f;
        }

        rdline, f;
        rdline, f;
      }
    } else {
      //Disk
      for (i=1 ; i <=  param.simu.n_zones ; i++) {
        read, f, val1 ; param.zones(i).geometry = val1 ;
        read, f, val1; param.zones(i).mass = val1;
        read, f, val1, val2; param.zones(i).Ho = val1; param.zones(i).Ro = val2;
        read, f, val1, val2, val3, val4;  param.zones(i).rin = val1; param.zones(i).rout = val2; param.zones(i).size_neb = val3 ; param.zones(i).edge = val4;
        if (param.zones(i).size_neb > param.map.size_neb) param.map.size_neb = param.zones(i).size_neb ;
        read, f, val1; param.zones(i).beta = val1;
        read, f, val1; param.zones(i).dens = val1;
      }
      rdline, f;
      rdline, f;

      // Cavity
      read, f, string1 ;  param.cavity.lcavity = string1;
      read, f, val1 ;  param.cavity.H0 = val1; param.cavity.R0 = val2 ;
      read, f, val1 ;  param.cavity.beta = val1;
      rdline, f;
      rdline, f;

      //Grain
      for (i=1 ; i <=  param.simu.n_zones ; i++) {
        read, f, val1, com1; param.zones(i).n_especes=val1;
        if (param.zones(i).n_especes > n_especes) { // on ne lit pas les valeurs
          for (j=1 ; j<= param.zones(i).n_especes ; j++) {
            read, f, string1, val1, val2;
            read, f, val1;
            read, f, val1, val2, val3, val4 ;
          }
        } else { // On lit les valeurs
          for (j=1 ; j<= param.zones(i).n_especes ; j++) {
            read, f, string1, val1, val2; param.dust_pop(i,j).file = string1;  param.dust_pop(i,j).porosity = val1; param.dust_pop(i,j).mass_fraction = val2;
            read, f, val1; param.dust_pop(i,j).type_grain = val1;
            read, f, val1, val2, val3, val4 ; param.dust_pop(i,j).amin=val1;  param.dust_pop(i,j).amax=val2;  param.dust_pop(i,j).aexp=val3; param.dust_pop(i,j).n_grains=val4;
          }
        } // test n_especes
      }


      if (param.dust_pop(1,1).type_grain == 2) {
        param.gridT.n_T_grain = sum(param.dust_pop.n_grains);
      } else {
        param.gridT.n_T_grain = 1;
      }
      rdline, f;
      rdline, f;

    }

    // mol
    read, f, val1, val2, val3 ; param.mol.vmax=val1 ; param.mol.vturb=val2 ; param.mol.n_speed = val3;
    read, f, string1, string2, string3, val1; param.mol.lpop=string1 ; param.mol.lpop_precise=string2 ; param.mol.lLTE=string3 ; param.mol.profile_width = val1;
    read, f, string1, val1; param.mol.molecule_file=string1 ; param.mol.level_max = val1;
    read, f, string1, val1, string2 ; param.mol.cst_abundance=string1 ; param.mol.abundance=val1 ; param.mol.abundance_file = string2;
    read, f, string1, val1; param.mol.ray_tracing=string1 ; param.mol.n_lines_ray_tracing=val1;


    param.simu.nmol = 1 ;
    if (param.mol(1).n_lines_ray_tracing == 1) {
      read, f, val1 ; param.mol(1).transition_numbers(1) = val1;
    } else if (param.mol(1).n_lines_ray_tracing == 2) {
      read, f, val1, val2 ; param.mol(1).transition_numbers(1:2) = [val1,val2];
    } else if (param.mol(1).n_lines_ray_tracing == 3) {
      read, f, val1, val2, val3 ; param.mol(1).transition_numbers(1:3) = [val1,val2,val3];
    } else if (param.mol(1).n_lines_ray_tracing == 4) {
      read, f, val1, val2, val3, val4 ; param.mol(1).transition_numbers(1:4) = [val1,val2,val3,val4];
    } else if (param.mol(1).n_lines_ray_tracing == 5) {
      read, f, val1, val2, val3, val4, val5 ; param.mol(1).transition_numbers(1:5) = [val1,val2,val3,val4,val5];
    }

    rdline, f;
    rdline, f;

    // Star
    read, f, val1, com1; param.simu.n_stars = val1;
    if (param.simu.n_stars == n_stars) {
      for (i=1 ; i <= param.simu.n_stars ; i++) {
        read, f, val1, val2, val3, val4, val5, val6, string1; param.stars(i).T = val1; param.stars(i).R = val2 ;  param.stars(i).Mass = val3;
        param.stars(i).x = val4 ; param.stars(i).y = val5 ; param.stars(i).z = val6 ;
        param.stars(i).is_bb = string1;
        read, f, param.stars(i).spectre ;
      }
    }
  } else if (abs(param.simu.version -2.09) < 1.0e-6) {
    rdline, f;

    //photon;
    read, f, val1, com1; param.phot.n1 = val1;
    read, f, val1, com1; param.phot.n2 = val1;
    read, f, val1, com1; param.phot.nlambda = val1;
    read, f, val1, com1; param.phot.nimage = val1;
    read, f, string1, val1; param.simu.lcheckpoint = string1 ; param.simu.checkpoint_time = val1;

    rdline, f;
    rdline, f;

    //Wavelength
    read, f, val1, val2, val3, com1;
    param.wave.nlamb1 = val1; param.wave.lambda_min = val2; param.wave.lambda_max = val3; lambda_min=val2 ; lambda_max=val3;
    read, f, string1, string2, string3 ; param.simu.ltemp = string1; param.simu.lsed = string2; param.simu.lsed_complete = string3;
    read, f, com1, com2; param.wave.file = com1;
    read, f, string1 ; param.simu.lem_disk_image = string1;
    read, f, string1, string2 ; param.simu.lsepar = string1; param.simu.lsepar_pola = string2;
    param.map.n_type=1;
    if (param.simu.lsepar=="T") param.map.n_type += 4;
    if (param.simu.lsepar_pola=="T") param.map.n_type += 3;
    read, f, val1, val2 ; param.simu.tau_seuil = val1 ; param.simu.wl_seuil = val2;

    rdline, f;
    rdline, f;

    //Grid
    read, f, val1; param.grid.geometry = val1;
    read, f, val1, val2, val3, val4, com1;
    param.grid.nrad=val1; param.grid.nz=val2; param.grid.n_az=val3 ; param.grid.nrad_in=val4;

    rdline, f;
    rdline, f;

    //map;
    read, f, val1, val2, val3, val4, val5, com1;
    param.map.nthet = val1;   param.map.nphi = val2
    param.map.nx = val3; param.map.ny=val4; param.map.zoom = val5;
    read, f, val1, val2, val3, string1; param.map.inc = val1; param.map.delta = val2; param.map.angle_interest =val3;param.map.lonly_bin_interest = string1;
    if (param.map.lonly_bin_interest=="T") {
      capt_inf = max(1,param.map.inc-param.map.delta);
      capt_sup = min(param.map.nthet,param.map.inc+param.map.delta);
      param.map.nthet = capt_sup-capt_inf+1;
    }
    read, f, val1 ; param.map.dist = val1;
    read, f, val1 ; param.map.PA = val1;


    rdline, f;
    rdline, f;

    // Scattering method
    read, f, val1 ; param.simu.scatt_method = val1;
    read, f, val1 ; param.simu.aniso_method = val1;

    rdline, f;
    rdline, f;

    // Symetries
    read, f, string1; param.simu.lsym_ima = string1 ;
    read, f, string1; param.simu.lsym_centrale = string1 ;
    read, f, string1; param.simu.lsym_axiale = string1 ;

    rdline, f;
    rdline, f;

    // Dust global properties
    read, f, val1 ; param.dust.gas_dust = val1;
    read, f, string1, val1, val2, com2; param.dust.lstrat = string1 ;  param.dust.strat = val1; param.dust.a_strat = val2;
    if (param.dust.lstrat=="T") {
      param.dust.strat_type=1 ;
    } else {
      param.dust.strat_type=0 ;
    }
    read, f, string1; param.simu.ldust_sublimation = string1 ;
    read, f, string1 , val1 ; param.simu.lchauff_int = string1 ; param.simu.viscosity = val1;
    read, f, val1, val2, val3; param.gridT.T_min=val1; param.gridT.T_max=val2; param.gridT.n_T=val3;

    read, f;
    read, f;

    //Nbre de zones
    read,f, val1 ; param.simu.n_zones = val1;

    read, f;
    read, f;

    if (n_zones!=param.simu.n_zones) { // on veut juste lire n_especes
      for (i=1 ; i <=  param.simu.n_zones ; i++) {
        rdline, f;
        rdline, f;
        rdline, f;
        rdline, f;
        rdline, f;
        rdline, f;

        rdline, f;
        rdline, f;
        read, f, val1, com1;
        for (j=1 ; j<= val1 ; j++) {
          rdline, f;
          rdline, f;
          rdline, f;
        }

        rdline, f;
        rdline, f;
      }
    } else {
      //Disk
      for (i=1 ; i <=  param.simu.n_zones ; i++) {
        read, f, val1 ; param.zones(i).geometry = val1 ;
        read, f, val1; param.zones(i).mass = val1;
        read, f, val1, val2; param.zones(i).Ho = val1; param.zones(i).Ro = val2;
        read, f, val1, val2, val3, val4;  param.zones(i).rin = val1; param.zones(i).rout = val2; param.zones(i).size_neb = val3 ; param.zones(i).edge = val4;
        if (param.zones(i).size_neb > param.map.size_neb) param.map.size_neb = param.zones(i).size_neb ;
        read, f, val1; param.zones(i).beta = val1;
        read, f, val1; param.zones(i).dens = val1;
      }
      rdline, f;
      rdline, f;


      //Grain
      for (i=1 ; i <=  param.simu.n_zones ; i++) {
        read, f, val1, com1; param.zones(i).n_especes=val1;
        for (j=1 ; j<= param.zones(i).n_especes ; j++) {
          read, f, string1, val1, val2; param.dust_pop(i,j).file = string1;  param.dust_pop(i,j).porosity = val1; param.dust_pop(i,j).mass_fraction = val2;
          read, f, val1; param.dust_pop(i,j).type_grain = val1;
          read, f, val1, val2, val3, val4 ; param.dust_pop(i,j).amin=val1;  param.dust_pop(i,j).amax=val2;  param.dust_pop(i,j).aexp=val3; param.dust_pop(i,j).n_grains=val4;
        }
      }
      if (param.dust_pop(1,1).type_grain == 2) {
        param.gridT.n_T_grain = sum(param.dust_pop.n_grains);
      } else {
        param.gridT.n_T_grain = 1;
      }
      rdline, f;
      rdline, f;

    }

    // mol
    read, f, val1, val2, val3 ; param.mol.vmax=val1 ; param.mol.vturb=val2 ; param.mol.n_speed = val3;
    read, f, string1, string2, string3, val1; param.mol.lpop=string1 ; param.mol.lpop_precise=string2 ; param.mol.lLTE=string3 ; param.mol.profile_width = val1;
    read, f, string1, val1; param.mol.molecule_file=string1 ; param.mol.level_max = val1;
    read, f, string1, val1, string2 ; param.mol.cst_abundance=string1 ; param.mol.abundance=val1 ; param.mol.abundance_file = string2;
    read, f, string1, val1; param.mol.ray_tracing=string1 ; param.mol.n_lines_ray_tracing=val1;
    if (param.mol.n_lines_ray_tracing(1) == n_lines_ray_tracing)  {
      if (param.mol.n_lines_ray_tracing(1) == 1) {
        read, f, val1 ; param.mol.transition_numbers = val1;
      } else if (param.mol.n_lines_ray_tracing(1) == 2) {
        read, f, val1, val2 ; param.mol.transition_numbers = [val1,val2];
      } else if (param.mol.n_lines_ray_tracing(1) == 3) {
        read, f, val1, val2, val3 ; param.mol.transition_numbers = [val1,val2,val3];
      } else if (param.mol.n_lines_ray_tracing(1) == 4) {
        read, f, val1, val2, val3, val4 ; param.mol.transition_numbers = [val1,val2,val3,val4];
      } else if (param.mol.n_lines_ray_tracing(1) == 5) {
        read, f, val1, val2, val3, val4, val5 ; param.mol.transition_numbers = [val1,val2,val3,val4,val5];
      } else {
        print, "0 mol" ;
      }
    }
    rdline, f;
    rdline, f;

    // Star
    read, f, val1, com1; param.simu.n_stars = val1;
    if (param.simu.n_stars == n_stars) {
      for (i=1 ; i <= param.simu.n_stars ; i++) {
        read, f, val1, val2, val3, val4, val5, val6, string1; param.stars(i).T = val1; param.stars(i).R = val2 ;  param.stars(i).Mass = val3;
        param.stars(i).x = val4 ; param.stars(i).y = val5 ; param.stars(i).z = val6 ;
        param.stars(i).is_bb = string1;
        read, f, param.stars(i).spectre ;
      }
    }
  }  else if (abs(param.simu.version -2.08) < 1.0e-6) {
    rdline, f;

    //photon;
    read, f, val1, com1; param.phot.n1 = val1;
    read, f, val1, com1; param.phot.n2 = val1;
    read, f, val1, com1; param.phot.nlambda = val1;
    read, f, val1, com1; param.phot.nimage = val1;
    read, f, string1, val1; param.simu.lcheckpoint = string1 ; param.simu.checkpoint_time = val1;

    rdline, f;
    rdline, f;

    //Wavelength
    read, f, val1, val2, val3, com1;
    param.wave.nlamb1 = val1; param.wave.lambda_min = val2; param.wave.lambda_max = val3; lambda_min=val2 ; lambda_max=val3;
    read, f, string1, string2, string3 ; param.simu.ltemp = string1; param.simu.lsed = string2; param.simu.lsed_complete = string3;
    read, f, com1, com2; param.wave.file = com1;
    read, f, string1 ; param.simu.lem_disk_image = string1;
    read, f, string1, string2 ; param.simu.lsepar = string1; param.simu.lsepar_pola = string2;
    param.map.n_type=1;
    if (param.simu.lsepar=="T") param.map.n_type += 4;
    if (param.simu.lsepar_pola=="T") param.map.n_type += 3;
    read, f, val1, val2 ; param.simu.tau_seuil = val1 ; param.simu.wl_seuil = val2;

    rdline, f;
    rdline, f;

    //Grid
    read, f, val1; param.grid.geometry = val1;
    read, f, val1, val2, val3, val4, com1;
    param.grid.nrad=val1; param.grid.nz=val2; param.grid.n_az=val3 ; param.grid.nrad_in=val4;

    rdline, f;
    rdline, f;

    //map;
    read, f, val1, val2, val3, val4, val5, com1;
    param.map.nthet = val1;   param.map.nphi = val2
    param.map.nx = val3; param.map.ny=val4; param.map.zoom = val5;
    read, f, val1, val2, val3, string1; param.map.inc = val1; param.map.delta = val2; param.map.angle_interest =val3;param.map.lonly_bin_interest = string1;
    if (param.map.lonly_bin_interest=="T") {
      capt_inf = max(1,param.map.inc-param.map.delta);
      capt_sup = min(param.map.nthet,param.map.inc+param.map.delta);
      param.map.nthet = capt_sup-capt_inf+1;
    }
    read, f, val1, com1; param.map.dist = val1;

    rdline, f;
    rdline, f;

    // Scattering method
    read, f, val1 ; param.simu.scatt_method = val1;
    read, f, val1 ; param.simu.aniso_method = val1;

    rdline, f;
    rdline, f;

    // Symetries
    read, f, string1; param.simu.lsym_ima = string1 ;
    read, f, string1; param.simu.lsym_centrale = string1 ;
    read, f, string1; param.simu.lsym_axiale = string1 ;

    rdline, f;
    rdline, f;

    // Dust global properties
    read, f, val1 ; param.dust.gas_dust = val1;
    read, f, string1, val1, val2, com2; param.dust.lstrat = string1 ;  param.dust.strat = val1; param.dust.a_strat = val2;
    if (param.dust.lstrat=="T") {
      param.dust.strat_type=1 ;
    } else {
      param.dust.strat_type=0 ;
    }
    read, f, string1; param.simu.ldust_sublimation = string1 ;
    read, f, string1 , val1 ; param.simu.lchauff_int = string1 ; param.simu.viscosity = val1;
    read, f, val1, val2, val3; param.gridT.T_min=val1; param.gridT.T_max=val2; param.gridT.n_T=val3;

    read, f;
    read, f;

    //Nbre de zones
    read,f, val1 ; param.simu.n_zones = val1;

    read, f;
    read, f;

    if (n_zones!=param.simu.n_zones) { // on veut juste lire n_especes
      for (i=1 ; i <=  param.simu.n_zones ; i++) {
        rdline, f;
        rdline, f;
        rdline, f;
        rdline, f;
        rdline, f;
        rdline, f;

        rdline, f;
        rdline, f;
        read, f, val1, com1;
        for (j=1 ; j<= val1 ; j++) {
          rdline, f;
          rdline, f;
          rdline, f;
        }

        rdline, f;
        rdline, f;
      }
    } else {
      //Disk
      for (i=1 ; i <=  param.simu.n_zones ; i++) {
        read, f, val1 ; param.zones(i).geometry = val1 ;
        read, f, val1; param.zones(i).mass = val1;
        read, f, val1, val2; param.zones(i).Ho = val1; param.zones(i).Ro = val2;
        read, f, val1, val2, val3, val4;  param.zones(i).rin = val1; param.zones(i).rout = val2; param.zones(i).size_neb = val3 ; param.zones(i).edge = val4;
        if (param.zones(i).size_neb > param.map.size_neb) param.map.size_neb = param.zones(i).size_neb ;
        read, f, val1; param.zones(i).beta = val1;
        read, f, val1; param.zones(i).dens = val1;
      }
      rdline, f;
      rdline, f;


      //Grain
      for (i=1 ; i <=  param.simu.n_zones ; i++) {
        read, f, val1, com1; param.zones(i).n_especes=val1;
        for (j=1 ; j<= param.zones(i).n_especes ; j++) {
          read, f, string1, val1, val2; param.dust_pop(i,j).file = string1;  param.dust_pop(i,j).dens = val1; param.dust_pop(i,j).mass_fraction = val2;
          read, f, val1; param.dust_pop(i,j).type_grain = val1;
          read, f, val1, val2, val3, val4 ; param.dust_pop(i,j).amin=val1;  param.dust_pop(i,j).amax=val2;  param.dust_pop(i,j).aexp=val3; param.dust_pop(i,j).n_grains=val4;
        }
      }
      if (param.dust_pop(1,1).type_grain == 2) {
        param.gridT.n_T_grain = sum(param.dust_pop.n_grains);
      } else {
        param.gridT.n_T_grain = 1;
      }
      rdline, f;
      rdline, f;

    }

    // mol
    read, f, val1, val2, val3 ; param.mol.vmax=val1 ; param.mol.vturb=val2 ; param.mol.n_speed = val3;
    read, f, string1, string2, string3, val1; param.mol.lpop=string1 ; param.mol.lpop_precise=string2 ; param.mol.lLTE=string3 ; param.mol.profile_width = val1;
    read, f, string1, val1; param.mol.molecule_file=string1 ; param.mol.level_max = val1;
    read, f, string1, val1, string2 ; param.mol.cst_abundance=string1 ; param.mol.abundance=val1 ; param.mol.abundance_file = string2;
    read, f, string1, val1; param.mol.ray_tracing=string1 ; param.mol.n_lines_ray_tracing=val1;
    if (param.mol.n_lines_ray_tracing == n_lines_ray_tracing)  {
      if (param.mol.n_lines_ray_tracing == 1) {
        read, f, val1 ; param.mol.transition_numbers = val1;
      } else if (param.mol.n_lines_ray_tracing == 2) {
        read, f, val1, val2 ; param.mol.transition_numbers = [val1,val2];
      } else if (param.mol.n_lines_ray_tracing == 3) {
        read, f, val1, val2, val3 ; param.mol.transition_numbers = [val1,val2,val3];
      } else if (param.mol.n_lines_ray_tracing == 4) {
        read, f, val1, val2, val3, val4 ; param.mol.transition_numbers = [val1,val2,val3,val4];
      } else if (param.mol.n_lines_ray_tracing == 5) {
        read, f, val1, val2, val3, val4, val5 ; param.mol.transition_numbers = [val1,val2,val3,val4,val5];
      }
    }
    rdline, f;
    rdline, f;

    // Star
    read, f, val1, com1; param.simu.n_stars = val1;
    if (param.simu.n_stars == n_stars) {
      for (i=1 ; i <= param.simu.n_stars ; i++) {
        read, f, val1, val2, val3, val4, val5, val6, string1; param.stars(i).T = val1; param.stars(i).R = val2 ;  param.stars(i).Mass = val3;
        param.stars(i).x = val4 ; param.stars(i).y = val5 ; param.stars(i).z = val6 ;
        param.stars(i).is_bb = string1;
        read, f, param.stars(i).spectre ;
      }
    }
  } else if (abs(param.simu.version -2.07) < 1.0e-6) {
    rdline, f;

    //photon;
    read, f, val1, com1; param.phot.n1 = val1;
    read, f, val1, com1; param.phot.n2 = val1;
    read, f, val1, com1; param.phot.nlambda = val1;
    read, f, val1, com1; param.phot.nimage = val1;
    read, f, string1, val1; param.simu.lcheckpoint = string1 ; param.simu.checkpoint_time = val1;

    rdline, f;
    rdline, f;

    //Wavelength
    read, f, val1, val2, val3, com1;
    param.wave.nlamb1 = val1; param.wave.lambda_min = val2; param.wave.lambda_max = val3; lambda_min=val2 ; lambda_max=val3;
    read, f, string1, string2, string3 ; param.simu.ltemp = string1; param.simu.lsed = string2; param.simu.lsed_complete = string3;

    read, f, com1, com2; param.wave.file = com1;
    read, f, string1 ; param.simu.lem_disk_image = string1;
    read, f, string1, string2 ; param.simu.lsepar = string1; param.simu.lsepar_pola = string2;
    param.map.n_type=1;
    if (param.simu.lsepar=="T") param.map.n_type += 4;
    if (param.simu.lsepar_pola=="T") param.map.n_type += 3;
    read, f, val1, val2 ; param.simu.tau_seuil = val1 ; param.simu.wl_seuil = val2;

    rdline, f;
    rdline, f;

    //Grid
    read, f, val1; param.grid.geometry = val1;
    read, f, val1, val2, val3, val4, com1;
    param.grid.nrad=val1; param.grid.nz=val2; param.grid.n_az=val3 ; param.grid.nrad_in=val4;

    rdline, f;
    rdline, f;

    //map;
    read, f, val1, val2, val3, val4, val5, com1;
    param.map.nthet = val1;   param.map.nphi = val2
    param.map.nx = val3; param.map.ny=val4; param.map.zoom = val5;
    read, f, val1, val2, val3, string1; param.map.inc = val1; param.map.delta = val2; param.map.angle_interest =val3;param.map.lonly_bin_interest = string1;
    if (param.map.lonly_bin_interest=="T") {
      capt_inf = max(1,param.map.inc-param.map.delta);
      capt_sup = min(param.map.nthet,param.map.inc+param.map.delta);
      param.map.nthet = capt_sup-capt_inf+1;
    }
    read, f, val1, com1; param.map.dist = val1;

    rdline, f;
    rdline, f;

    // Scattering method
    read, f, val1 ; param.simu.scatt_method = val1;
    read, f, val1 ; param.simu.aniso_method = val1;

    rdline, f;
    rdline, f;

    // Symetries
    read, f, string1; param.simu.lsym_ima = string1 ;
    read, f, string1; param.simu.lsym_centrale = string1 ;
    read, f, string1; param.simu.lsym_axiale = string1 ;

    rdline, f;
    rdline, f;

    // Dust global properties
    read, f, val1 ; param.dust.gas_dust = val1;
    read, f, string1, val1, com2; param.dust.lstrat = string1 ;  param.dust.strat = val1;
    if (param.dust.lstrat=="T") {
      param.dust.strat_type=1 ;
    } else {
      param.dust.strat_type=0 ;
    }
    read, f, string1; param.simu.ldust_sublimation = string1 ;
    read, f, string1 , val1 ; param.simu.lchauff_int = string1 ; param.simu.viscosity = val1;
    read, f, val1, val2, val3; param.gridT.T_min=val1; param.gridT.T_max=val2; param.gridT.n_T=val3;

    read, f;
    read, f;

    //Nbre de zones
    read,f, val1 ; param.simu.n_zones = val1;

    read, f;
    read, f;

    if (n_zones!=param.simu.n_zones) { // on veut juste lire n_especes
      for (i=1 ; i <=  param.simu.n_zones ; i++) {
        rdline, f;
        rdline, f;
        rdline, f;
        rdline, f;
        rdline, f;
        rdline, f;

        rdline, f;
        rdline, f;
        read, f, val1, com1;
        for (j=1 ; j<= val1 ; j++) {
          rdline, f;
          rdline, f;
          rdline, f;
        }

        rdline, f;
        rdline, f;
      }
    } else {
      //Disk
      for (i=1 ; i <=  param.simu.n_zones ; i++) {
        read, f, val1 ; param.zones(i).geometry = val1 ;
        read, f, val1; param.zones(i).mass = val1;
        read, f, val1, val2; param.zones(i).Ho = val1; param.zones(i).Ro = val2;
        read, f, val1, val2, val3, val4;  param.zones(i).rin = val1; param.zones(i).rout = val2; param.zones(i).size_neb = val3 ; param.zones(i).edge = val4;
        if (param.zones(i).size_neb > param.map.size_neb) param.map.size_neb = param.zones(i).size_neb ;
        read, f, val1; param.zones(i).beta = val1;
        read, f, val1; param.zones(i).dens = val1;
      }
      rdline, f;
      rdline, f;


      //Grain
      for (i=1 ; i <=  param.simu.n_zones ; i++) {
        read, f, val1, com1; param.zones(i).n_especes=val1;
        for (j=1 ; j<= param.zones(i).n_especes ; j++) {
          read, f, string1, val1, val2; param.dust_pop(i,j).file = string1;  param.dust_pop(i,j).dens = val1; param.dust_pop(i,j).mass_fraction = val2;
          read, f, val1; param.dust_pop(i,j).type_grain = val1;
          read, f, val1, val2, val3, val4 ; param.dust_pop(i,j).amin=val1;  param.dust_pop(i,j).amax=val2;  param.dust_pop(i,j).aexp=val3; param.dust_pop(i,j).n_grains=val4;
        }
      }
      if (param.dust_pop(1,1).type_grain == 2) {
        param.gridT.n_T_grain = sum(param.dust_pop.n_grains);
      } else {
        param.gridT.n_T_grain = 1;
      }
      rdline, f;
      rdline, f;

    }

    // mol
    read, f, val1, val2, val3 ; param.mol.vmax=val1 ; param.mol.vturb=val2 ; param.mol.n_speed = val3;
    read, f, string1, string2, string3, val1; param.mol.lpop=string1 ; param.mol.lpop_precise=string2 ; param.mol.lLTE=string3 ; param.mol.profile_width = val1;
    read, f, string1, val1; param.mol.molecule_file=string1 ; param.mol.level_max = val1;
    read, f, string1, val1, string2 ; param.mol.cst_abundance=string1 ; param.mol.abundance=val1 ; param.mol.abundance_file = string2;
    read, f, string1, val1; param.mol.ray_tracing=string1 ; param.mol.n_lines_ray_tracing=val1;
    if (param.mol.n_lines_ray_tracing == n_lines_ray_tracing)  {
      if (param.mol.n_lines_ray_tracing == 1) {
        read, f, val1 ; param.mol.transition_numbers = val1;
      } else if (param.mol.n_lines_ray_tracing == 2) {
        read, f, val1, val2 ; param.mol.transition_numbers = [val1,val2];
      } else if (param.mol.n_lines_ray_tracing == 3) {
        read, f, val1, val2, val3 ; param.mol.transition_numbers = [val1,val2,val3];
      } else if (param.mol.n_lines_ray_tracing == 4) {
        read, f, val1, val2, val3, val4 ; param.mol.transition_numbers = [val1,val2,val3,val4];
      } else if (param.mol.n_lines_ray_tracing == 5) {
        read, f, val1, val2, val3, val4, val5 ; param.mol.transition_numbers = [val1,val2,val3,val4,val5];
      }
    }
    rdline, f;
    rdline, f;

    // Star
    read, f, val1, com1; param.simu.n_stars = val1;
    if (param.simu.n_stars == n_stars) {
      for (i=1 ; i <= param.simu.n_stars ; i++) {
        read, f, val1, val2, val3, val4, val5, val6, string1; param.stars(i).T = val1; param.stars(i).R = val2 ;  param.stars(i).Mass = val3;
        param.stars(i).x = val4 ; param.stars(i).y = val5 ; param.stars(i).z = val6 ;
        param.stars(i).is_bb = string1;
        read, f, param.stars(i).spectre ;
      }
    }
  } else if (abs(param.simu.version -2.05) < 1.0e-6) {
    rdline, f;

    //photon;
    read, f, val1, com1; param.phot.n1 = val1;
    read, f, val1, com1; param.phot.n2 = val1;
    read, f, val1, com1; param.phot.nlambda = val1;
    read, f, val1, com1; param.phot.nimage = val1;
    read, f, val1, com1; // mol
    rdline, f; // skip checkpointing

    rdline, f;
    rdline, f;

    //Wavelength
    read, f, val1, val2, val3, com1;
    param.wave.nlamb1 = val1; param.wave.lambda_min = val2; param.wave.lambda_max = val3; lambda_min=val2 ; lambda_max=val3;
    read, f, string1, string2, string3 ; param.simu.ltemp = string1; param.simu.lsed = string2; param.simu.lsed_complete = string3;
    read, f, com1, com2; param.wave.file = com1;
    rdline, f; //SKIP consider disk emission in images
    read, f, string1, string2; param.simu.lsepar = string1 ; param.simu.lsepar_pola = string2;
    param.map.n_type=1;
    if (param.simu.lsepar=="T") param.map.n_type += 4;
    if (param.simu.lsepar_pola=="T") param.map.n_type += 3;
    rdline, f; //SKIP tau_max @lambda

    rdline, f;
    rdline, f;

    //Grid
    read, f, val1; // type grille
    read, f, val1, val2, val3, val4, com1;
    param.grid.nrad=val1; param.grid.nz=val2; param.grid.nrad_in=val4;

    rdline, f;
    rdline, f;

    //map;
    read, f, val1, val2, val3, val4, val5, com1;
    param.map.nthet = val1;   param.map.nphi = val2
    param.map.nx = val3; param.map.ny=val4; param.map.zoom = val5;
    read, f, val1, val2, val3, string1; param.map.inc = val1; param.map.delta = val2; param.map.lonly_bin_interest = string1;
    if (param.map.lonly_bin_interest=="T") param.map.nthet = 1;
    read, f, val1, com1; param.map.dist = val1;

    rdline, f;
    rdline, f;


    // Scattering method
    rdline, f;
    rdline, f;

    rdline, f;
    rdline, f;

    // Symetries
    rdline, f;
    rdline, f;
    rdline, f;

    rdline, f;
    rdline, f;

    // Dust global properties
    read, f, com1, val1, com2;  param.dust.strat = val1;
    read, f ;// skip sublimation radius
    read, f; // skip viscous heating
    read, f, val1, val2, val3; param.gridT.T_min=val1; param.gridT.T_max=val2; param.gridT.n_T=val3;

    read, f;
    read, f;

    //Nbre de zones
    read,f, val1 ;
    read, f;
    read, f;

    //Disk
    read, f, val1; param.zones.geometry=val1;
    read, f, val1; param.zones.mass = val1;
    read, f, val1, val2; param.zones.Ho = val1; param.zones.Ro = val2;
    read, f, val1, val2, val3, val4;  param.zones.rin = val1; param.zones.rout = val2; param.zones.size_neb = val3 ; param.zones.edge = val4;
    param.map.size_neb = val3;
    read, f, val1; param.zones.beta = val1;
    read, f, val1; param.zones.dens = val1;

    rdline, f;
    rdline, f;

    //Grain
    read, f, val1, com1; // skip n_especes 1 pour le moment
    param.zones(1).n_especes=val1;
    for (i=1 ; i<= param.zones(1).n_especes ; i++) {
      read, f, com1, val1, val2; param.dust_pop(1,i).file = com1;  param.dust_pop(1,i).dens = val1;
      read, f, val1;
      param.dust_pop(1,i).type_grain = val1
      read, f, val1, val2, val3, val4 ;
      param.dust_pop(1,i).amin=val1;  param.dust_pop(1,i).amax=val2;  param.dust_pop(1,i).aexp=val3; param.dust_pop(1,i).n_grains=val4;
    }

    rdline, f;
    rdline, f;

    // mol
    read, f, val1, val2, val3 ; param.mol.vmax=val1 ; param.mol.vturb=val2 ; param.mol.n_speed = val3;
    read, f, string1, string2, string3, string4, val1;
    param.mol.lpop="F" ; param.mol.lpop_precise="F" ; param.mol.lLTE=string1; param.mol.profile_width=val1;
    read, f, string1; param.mol.molecule_file = string1;
    read, f, val1, val2 ;  param.mol.abundance = val2; param.dust.gas_dust = val1;
    param.mol.cst_abundance="T" ; param.mol.abundance_file = "abundance.fits.gz";
    param.mol.ray_tracing = "T" ; param.mol.n_lines_ray_tracing=1 ; param.mol.transition_numbers(1) = 1;



    rdline, f;
    rdline, f;

    // Star
    read, f, val1, com1; param.simu.n_stars = val1;
    read, f, val1, val2, val3; param.stars(1).T = val1;
    param.stars(1).R = val2 ;  param.stars(1).Mass = val3;

  }  else if (abs(param.simu.version -2.04) < 1.0e-6) {
    rdline, f;

    //photon;
    read, f, val1, com1; param.phot.n1 = val1;
    read, f, val1, com1; param.phot.n2 = val1;
    read, f, val1, com1; param.phot.nlambda = val1;
    read, f, val1, com1; param.phot.nimage = val1;
    read, f, val1, com1; // mol
    rdline, f; // skip checkpointing

    rdline, f;
    rdline, f;

    //Wavelength
    read, f, val1, val2, val3, com1;
    param.wave.nlamb1 = val1; param.wave.lambda_min = val2; param.wave.lambda_max = val3; lambda_min=val2 ; lambda_max=val3;
    read, f, string1, string2, string3 ; param.simu.ltemp = string1; param.simu.lsed = string2; param.simu.lsed_complete = string3;
    read, f, com1, com2; param.wave.file = com1;
    rdline, f; //SKIP consider disk emission in images
    rdline, f; //SKIP separation differentes contributions
    rdline, f; //SKIP tau_max @lambda

    param.map.n_type=4;

    rdline, f;
    rdline, f;

    //Grid
    read, f, val1, val2, val3, val4, com1;
    param.grid.nrad=val1; param.grid.nz=val2; param.grid.nrad_in=val4;

    rdline, f;
    rdline, f;

    //map;
    read, f, val1, val2, val3, val4, val5, com1;
    param.map.nthet = val1;   param.map.nphi = val2
    param.map.nx = val3; param.map.ny=val4; param.map.zoom = val5;
    read, f, val1, val2, val3, com1; param.map.inc = val1; param.map.delta = val2;
    read, f, val1, com1; param.map.dist = val1;

    rdline, f;
    rdline, f;


    // Scattering method
    rdline, f;
    rdline, f;

    rdline, f;
    rdline, f;

    // Symetries
    rdline, f;
    rdline, f;
    rdline, f;

    rdline, f;
    rdline, f;

    // Dust global properties
    read, f, com1, val1, com2;  param.dust.strat = val1;
    read, f ;// skip sublimation radius
    read, f; // skip viscous heating
    read, f, val1, val2, val3; param.gridT.T_min=val1; param.gridT.T_max=val2; param.gridT.n_T=val3;

    read, f;
    read, f;

    //Disk
    read, f, val1; param.zones.mass = val1;
    read, f, val1, val2; param.zones.Ho = val1; param.zones.Ro = val2;
    read, f, val1, val2, val3, val4;  param.zones.rin = val1; param.zones.rout = val2; param.map.size_neb = val3 ; param.zones.edge = val4;
    read, f, val1; param.zones.beta = val1;
    read, f, val1; param.zones.dens = val1;

    rdline, f;
    rdline, f;

    //Grain
    read, f, val1, com1; // skip n_especes 1 pour le moment
    param.zones(1).n_especes=val1;
    for (i=1 ; i<= val1 ; i++) {
      read, f, com1, val1, val2; param.dust_pop(1,i).file = com1;  param.dust_pop(1,i).dens = val1;
      read, f, val1;
      param.dust_pop(1,i).type_grain = val1
      read, f, val1, val2, val3, val4 ;
      param.dust_pop(1,i).amin=val1;  param.dust_pop(1,i).amax=val2;  param.dust_pop(1,i).aexp=val3; param.dust_pop(1,i).n_grains=val4;
    }

    rdline, f;
    rdline, f;

    // mol
    rdline, f;
    rdline, f;
    rdline, f;

    rdline, f;
    rdline, f;

    // Star
    read, f, val1, com1; param.simu.n_stars = val1;
    read, f, val1, val2, val3; param.stars(1).T = val1;
    param.stars(1).R = val2 ;  param.stars(1).Mass = val3;

  }  else if (abs(param.simu.version -2.03) < 1.0e-6) {
    rdline, f;

    //photon;
    read, f, val1, com1; param.phot.n1 = val1;
    read, f, val1, com1; param.phot.n2 = val1;
    read, f, val1, com1; param.phot.nlambda = val1;
    read, f, val1, com1; param.phot.nimage = val1;
    rdline, f; // skip checkpointing

    rdline, f;
    rdline, f;

    //Wavelength
    read, f, val1, val2, val3, com1;
    param.wave.nlamb1 = val1; param.wave.lambda_min = val2; param.wave.lambda_max = val3; lambda_min=val2 ; lambda_max=val3;
    read, f, string1, string2, string3; param.simu.ltemp=string1 ; param.simu.lsed=string2 ; param.simu.lsed_complete=string3;//skip  ltemp, lsed, lsed_complete
    read, f, com1, com2; param.wave.file = com1;
    rdline, f; //SKIP consider disk emission in images
    rdline, f; //SKIP separation differentes contributions
    rdline, f; //SKIP tau_max @lambda

    rdline, f;
    rdline, f;

    //Grid
    read, f, val1, val2, val3, val4, com1;
    param.grid.nrad=val1; param.grid.nz=val2; param.grid.nrad_in=val4;

    rdline, f;
    rdline, f;

    //map;
    read, f, val1, val2, val3, val4, val5, com1;
    param.map.nthet = val1;   param.map.nphi = val2
    param.map.nx = val3; param.map.ny=val4; param.map.zoom = val5;
    read, f, val1, val2, val3, com1; param.map.inc = val1; param.map.delta = val2;
    read, f, val1, com1; param.map.dist = val1;

    rdline, f;
    rdline, f;


    // Scattering method
    rdline, f;
    rdline, f;

    rdline, f;
    rdline, f;

    // Symetries
    rdline, f;
    rdline, f;
    rdline, f;

    rdline, f;
    rdline, f;

    // Dust global properties
    read, f, com1, val1, com2;  param.dust.strat = val1;
    read, f ;// skip sublimation radius
    read, f; // skip viscous heating
    read, f, val1, val2, val3; param.gridT.T_min=val1; param.gridT.T_max=val2; param.gridT.n_T=val3;

    read, f;
    read, f;

    //Disk
    read, f, val1; param.zones.mass = val1;
    read, f, val1, val2; param.zones.Ho = val1; param.zones.Ro = val2;
    read, f, val1, val2, val3, val4;  param.zones.rin = val1; param.zones.rout = val2; param.map.size_neb = val3 ; param.zones.edge = val4;
    read, f, val1; param.zones.beta = val1;
    read, f, val1; param.zones.dens = val1;

    rdline, f;
    rdline, f;

    //Grain
    read, f, val1, com1; // skip n_especes 1 pour le moment
    param.zones(1).n_especes=val1;
    for (i=1 ; i<= val1 ; i++) {
      read, f, com1, val1, val2; param.dust_pop(1,i).file = com1;  param.dust_pop(1,i).dens = val1;
      read, f, val1;
      param.dust_pop(1,i).type_grain = val1
      read, f, val1, val2, val3, val4 ;
      param.dust_pop(1,i).amin=val1;  param.dust_pop(1,i).amax=val2;  param.dust_pop(1,i).aexp=val3; param.dust_pop(1,i).n_grains=val4;
    }

    rdline, f;
    rdline, f;

    read, f, val1, com1; param.simu.n_stars = val1;
    read, f, val1, val2, val3; param.stars(1).T = val1;
    param.stars(1).R = val2 ;  param.stars(1).Mass = val3;

  } else if (abs(param.simu.version -2.02) < 1.0e-6) {
    rdline, f;

    //photon;
    read, f, val1, com1; param.phot.n1 = val1;
    read, f, val1, com1; param.phot.n2 = val1;
    read, f, val1, com1; param.phot.nlambda = val1;
    read, f, val1, com1; param.phot.nimage = val1;
    rdline, f; // skip checkpointing

    rdline, f;
    rdline, f;

    //Wavelength
    read, f, val1, com1; param.wave.nlamb1 = val1;
    read, f, string1, string2, string3; param.simu.ltemp=string1; param.simu.lsed=string2; param.simu.lsed_complete=string3;
    read, f, com1, com2; param.wave.file = com1;
    rdline, f; //SKIP consider disk emission in images
    read, f, string1; param.simu.lsepar = string1;
    param.map.n_type = 4;
    if (param.simu.lsepar=="T") param.map.n_type += 4;
    rdline, f; //SKIP tau_max @lambda

    rdline, f;
    rdline, f;

    //Grid
    read, f, val1, val2, val3, val4, com1;
    param.grid.nrad=val1; param.grid.nz=val2; param.grid.nrad_in=val4;

    rdline, f;
    rdline, f;

    //map;
    read, f, val1, val2, val3, val4, val5, com1;
    param.map.nthet = val1;   param.map.nphi = val2
    param.map.nx = val3; param.map.ny=val4; param.map.zoom = val5;
    read, f, val1, val2, val3, com1; param.map.inc = val1; param.map.delta = val2;
    read, f, val1, com1; param.map.dist = val1;

    rdline, f;
    rdline, f;


    // Scattering method
    rdline, f;
    rdline, f;

    rdline, f;
    rdline, f;

    // Symetries
    rdline, f;
    rdline, f;
    rdline, f;

    rdline, f;
    rdline, f;

    // Dust global properties
    read, f, com1, val1, com2;  param.dust.strat = val1;
    read, f ;// skip sublimation radius
    read, f; // skip LTE
    read, f; // skip viscous heating
    read, f, val1, val2, val3; param.gridT.T_min=val1; param.gridT.T_max=val2; param.gridT.n_T=val3;

    read, f;
    read, f;

    //Disk
    read, f, val1; param.zones.mass = val1;
    read, f, val1, val2; param.zones.Ho = val1; param.zones.Ro = val2;
    read, f, val1, val2, val3, val4;  param.zones.rin = val1; param.zones.rout = val2; param.map.size_neb = val3 ; param.zones.edge = val4;
    read, f, val1; param.zones.beta = val1;
    read, f, val1; param.zones.dens = val1;

    rdline, f;
    rdline, f;

    //Grain
    read, f; // skip n_especes 1 pour le moment
    read, f, com1, val1, val2; param.dust_pop(1,1).file = com1;  param.dust_pop(1,1).dens = val1;
    read, f, val1, val2, val3, val4 ;
    param.dust_pop(1,1).amin=val1;  param.dust_pop(1,1).amax=val2;  param.dust_pop(1,1).aexp=val3; param.dust_pop(1,1).n_grains=val4;

    rdline, f;
    rdline, f;

    read, f, val1, com1; param.simu.n_stars = val1;
    read, f, val1, val2, val3; param.stars(1).T = val1;
    param.stars(1).R = val2 ;  param.stars(1).Mass = val3;

  }  else if (abs(param.simu.version -2.01) < 1.0e-6) {
    rdline, f;

    //photon;
    read, f, val1, com1; param.phot.n1 = val1;
    read, f, val1, com1; param.phot.n2 = val1;
    read, f, val1, com1; param.phot.nlambda = val1;
    read, f, val1, com1; param.phot.nimage = val1;
    rdline, f; // skip checkpointing

    rdline, f;
    rdline, f;

    //Wavelength
    read, f, val1, com1; param.wave.nlamb1 = val1;
    read, f, string1, string2, string3; param.simu.ltemp=string1; param.simu.lsed=string2; param.simu.lsed_complete=string3;
    read, f, com1, com2; param.wave.file = com1;
    rdline, f; //SKIP consider disk emission in images
    rdline, f; //SKIP separation differentes contributions
    rdline, f; //SKIP tau_max @lambda

    rdline, f;
    rdline, f;

    //Grid
    read, f, val1, val2, val3, com1;
    param.grid.nrad=val1; param.grid.nz=val2; param.grid.nrad_in=val3;

    rdline, f;
    rdline, f;

    //map;
    read, f, val1, val2, val3, val4, val5, com1;
    param.map.nthet = val1;   param.map.nphi = val2;
    param.map.nx = val3; param.map.ny=val4; param.map.zoom = val5;
    read, f, val1, val2, val3, com1; param.map.inc = val1; param.map.delta = val2;
    read, f, val1, com1; param.map.dist = val1;
    rdline, f;
    rdline, f;


    // Scattering method
    rdline, f;
    rdline, f;

    rdline, f;
    rdline, f;

    // Symetries
    rdline, f;
    rdline, f;
    rdline, f;

    rdline, f;
    rdline, f;

    // Dust global properties
    read, f, com1, val1, com2;  param.dust.strat = val1;
    read, f ;// skip sublimation radius
    read, f; // skip LTE
    read, f; // skip viscous heating
    read, f; // skip Temperature

    read, f;
    read, f;

    //Disk
    read, f, val1; param.zones.mass = val1;
    read, f, val1, val2; param.zones.Ho = val1; param.zones.Ro = val2;
    read, f, val1, val2, val3, val4;  param.zones.rin = val1; param.zones.rout = val2; param.zones.edge = val4;
    read, f, val1; param.zones.beta = val1;
    read, f, val1; param.zones.dens = val1;

    rdline, f;
    rdline, f;

    //Grain
    read, f; // skip n_especes 1 pour le moment
    read, f, com1, val1, val2; param.dust_pop(1,1).file = com1;  param.dust_pop(1,1).dens = val1;
    read, f, val1, val2, val3, val4 ;
    param.dust_pop(1,1).amin=val1;  param.dust_pop(1,1).amax=val2;  param.dust_pop(1,1).aexp=val3;

    rdline, f;
    rdline, f;

    read, f, val1, com1; param.simu.n_stars = val1;
    read, f, val1, val2, val3; param.stars(1).T = val1;
    param.stars(1).R = val2 ;  param.stars(1).Mass = val3;

  }  else if (abs(param.simu.version -2.0) < 1.0e-6) {
    rdline, f;

    //photon;
    read, f, val1, com1; param.phot.n1 = val1;
    read, f, val1, com1; param.phot.n2 = val1;
    read, f, val1, com1; param.phot.nlambda = val1;
    rdline, f; // skip checkpointing

    rdline, f;
    rdline, f;

    //Wavelength
    read, f, val1, com1; param.wave.nlamb1 = val1;
    rdline, f; //skip  ltemp, lsed, lsed_complete
    read, f, com1, com2; param.wave.file = com1;
    rdline, f; //SKIP consider disk emission in images
    rdline, f; //SKIP separation differentes contributions
    rdline, f; //SKIP tau_max @lambda

    rdline, f;
    rdline, f;

    //Grid
    read, f, val1, val2, val3, com1;
    param.grid.nrad=val1; param.grid.nz=val2; param.grid.nrad_in=val3;

    rdline, f;
    rdline, f;

    //map;
    read, f, val1, val2, val3, val4, val5, com1;
    param.map.nthet = val1;   param.map.nphi = val2;
    param.map.nx = val3; param.map.ny=val4; param.map.zoom = val5;
    read, f, val1, val2, val3, com1; param.map.inc = val1; param.map.delta = val2;
    read, f, val1, com1; param.map.dist = val1;
    rdline, f;
    rdline, f;


    // Scattering method
    rdline, f;
    rdline, f;

    rdline, f;
    rdline, f;

    // Symetries
    rdline, f;
    rdline, f;
    rdline, f;

    rdline, f;
    rdline, f;

    // Dust global properties
    read, f, com1, val1, com2;  param.dust.strat = val1;
    read, f ;// skip sublimation radius
    read, f; // skip LTE
    read, f; // skip viscous heating
    read, f; // skip T_min
    read, f; // skip T_max

    read, f;
    read, f;

    //Disk
    read, f, val1; param.zones.mass = val1;
    read, f, val1, val2; param.zones.Ho = val1; param.zones.Ro = val2;
    read, f, val1, val2, val3, val4;  param.zones.rin = val1; param.zones.rout = val2; param.zones.edge = val4;
    read, f, val1; param.zones.beta = val1;
    read, f, val1; param.zones.dens = val1;

    rdline, f;
    rdline, f;

    //Grain
    read, f; // skip n_especes 1 pour le moment
    read, f, com1, val1, val2; param.dust_pop(1,1).file = com1;  param.dust_pop(1,1).dens = val1;
    read, f, val1, val2, val3, val4 ;
    param.dust_pop(1,1).amin=val1;  param.dust_pop.amax(1,1)=val2;  param.dust_pop(1,1).aexp=val3;

    rdline, f;
    rdline, f;

    read, f, val1, com1; param.simu.n_stars = val1;
    read, f, val1, val2, val3; param.stars(1).T = val1;
    param.stars(1).R = val2 ;  param.stars(1).Mass = val3;

  } else if (abs(param.simu.version - 1.9) < 1.e-6) {
    rdline, f;

    //photon;
    read, f, val1, com1; param.phot.n1 = val1;
    read, f, val1, com1; param.phot.n2 = val1;
    read, f, val1, com1; param.phot.nlambda = val1;
    rdline, f; // skip checkpointing

    rdline, f;
    rdline, f;

    //Wavelength
    read, f, val1, com1; param.wave.nlamb1 = val1;
    rdline, f; //skip  ltemp, lsed, lsed_complete
    read, f, com1, com2; param.wave.file = com1;
    rdline, f; //SKIP consider disk emission in images
    rdline, f; //SKIP separation differentes contributions
    rdline, f; //SKIP tau_max @lambda

    rdline, f;
    rdline, f;

    //Grid
    read, f, val1, val2, val3, com1;
    param.grid.nrad=val1; param.grid.nz=val2; param.grid.nrad_in=val3;

    rdline, f;
    rdline, f;

    //map;
    read, f, val1, val2, val3, val4, val5, com1;
    param.map.nthet = val1;   param.map.nphi = val2;
    param.map.nx = val3; param.map.ny=val4; param.map.zoom = val5;
    read, f, val1, val2, val3, com1; param.map.inc = val1; param.map.delta = val2;
    read, f, val1, com1; param.map.dist = val1;
    rdline, f;
    rdline, f;


    // Scattering method
    rdline, f;
    rdline, f;

    rdline, f;
    rdline, f;

    // Symetries
    rdline, f;
    rdline, f;
    rdline, f;

    rdline, f;
    rdline, f;

    // Dust global properties
    read, f, com1, val1, com2;  param.dust.strat = val1;
    read, f ;// skip sublimation radius
    read, f; // skip LTE
    read, f; // skip viscous heating
    read, f; // skip T_min
    read, f; // skip T_max

    read, f;
    read, f;

    //Grain
    read, f; // skip n_especes 1 pour le moment
    read, f, com1, val1, val2; param.dust_pop(1,1).file = com1;  param.dust_pop(1,1).dens = val1;
    read, f, val1, val2, val3, val4 ;
    param.dust_pop(1,1).amin=val1;  param.dust_pop(1,1).amax=val2;  param.dust_pop(1,1).aexp=val3;

    rdline, f;
    rdline, f;

    //Disk
    read, f, val1; param.zones.mass = val1;
    read, f, val1, val2; param.zones.Ho = val1; param.zones.Ro = val2;
    read, f, val1, val2; param.zones.rout = val1;
    read, f, val1, val2; param.zones.rin = val1; param.zones.edge = val2;
    read, f, val1; param.zones.beta = val1;
    read, f, val1; param.zones.dens = val1;

    rdline, f;
    rdline, f;

    read, f, val1, com1; param.simu.n_stars = val1;
    read, f, val1, val2, val3; param.stars(1).T = val1;
    param.stars(1).R = val2 ;  param.stars(1).Mass = val3;

  }  else {
    "Version inconnue du fichier de parametres";
    return [];
  }

  return param;
}

func _read_nlambda(file, param) {

  extern nlamb ;
  val1 = int() ;

  if (param.simu.lsed_complete=="T") {
    param.wave.nlamb =  param.wave.nlamb1;
    nlamb = param.wave.nlamb ;
    if (nlamb /=  param.wave.nlamb) {
      nlamb = param.wave.nlamb ;
      include, struct_file, 1 ;
    }
  } else {
    if (param.simu.lsed=="T") {
      wave_dir = string();
      f = popen("dirname "+file,0);
      read, f, wave_dir;
      close, f;

      if (strpart(wave_dir,-1:0) == "th") {
        f = popen("wc -l "+wave_dir+"/"+param.wave.file,0);
        read, f, val1;
        param.wave.nlamb=val1;
        close, f;
      }
    }
  }
  return param ;
}


func read_params(file,resol_x,resol_y,zoom) {

  param = _read_params(file) ;
  param = _read_nlambda(file,param) ; // in case _read_params returns before finishing reading the file

  if (is_void(param)) return [];

  // 1er test sur n_zones
  if (param.simu.n_zones!=n_zones) {
    n_zones = param.simu.n_zones;
    include, struct_file, 1;
    param = _read_params(file);
  }

  // tests sur les autres parametres
  if (param.map.nthet!=nthet || param.map.nphi!=nphi || param.wave.nlamb1!=nlamb1 || param.wave.nlamb!=nlamb || param.simu.n_especes_max!=n_especes || param.simu.n_stars!=n_stars || param.grid.nrad!=nrad || param.grid.nz!=nz || param.grid.nrad_in!=nrad_in || param.grid.n_az!=n_az || param.map.nx!=nx || param.map.ny!=ny || param.map.n_type!=n_type || param.gridT.n_T_grain!=n_T_grain || param.map.RT_ntheta!=RT_ntheta  || param.simu.nmol!=nmol || (param.mol.n_speed!=n_speed)(1) || (n_lines_ray_tracing!=param.mol.n_lines_ray_tracing)(1) || param.simu.n_components_max!=n_components || param.simu.n_zones!=n_zones) {


    nlamb1 =   param.wave.nlamb1;
    nlamb = param.wave.nlamb;
    nphi   = param.map.nphi;
    nthet  =  param.map.nthet;
    RT_ntheta  =  param.map.RT_ntheta;
    nz =  param.grid.nz;
    nrad = param.grid.nrad;
    n_az = param.grid.n_az;
    n_rad_in = param.grid.nrad_in;
    n_zones = param.simu.n_zones ;
    n_especes = param.simu.n_especes_max ;
    n_components = param.simu.n_components_max ;
    n_stars = param.simu.n_stars;

    n_speed = param.mol.n_speed(1);
    n_lines_ray_tracing = param.mol.n_lines_ray_tracing(1);

    nx = param.map.nx;
    ny = param.map.ny;

    n_type = param.map.n_type;
    n_T_grain = param.gridT.n_T_grain;
    n_speed = param.mol.n_speed(1) ;
    nmol = param.simu.nmol ;

    // Parametres passes en option de mcfost
    if (!is_void(resol_x)) {
      param.map.nx = resol_x ;
      param.map.ny = resol_y ;
      nx = resol_x ;
      ny = resol_y ;
    }
    if (!is_void(zoom)) {
      param.map.zoom = zoom ;
    }

    include, struct_file, 1;
    param = _read_params(file);

    // Parametres passes en option de mcfost
    if (!is_void(resol_x)) {
      param.map.nx = resol_x ;
      param.map.ny = resol_y ;
      nx = resol_x ;
      ny = resol_y ;
    }
    if (!is_void(zoom)) {
      param.map.zoom = zoom ;
    }
  }

  return param;
}

/*func _read_T (dir, nrad, nz) {
  file = dir+data_dir+temperature_file;
  f = open(file);
  i = j = array(int, nrad*nz);
  T1 = array(double, nrad*nz);
  read, f, i, j, T1;
  T = array(double, nz , nrad);
  T(*) = T1;
  return transpose(T);
  }*/

func _read_T (dir, nrad, nz, v) {
  file = dir+data_dir+temperature_file;
  if ( !open(file, "r", 1)) {
    if (v) "the file "+file+" does not exist";
    return [];
  }
  T =cfitsRead(file);
  return T;
}


/*func _read_dim_sed(dir) {
  file = dir+data_dir+sed2_file;
  f = open (file);
  ninc = nphi = nlamb = int();
  read, f, ninc, nphi, nlamb;
  return [ninc, nphi, nlamb];
  }*/



/*func _read_sed1(dir, &lamb) {
  file = dir+data_dir+sed1_file;

  if ( !open(file, "r", 1)) {
    error, "the file "+file+" does not exist";
    return [];
  }
  f = open (file);
  ninc = nphi = nlamb = int();
  read, f, ninc, nphi, nlamb;

  sed1 = array(double, nlamb, ninc, nphi, nsed1);
  for (j=1; j<=ninc ; j++) {
    for (k=1; k<=nphi ; k++) {
      rdline, f;

      flux1 = flux2 = lamb = array(double,nlamb);

      read, f, lamb ,  flux1 , flux2 ;
      tmp_flux = [flux1 , flux2];

      lamb     = lamb   * lamb_adjust;
      sed1(,j,k,) = tmp_flux( ,vect_n_sed1)  * flux_adjust;

    }
  }
  close, f;
  return sed1( ,vect_n_inc, vect_n_phi,);
  }*/


func _read_sed1(dir, v) {
  file = dir+data_dir+sed1_file;

  if ( !open(file, "r", 1)) {
    if (v) "the file "+file+" does not exist";
    return [];
  }

  sed1 = cfitsRead(file);
  return sed1;
}


func _calc_lambda1(n_lambda) {

  extern lambda_min, lambda_max ;

  tab_lambda= array(double, n_lambda) ;

  //delta_lambda = (lambda_max/lambda_min)**(1.0/real(n_lambda))
   delta_lambda =  exp((1.0/double(n_lambda)) * log(lambda_max/lambda_min));

   tab_lambda(1)=lambda_min*sqrt(delta_lambda);
   for (i=2 ;  i<=n_lambda ; i++) {
     tab_lambda(i)= tab_lambda(i-1)*delta_lambda;
   }
   return tab_lambda
}

func _calc_rt_incl(P) {

  if (P.map.RT_ntheta==1) {
    return [P.map.RT_imin] ;
  } else {
    cos_min = cos(P.map.RT_imin / 180.0 * pi) ;
    cos_max = cos(P.map.RT_imax / 180.0 * pi) ;
    if (P.map.lRT_i_centered=="T") {
         return acos(  cos_min + (indgen(P.map.RT_ntheta) - 0.5)/P.map.RT_ntheta * (cos_max - cos_min) ) /pi * 180.0 ;
       } else {
         return acos(  cos_min + (indgen(P.map.RT_ntheta) - 1.0)/P.map.RT_ntheta * (cos_max - cos_min) ) /pi * 180.0 ;
       }
  }
}

/*func _read_sed2(dir, &lamb) {
  file = dir+data_dir+sed2_file;

  if ( !open(file, "r", 1)) {
    error, "the file "+file+" does not exist";
    return [];
  }
  f = open (file);
  ninc = nphi = nlamb = int();
  read, f, ninc, nphi, nlamb;

  sed = array(double, nlamb, ninc, nphi, nsed);
  for (j=1; j<=ninc ; j++) {
    for (k=1; k<=nphi ; k++) {
      rdline, f;

      flux1 = flux2 = flux3 = flux4 = flux5 = flux6 = flux7 = flux8 = flux9 =
      flux10 = lamb = array(double,nlamb);

      read, f, lamb ,  flux1 , flux2 , flux3 , flux4 , flux5 , flux6 , flux7 , flux8 , flux9 ,flux10;
      tmp_flux = [flux1 , flux2 , flux3 , flux4 , flux5 , flux6 , flux7 , flux8, flux9 ,flux10];

      lamb     = lamb   * lamb_adjust;
      sed(,j,k,) = tmp_flux( ,vect_n_flux)  * flux_adjust;

    }
  }
  close, f;
  return sed( ,vect_n_inc, vect_n_phi,);
  } */


func _read_sed2(dir,v) {
  file = dir+data_dir+sed2_file;

  if ( !open(file, "r", 1)) {
    if (v)  "the file "+file+" does not exist";
    return [];
  }

  sed2=cfitsRead(file);

  return sed2;
}

func _read_sed_rt(dir,v) {
  file = dir+data_dir+sed_rt_file;

  if ( !open(file, "r", 1)) {
    if (v)  "the file "+file+" does not exist";
    return [];
  }

  sed_rt=cfitsRead(file);
  return sed_rt;
}


func _read_lambda2(file) {

  /*
  n_lambda2 = int();
  val1 = float();
  com1 = string();

  cmd = "wc -l "+file+" > nbr_lignes.txt";
  system,  cmd;
  f = open("nbr_lignes.txt");
  read, f, n_lambda2;
  close, f;
  system, "rm nbr_lignes.txt";

  lambda2 = array(double, n_lambda2);
  f = open(file);
  for (i=1 ; i<=n_lambda2 ; i++) {
  read, f, val1;
  lambda2(i) = val1;
  }
  close, f;
  */

  f = open(file,"r") ;
  bmark = bookmark(f) ;

  // On compte les lignes
  for (n=1; (line=rdline(f)) && line!="";n++) ;
  n_lines = n-1 ;

  // On lit le tableau
  lambda2 = array(double,n_lines) ;
  backup, f, bmark ;
  read, f, lambda2 ;

  close, f ;

  return lambda2;
}

// Ouverture d'une liste de modèles
func _open_mcfost_sed(dir,v) {
  extern nthet, nphi, nlamb, nz, nrad, vect_n_phi, vect_n_inc ;

  // Recherche du nom du fichier de parametres
  dir2 = dir+data_dir;
  f = popen("ls "+dir2+"*.par* | wc -l",0) ;
  n_param = int(); read, f, n_param;
  close, f;
  if (n_param != 1) {
    "Error : I found several parameter files";
    return [];
  }
  f= popen("ls "+dir2+"*.par*",0);
  file = string() ; read, f, file;
  close, f;

  param = read_params(file); if (is_void(param)) return [];

  model = McfostSED(dir = dir);
  model.P = param;

  sed1 = _read_sed1(dir,v);
  sed2 = _read_sed2(dir,v);

  if (IsFile(dir2+"sed_rt.fits.gz") == 1) {
    sed_rt = _read_sed_rt(dir,v) ;
    model.P.map.RT_incl = _calc_rt_incl(model.P) ;
  }
  close, f ;

  T = _read_T (dir, model.P.grid.nrad, model.P.grid.nz,v);

  // calcul lambda step1
  model.lamb1 = _calc_lambda1(model.P.wave.nlamb1);
  if ((nlamb==nlamb1) && (!is_void(sed2)))  model.lamb=model.lamb1; // Pb si nlamb change

  if ((param.simu.lsed_complete=="F") && (param.simu.lsed=="T")) {
    wave_dir = string();
    f = popen("dirname "+file,0);
    read, f, wave_dir;
    close, f;

    model.lamb = _read_lambda2(wave_dir+"/"+param.wave.file);
  }


  if (!is_void(sed1)) {      // Bug ici !!!!
    model.sed1=sed1;
    sed1=[];
  }

  if (!is_void(sed2)) {
    model.sed = sed2;
  }

  if (!is_void(sed_rt)) {
    model.sed_rt = sed_rt;
    sed_rt = [];
  }

  if (!is_void(T)) {
    //model.T  = T(,,-:1:1,) ; // Bug ici  ??
    T=[];
  }
  return model;
}

func _open_mcfost_image(dir,v,resol_x,resol_y,zoom) {
  extern nthet, nphi, nx, ny, n_type, vect_n_phi, vect_n_inc ;

  // Recherche du nom du fichier de parametres
  dir2 = dir+data_dir;
  f = popen("ls "+dir2+"*.par* | wc -l",0);
  n_param = int(); read, f, n_param;
  if (n_param != 1) {
    "Error : multiple parameter files";
    return [];
  }
  close, f ;
  f= popen("ls "+dir2+"*.par*",0);
  file = string() ; read, f, file;
  close, f;
  param = read_params(file,resol_x,resol_y,zoom); if (is_void(param)) return [];
  model = McfostImage(dir = dir);

  model.P = param;

  if ( !open(dir+data_dir+image_file, "r", 1)) {
    if (v) "the file "+dir+data_dir+image_file+" does not exist" ;
    //    return [] ;
  } else {
    image = cfitsRead(dir+data_dir+image_file);
    dim_image = dimsof(image) ;
    model.image = image;
  }

  if ( !open(dir+data_dir+image_rt_file, "r", 1)) {
    if (v) "the file "+dir+data_dir+image_rt_file+" does not exist" ;
    //  return [] ;
  } else {
    image = cfitsRead(dir+data_dir+image_rt_file);
    dim_image = dimsof(image) ;
    model.image_rt = image;
  }

  // Il manque lambda
  return model;
}

func _open_mcfost_spectre(dir,v) {

  // Recherche du nom du fichier de parametres
  dir2 = dir+data_dir;
  f = popen("ls "+dir2+"*.par* | wc -l",0);
  n_param = int(); read, f, n_param;
  close, f ;

  f = popen("ls "+dir2+"*.lines* 2> /dev/null | wc -l",0);
  n_lines = int(); read, f, n_lines;
  close, f ;

  if (n_param + n_lines != 1) {
    "Error : multiple parameter files";
    return [];
  }

  model = McfostSpectre(dir = dir);

  if (n_param) {
    f= popen("ls "+dir2+"*.par*",0);
    file = string() ; read, f, file;
    close, f;

    param = read_params(file,resol_x,resol_y,zoom); if (is_void(param)) return [];
    model = McfostSpectre(dir = dir); // on redefinit la structure

    model.P = param;
  }


  if ( !open(dir+data_dir+spectre_file, "r", 1)) {
    if (v) "the file "+dir+data_dir+spectre_file+" does not exist" ;
    return [] ;
  }

  fits = cfitsio_open(dir+data_dir+spectre_file, "r");

  Fline = cfitsio_read_image(fits);
  model.Fline = Fline;
  cfitsio_goto_hdu,fits,2;
  Fcont = cfitsio_read_image(fits);
  model.Fcont = Fcont;
  cfitsio_goto_hdu,fits,3;
  ifreq = cfitsio_read_image(fits);
  model.ifreq = ifreq;
  cfitsio_goto_hdu,fits,4;
  freq = cfitsio_read_image(fits);
  model.freq = freq;
  cfitsio_goto_hdu,fits,5;
  velocity = cfitsio_read_image(fits);
  model.velocity = velocity;
  cfitsio_close, fits ;

  return model ;
}

func open_mcfost(dir,type,v=,resol_x=,resol_y=,zoom=) {
  /* DOCUMENT
     arguments : dir, type, v=, resol_x=, resol_y=, zoom=
  SEE ALSO:
  */

  extern data_dir;
  data_dir="/data_"+type+"/";

  if (is_void(v)) v=0 ; // verbose

  // Liste de molecules reconnues par yorick
  molecules = ["C+","O","CO","o-H2O","p-H2O","13C16O","CS","C18O"] ;

  if (type == "th") {
    file = sed2_file;
  } else if (numberof(where(type == molecules))){
    mol_name = type ;
    file = spectre_file ;
  } else {
     file = image_file;
  };

  if (!dimsof(dir)(1)) dir = [dir];
  size = numberof(dir);
  test_dir = array(0, size);

  // control if all files are OK
  for (i=1; i<=size; i++) {
    if (open(dir(i)+data_dir+"/"+file ,"r" , 1))
      test_dir(i) = 1;
    else {
      if (v) write, "File : "+dir(i)+data_dir+"/"+file+" does not exist";
      test_dir(i) = 1;
    }
  }
  dir = dir(where(test_dir));
  size= numberof(dir);

  model = [];
  if (type == "th") {
    for (i=1; i<=size; i++) {
      mod = _open_mcfost_sed(dir(i),v) ;
      if (!is_void(mod)) grow, model ,  mod;
    }

  } else if (numberof(where(type == molecules))){
    for (i=1; i<=size; i++) {
      mod = _open_mcfost_spectre(dir(i),v) ;
      if (!is_void(mod)) {
        mod.mol_name = type ;
        grow, model, mod;
      }
    }

  } else {
    for (i=1; i<=size; i++)  {
      mod = _open_mcfost_image(dir(i),v,resol_x,resol_y,zoom);
      if (!is_void(mod)) {
        grow, model , mod;
        lamb=array(float);  buf= sread(type,lamb);  model.lamb=lamb ;
      }
    }
  }

  if (is_void(model)) {
    return [];
  } else {
    //model.grid = grille(model);
  }

  //if (sizeof(model)) model.dir = dir;
  if (size==1) return model(1);
  return model;
}

func grille(model) {
/* DOCUMENT
     Calcule la grille utilisée par MCFOST
     Donne le rayon du milieu des cellules et l'echelle de hauteur correspondante
   SEE ALSO:
 */

  param = model.P
  dimsP=dimsof(param) ;
  if (dimsP(1) > 1) {
    lmax = dimsP(2);
  } else {
    lmax=1;
  };

  for (l=1 ; l<=lmax ; l++) {
    n_rad = param(l).grid.nrad;
    n_rad_in = param(l).grid.nrad_in;
    nz = param(l).grid.nz;


    r_lim_inf=array(double,n_rad+1);

    rout=max(param(l).zones.rout);
    rmin=min(param(l).zones.rin - 5*param(l).zones.edge);
    exp_beta=param(l).zones.beta;
    surf = param(l).zones.dens;

    ln_delta_r = (1.0/double(n_rad-n_rad_in+1))*log(rout/rmin);
    delta_r = exp(ln_delta_r);

    ln_delta_r_in = (1.0/double(n_rad_in))*log(delta_r);
    delta_r_in = exp(ln_delta_r_in);

    puiss=min(1+surf-exp_beta);

    r_lim_inf(1)=min(rmin);
      // Subdivision premiere cellule
      if (puiss(1) == 0.0) {
        for (i=2 ;  i<=n_rad_in+1; i++) {
          r_lim_inf(i) = exp(log(rmin) - (log(rmin)-log(rmin*delta_r))*(2.0^(i-1)-1.0)/(2.0^n_rad_in-1.0));
        }
      } else {
        for (i=2 ;  i<=n_rad_in+1; i++) {
          r_lim_inf(i) = (rmin^puiss - (rmin^puiss-(rmin*delta_r)^puiss)*(2.0^(i)-1.0)/(2.0^(n_rad_in+1)-1.0))^(1.0/puiss);
        }
      }

    // Autres cellules
    for (i=n_rad_in+2 ;  i<= n_rad+1 ; i++) {
      r_lim_inf(i) = r_lim_inf(i-1) * delta_r;
    }

    model(l).grid.R_grid = 0.5 * (r_lim_inf(1:n_rad) + r_lim_inf(2:n_rad+1));

    // echelle de hauteur
    model(l).grid.H_grid = param(l).zones(1).Ho * (model(l).grid.R_grid()/param(l).zones(1).Ro)^exp_beta(1);

  }
  return model.grid;
}

func eq_hydro(model,M_star,MAX=) {
/* DOCUMENT
   Calcule l'échelle de hauteur à partir de la temperature dans le plan median
   SEE ALSO:
 */
  if (is_void(MAX)) MAX = 0 ;

  mu = 2.3e-3 / SI.Na; // masse molaire moyenne : 2.3 selon Walker 2004
  //  model.h_hydro = 30.0*sqrt((model.T(,1)/15.5)*(0.045/M_star)*(model.grid.R_grid/100.)^3);

  if (MAX == 0) {
    model.h_hydro = sqrt( (SI.kB * model.T(,1) * (model.grid.R_grid*SI.AU)^3) / (SI.Grav * M_star * SI.Msun * mu) ) / SI.AU;
  } else {
    model.h_hydro = sqrt( (SI.kB * model.T(,max) * (model.grid.R_grid*SI.AU)^3) / (SI.Grav * M_star * SI.Msun * mu) ) / SI.AU;
  }

  return model;
}


func H2T(model,M_star) {

/* DOCUMENT
     Calcule la temperature requise pour produire l'echelle de hauteur du modele
   SEE ALSO:
 */
  mu = 2.3e-3 / SI.Na; // masse molaire moyenne : 2.3 selon Walker 2004
  T = (model.grid.H_grid * SI.AU)^2  * (SI.Grav * M_star * SI.Msun * mu)  / (SI.kB * (model.grid.R_grid*SI.AU)^3)

  return T;
}


func sed_verif(model,win,noclean=) {
  win_tmp = window();

  if (is_void(noclean)) noclean=0

  if (is_void(win)) win=20;
  window, win;
  if (!noclean) fma;
  logxy, 1, 1 ; checklog;
  plg, model.sed(,sum,1) , model.lamb, color="red";
  plg, model.sed1(,sum) , model.lamb1, type="dash";

  window, win_tmp;
}


func incl(n) {

  for (i=1 ; i<=n ; i++) {
    print(n-i+1,acos(float(i-1)/float(n))*180/pi,acos(float(i)/float(n))*180/pi, acos((float(i)-0.5)/float(n))*180/pi);
  }

}

//===========================3D functions========================
func R_span (param) {
  param = param;
  if (nameof(structof(param))=="McfostSED") param = param.P;

  return spanl(param.zones.rin , param.zones.rout, param.grid.nrad+1);
}
func R_grid (param) {
  param = param;
  if (nameof(structof(param))=="McfostSED") param = param.P;

  return R_span(param)(, -:1:param.grid.nz+1);
}
func H_grid (param) {
    param = param;
  if (nameof(structof(param))=="McfostSED") param = param.P;
  R = R_grid(param);
  return span(0 , param.zones.Ho, param.grid.nz+1)(-:1:param.grid.nrad+1, ) *
    (R/param.zones.Ro)^param.zones.beta;
}

func disk_model_3d (model, quart= , edges=, index=, Hfact=, alpha=, Tmax=, Tmin=) {

  nsamp_R = model.P.grid.nrad + 1;
  nsamp_H = model.P.grid.nz + 1;
  if (is_void(edges)) edges = 0;
  if (is_void(quart)) quart=pi;
  if (is_void(Hfact)) Hfact =1;
  if (is_void(alpha)) {
    nsamp_alpha = 30;
    alphamin = -pi +quart/2.;
    alphamax =  pi -quart/2.;
    alpha = span(alphamin, alphamax , nsamp_alpha)
  }
  else nsamp_alpha = numberof(alpha);
  if (dimsof(alpha)(sum)(1)<=2) {
    alpha = [alpha(1), alpha(1)];
    nsamp_alpha=2;
  }

  index = indgen(1:nsamp_R)(index);

  R = R_grid(model.P) (index,,-:1:nsamp_alpha);
  H = H_grid(model.P) (index,,-:1:nsamp_alpha)*Hfact;
  alpha = alpha(-:1:nsamp_R,-:1:nsamp_H,)(index,,);

  data = model.T(index(1:-1),,1,-:1:nsamp_alpha-1);

  if (!quart) palpha =0;
  xyz =  cylindre3d(data, R, H , alpha, , pz=[0], palpha=palpha);

  clear3;
  plwf2, _(data, data)   , _(xyz,xyz*[1,1,-1])  , mode=0, edges=edges, cmin=Tmin, cmax=Tmax;
}

func disk3d (rin, rout, beta, Ho, quart=, dup=, edges=) {
  if (is_void(edges)) edges = 0;
  if (is_void(quart)) quart=0;
  nsamp_R = 30;
  nsamp_H = 10;
  nsamp_alpha = 45;

  alphamin = -pi +quart/2.;
  alphamax =  pi -quart/2.;

  Ro = 100;
  R = span(rin, rout, nsamp_R)(,-:1:nsamp_H,-:1:nsamp_alpha);
  H = span(-Ho, Ho, nsamp_H)(-:1:nsamp_R, ,-:1:nsamp_alpha) * (R/Ro)^beta;
  //H = span(-Ho, Ho, nsamp_H)(-:1:nsamp_R, ,-:1:nsamp_alpha);
  alpha = span(alphamin, alphamax , nsamp_alpha)(-:1:nsamp_R,-:1:nsamp_H,);

  //H = H + sign(H)*(alpha>(-pi/4)) * (alpha<(pi/4)) * cos(2*alpha)* (R/rin)^2.5 *.1;

  data2d = (rout - R(zcen,zcen,zcen)) *  abs(H(zcen,zcen,1)/H(zcen,0,1));

  data = data2d;
  if (!quart) palpha =0;
  xyz =  cylindre3d(data, R, H , alpha, palpha=palpha);


  clear3;
  plwf2, data   , xyz  , mode=0, edges=edges;
}


func cylindre3d (&data2d , R , H , alpha , inc, quart=, star= , nsamp_alpha=, dup=,
                 pz=,palpha=,pr=) {
  if (is_void(nsamp_alpha)) nsamp_alpha = 30;
  if (is_void(dup)) dup=1;
  if (is_void(pz)) pz = [1,0];
  if (is_void(pr)) pr = [1,0];
  if (is_void(palpha)) palpha = [1,0];
  if (!dimsof(pz)(1) && pz==0) pz =[];
  if (!dimsof(pr)(1) && pr==0) pr =[];
  if (!dimsof(palpha)(1) && palpha==0) palpha =[];

  if (is_void(data2d)) data2d = array(1,dimsof(R));
  dim_data = dimsof(data2d);
  dim_R   = dimsof(R);
  dim_H   = dimsof(H);
  if (dim_data(1) != dim_R(1) || dim_data(1) != dim_H(1) ||
      anyof(dim_data != dim_R-[0,1,1,1]) ||
       anyof(dim_data != dim_H-[0,1,1,1]) )
    error, "data2d must be a (n-1,n-1) 2D array if R and H is a (n,n) array";


  //==================Surfaces========================
  xyz_surf = data_surf = [];
  for (i=1 ; i<=numberof(pz) ;i++) {
    x_surf = R(,pz(i),) * cos(alpha(,pz(i),));
    y_surf = R(,pz(i),) * sin(alpha(,pz(i),));
    z_surf = H(,pz(i),);
    xyz_surf  = grow(xyz_surf, mesh2polygon(x_surf ,y_surf ,z_surf));
    data_surf = grow(data_surf , data2d(,pz(i),)(*));
  }
  //=================Suround==========================
xyz_suround = data_suround = [];
 for (i=1; i<=numberof(pr) ;i++) {
  x_suround = R(pr(i),,) *  cos(alpha(pr(i),,));
  y_suround = R(pr(i),,) *  sin(alpha(pr(i),,));
  z_suround = H(pr(i),,);
  xyz_suround = grow(xyz_suround, mesh2polygon(x_suround, y_suround, z_suround ));
  data_suround = grow(data_suround , data2d(pr(i),,)(*));
 }
 //======================Cut=================
 xyz_cut = data_cut = [];
  for (i=1 ; i<=numberof(palpha); i++) {
    x_cut = R(,,palpha(i))  * cos(alpha(,,palpha(i)));
    y_cut = R(,,palpha(i))  * sin(alpha(,,palpha(i)));
    z_cut = H(,,palpha(i));
    xyz_cut = grow(xyz_cut,  mesh2polygon(x_cut,y_cut,z_cut));
    data_cut = grow( data_cut, data2d(,,palpha(i))(*));
  }

  //================grow all===================
    xyz  = _(xyz_surf ,xyz_suround , xyz_cut);
    data2d = _(data_surf,data_suround, data_cut);

    return xyz;
}


func mesh2polygon (x,y,z) {
N = numberof(x(1,1,));
xyz = [];
for (i=1; i<=N; i++) {
ncorner =4;
tmp_xyz =   transpose([x(,,i),y(,,i),z(,,i)],2);
ni    = dimsof(x)(2);
nj    = dimsof(x)(3);
list  = indgen(1:ni-1)+ni*indgen(0:nj-2)(-,);
_xyz  = tmp_xyz(,([0,1,ni+1,ni]+list(-,))(*));
tmp_xyz = array(double,3,ncorner,(ni-1)*(nj-1));
tmp_xyz(,*) = _xyz;
xyz = grow(xyz, tmp_xyz);
}
  return xyz;
}


func interfero(model, incl, dist=, hor=, color=, type=, width=, compare_em=,title=,lamb=,compare_scatt=,champ=,gauss=,rt=) {
/* DOCUMENT
   interfero(model, incl, dist=, hor=, color=, type=, width=, compare_em=,title=,lamb=,compare_scatt=,champ=,gauss=,rt=)
   calcule visibilite d'un model mcfost suivant l'inclinaison (axe vertical)
   pix_size en sec
   lambda en microns
   dist en pc
   par defaut coupe verticale sauf si hor=1
   SEE ALSO:
 */

  if ( nameof(structof(model)) != "McfostImage") {
    print, "Error, McfostImage needed in func interfero" ;
    return [] ;
  }

  if (is_void(color)) color="black";
  if (is_void(type)) type=1;
  if (is_void(width)) width=1;
  if (is_void(compare_em)) compare_em=0;
  if (is_void(title)) title=1;
  if (is_void(compare_scatt)) compare_scatt=0;
  lambda = model.lamb;
  if (is_void(lamb)) lamb=lambda;
  if (!is_void(dist)) model.P.map.dist = dist;
  if (is_void(gauss)) gauss=0;
  if (is_void(rt)) rt=0;

  // Field of view
  write, "Field of view of the map = ", (model.P.map.size_neb*2.0)  / (model.P.map.zoom * model.P.map.dist) * 1000, "mas";

  // taille pixel image en AU
  pix_size_AU = (model.P.map.size_neb*2.0) / (model.P.map.ny*model.P.map.zoom);
  pix_size = pix_size_AU / model.P.map.dist ;      // taille pixel image en sec

  ou = 1;
  if (!is_void(champ)) {
    size=model.P.map.ny;
    x=indgen(size)(,-:1:size) - (size/2+1);
    y=indgen(size)(-:1:size,) - (size/2+1);
    distance = abs(x,y);

    ou = (distance * pix_size < 0.5 * champ) ;

    if (gauss==1) {
      FWHM = champ / pix_size; // OK : le FWHM est equivalent a la largeur de la porte !!! Cool !
      sigma = FWHM / (2*sqrt(2*log(2))); // en sec
      ou = gauss_kernel(size, sigma) ;
    }
    write, "Applying a field of view of ", champ, "as" ;
  }

  // images
  if (rt) {
    im = model.image_rt(,,incl,1,1) * ou;
    image = model.image_rt * ou(,,-);
  } else {
    im = model.image(,,incl,1,1) * ou;
    image = model.image * ou(,,-);
  }

  if (dimsof(image)(0) > 5) {
    im_scatt = image(,,incl,6,1) ;
    im_em = image(,,incl,7,1) ;
    im_em_scatt = image(,,incl,8,1) ;
    im_star =  image(,,incl,5,1) ;
    im_star(where(im_star < 1.0e-18)) = 0.0;
  }


  //  stat, im ;
  //stat, im_star + im_scatt + im_em + im_em_scatt; // ---> OK

  // limite 50 mas


  /*  window, 0; fma;
  palette, "heat.gp";
  pli, log10(im+1.0e-30);
  limits;
  */

  // fft
  fim=(fft(im,1));
  if (dimsof(image)(0) > 4) {
    fim_scatt = (fft(im_scatt, 1));
    fim_em = (fft(im_em, 1));
    fim_em_scatt = (fft(im_em_scatt,1));
    fim_star = (fft(im_star,1));
  }


  //stat, abs(fim) ;
  //stat, abs(fim_star + fim_scatt + fim_em + fim_em_scatt); // --> OK


  /*
  window, 10; fma;
  pli, log10((abs(roll(fim)))+1.0e-30);
  limits;
  */


  // lignes de base (en m)
  size = dimsof(im)(2);
  center = size/2+1;

  pix_size = pix_size/3600. * pi/180.; // taille pixel en radian

  pix_fft = 1.0/pix_size; // taille de l'image uv (pas du pixel !!) en unite de lambda
  pix=lambda*1e-6*pix_fft;
  B=span(0,pix/2,size/2); // m

  // visibilités
  if (hor==1) {
    vis = (abs(fim))(1:size/2,1);
    if (dimsof(image)(0) > 4) {
      vis_em = (abs(fim_em))(1:size/2,1);
      vis_scatt = (abs(fim_scatt))(1:size/2,1);
      vis_em_scatt = (abs(fim_em_scatt))(1:size/2,1);
      vis_star = (abs(fim_star))(1:size/2,1);
    }
  } else {
    vis = (abs(fim))(1,1:size/2);
    if (dimsof(image)(0) > 5) {
      vis_em = (abs(fim_em))(1,1:size/2);
      vis_scatt = (abs(fim_scatt))(1,1:size/2);
      vis_em_scatt = (abs(fim_em_scatt))(1,1:size/2);
      vis_star = (abs(fim_star))(1,1:size/2);
    }
  }



  // plot visibilités
  if (compare_scatt) {
    plg, vis/vis(1), B;
    plg, vis_em/vis(1), B, color="red";
    plg, vis_scatt/vis(1), B, color="blue";
    plg, vis_em_scatt/vis(1), B, color="green";
    plg, vis_star/vis(1), B, color="cyan";
    //    plg, (vis_em+vis_em_scatt+vis_star+vis_scatt)/vis(1), B, color="magenta"; // ! debile
    //    plg, ((abs(fim_star + fim_scatt + fim_em + fim_em_scatt))(1,1:size/2)) / vis(1), B, color="cyan"

    xytitles, "B", "V"  ;
    limits ;
    limits, , , 0., 1.;
    return
  }

  // V^2
  //  window, 2 ;
  //plg, (vis)^2, B, color=color, type=type, width=width;
  plg, (vis/vis(1))^2, B, color=color, type=type, width=width;

  if (dimsof(lamb)(1) > 0) {
    n_lamb = dimsof(lamb)(2);
    B_wide = B * lamb(-,)/lambda;
    vis_wide = array(double,dimsof(B)(2),n_lamb);

    for (i=1 ; i<=n_lamb ; i++) {
      vis_wide(,i) = interp(vis, B_wide(,i), B);
    }
    vis =  vis_wide(,avg);

    plg, (vis/vis(1))^2, B, color="green", type=type, width=width;
  }

  if (compare_em==1) {
    // visibilité en supposant que tous l'excès vient de l'émission thermique
    // Attention il faut prendre tout l'excès dans le champ de vue de l'interféromètre !!!
    if (hor==1) {
      vis_pb = (abs(fim_star + fim_em * sum(im_em+im_scatt+im_em_scatt)/sum(im_em)))(1:size/2,1);
    } else {
      vis_pb = (abs(fim_star + fim_em * sum(im_em+im_scatt+im_em_scatt)/sum(im_em)))(1,1:size/2);
    }

    plg, (vis_pb/vis_pb(1))^2, B, color="red", width=width, type=type+1;
  }

  if (compare_em==2) {
    // visibilité en ne prenant que la contribution de l'emission thermique (sans rescaling)
    // donc ON CHANGE L ERAPPORT DE FLUX !!!
    // Attention il faut prendre tout l'excès dans le champ de vue de l'interféromètre !!!
    if (hor==1) {
      vis_pb = (abs(fim_star + fim_em))(1:size/2,1);
    } else {
      vis_pb = (abs(fim_star + fim_em))(1,1:size/2);
    }

    plg, (vis_pb/vis_pb(1))^2, B, color="red", width=width, type=type+1;
  }


  if (title) xytitles, "Baseline (m)", "V^2^"

  limits ;
  limits, , , 0., 1.;

}


func calc_vis(model, incl, &baseline, &vis, dist=, hor=, lamb=, average=, rt=, Jy=, klambda=, Mlambda=, champ=, gauss=,n_padding=) {
/* DOCUMENT : calc_vis(model, incl, &baseline, &vis, dist=, hor=, lamb=, average=, rt=, Jy=, klambda=, Mlambda=, champ=, gauss=, n_padding=)
     calcule visibilite d'un model mcfost suivant l'inclinaison (axe vertical)
     pix_size en sec
     lambda en microns
     dist en pc
     par defaut coupe verticale sauf si hor=1
     vis et baseline en W.m^-2 et m
     vis en Jy avec option Jy=1
     baseline en klambda avec option klambda=1
     baseline en Mlambda avec option Mlambda=1

   SEE ALSO:
 */

  if ( nameof(structof(model)) != "McfostImage") {
    print, "Error, McfostImage needed in func interfero";
    return [];
  }


  lambda = model.lamb;
  if (is_void(lamb)) lamb=lambda;
  if (is_void(rt)) rt=0;
  if (is_void(Jy)) Jy=0;
  if (is_void(klambda)) klambda=0;

  // taille pixel image en AU
  pix_size = (model.P.map.size_neb*2.0) / (model.P.map.ny*model.P.map.zoom);
  //model.P.map.size_neb

  // taille pixel image en sec
  if (is_void(dist)) dist=model.P.map.dist ;
  pix_size = pix_size / dist;
  im_size = (model.P.map.size_neb*2.0) / (model.P.map.zoom * dist);

  ou = 1;
  if (!is_void(champ)) {
    size=model.P.map.ny;
    x=indgen(size)(,-:1:size) - (size/2+1);
    y=indgen(size)(-:1:size,) - (size/2+1);
    distance = abs(x,y);

    ou = (distance * pix_size < 0.5 * champ) ;

    if (gauss==1) {
      FWHM = champ / pix_size; // OK : le FWHM est equivalent a la largeur de la porte !!! Cool !
      sigma = FWHM / (2*sqrt(2*log(2))); // en sec
      ou = gauss_kernel(size, sigma) ;
    }
    // write, "Applying a field of view of ", champ, "as" ;

    if (champ > 0.5 * im_size) {
      write, "WARNING : image seems small to aply the filed of view accurately" ;
      write, im_size, champ ;
    }
  }

  // images
  if (rt) {
    im = model.image_rt(,,incl,1,1) * ou;
  } else {
    im = model.image(,,incl,1,1) * ou;
  }

  // images
  if (rt==1) {
    im = model.image_rt(,,incl,1,1);
  } else {
    im = model.image(,,incl,1,1);
  }

  // speed up Fourier transform by selecting optimized size ;
  if (is_void(n_padding)) n_padding=0 ;
  n = dimsof(im)(2) ;
  n_fft = fft_best_size(n+ 2*n_padding) ; n_extra = (n_fft -n)/2 ;

  write, "Padding with", n_extra,"extra pixels for fft" ;

  im2 = array(float,n_fft,n_fft) ;
  im2(n_extra+1:-n_extra,n_extra+1:-n_extra) = im ;
  im = im2 ; im2 = [] ;

  // fft
  fim=(fft(roll(im),1));

  // lignes de base
  size = n_fft ;
  center = size/2+1;

  pix_size = pix_size/3600. * pi/180.;
  pix_fft = 1.0/pix_size;
  pix=lambda*1e-6*pix_fft;
  baseline=span(0,pix/2,size/2);

  // visibilités
  if (hor==1) {
    vis = fim(1:size/2,1);
  } else {
    vis = fim(1,1:size/2);
  }

  if (average==1) {
    vis = 0.5 * (fim(1:size/2,1) + fim(1,1:size/2));
  }

  if (Jy==1) {
    vis = vis / (SI.c * 1e-26 / (lamb * 1e-6)) ;
  }

  if (klambda==1) {
    baseline = baseline / (lamb * 1e-3) ;
  }

  if (Mlambda==1) {
    baseline = baseline / (lamb * 1e-6) ;
  }

  return fim / (SI.c * 1e-26 / (lamb * 1e-6));

}

func plot_contrib(model,incl,win=,sm=,dist=,cumul=,width=) {
/* DOCUMENT : plot_contrib(model,incl,win=,sm=,dist=,cumul=)

   SEE ALSO:
 */

  if ( nameof(structof(model)) != "McfostImage") {
    print, "Error, McfostImage needed in func plot_contrib";
    return [];
  }

  if (is_void(width)) width=1


  if (!is_void(win)) {
    window, win;
  } else {
    win=1
  }

  // taille pixel image en AU
  pix_size = (model.P.map.size_neb*2.0) / (model.P.map.ny*model.P.map.zoom);
  // taille pixel image en sec
  if (is_void(dist)) dist=model.P.map.dist
  pix_size = pix_size / dist


  size = dimsof(model.image)(2);
  center = size/2+1;

  // abscisse en arcsec
  n_pix = size -center +1;
  x=(span(0,(n_pix)*pix_size,n_pix));

  if (is_void(sm)) sm=0;
  if (is_void(cumul)) cumul=0;

  start=4;
  if (model.P.simu.lsepar_pola == "F")  start=1;

  if (!cumul) {
    plh, smooth((model.image(center:0,center,incl,1,1)),sm), x, width=width;
    plg, smooth((model.image(center:0,center,incl,1,start+1)),sm), x, color="blue", width=width;
    plg, smooth((model.image(center:0,center,incl,1,start+1)),sm), x, color="red", width=width;
    plg, smooth((model.image(center:0,center,incl,1,start+1)),sm), x, color="green", width=width;
  } else {
    frac_etoile = model.image(center,center,incl,1,1)/ sum(model.image(,,incl,1,1));
    maxi = (smooth((model.image(center:0,center,incl,1,1)),sm)*x^2)(sum) ;

    ou = where( (smooth((model.image(center:0,center,incl,1,1)),sm)*x^2)(psum)/maxi * (1.0 - frac_etoile) >  1e-3* frac_etoile(-:1:n_pix) ) ;
    not_ou = where( (smooth((model.image(center:0,center,incl,1,1)),sm)*x^2)(psum)/maxi * (1.0 - frac_etoile) <  1e-3* frac_etoile(-:1:n_pix) ) ;

    y_tmp = ((model.image(center:0,center,incl,1,1)*x^2)(psum)/maxi * (1.0 - frac_etoile) + frac_etoile(-:1:n_pix)) ;
    y = ((smooth((model.image(center:0,center,incl,1,1)),sm)*x^2)(psum)/maxi * (1.0 - frac_etoile) + frac_etoile(-:1:n_pix)) ;


    plg, y(ou) , x(ou) , width=width;
    plg, y(not_ou) , x(not_ou) , width=width;

    // segments vertical et horizontal
    y2 = [y(max(not_ou)),y(min(ou))];
    x2 = [x(min(ou)),x(min(ou))];
    plg, y2, x2, width=width;

    y2 = [y(max(not_ou)),y(max(not_ou))];
    x2 = [x(max(not_ou)),x(min(ou))];
    plg, y2, x2, width=width;

    y_tmp = ((model.image(center:0,center,incl,1,start+1))*x^2)(psum)/maxi * (1.0 - frac_etoile);
    y = (smooth((model.image(center:0,center,incl,1,start+1)),sm)*x^2)(psum)/maxi * (1.0 - frac_etoile);
    //   y(1:10) = y_tmp(1:10) ;

    plg, y, x, color="blue", type=4, width=width;
    plg, (smooth((model.image(center:0,center,incl,1,start+2)),sm)*x^2)(psum)/maxi * (1.0 - frac_etoile), x, color="red",type=3, width=width;
    plg, (smooth((model.image(center:0,center,incl,1,start+3)),sm)*x^2)(psum)/maxi * (1.0 - frac_etoile), x, color="green", type=2, width=width;
    plg, frac_etoile(-:1:n_pix), x, color="magenta", type=5, width=width;
    //(smooth((model.image(center:0,center,incl,1,1)),sm)*x^2)(psum)/maxi * (1.0 - frac_etoile)
  }

  logxy, 0, 1;
  limits;


  /*  window, win+1;
  fma;
  palette, "heat.gp";
  gs_nm, win+1, 2, 3;
  plsys, 1;
  pli, log10(model.image(,,incl,1,1)+1.0e-30);
  limits; logxy, 0, 0;
  plsys, 2;
  pli, log10(model.image(,,incl,1,5)+1.0e-30);
  limits;  logxy, 0, 0;
  plsys, 3;
  pli, log10(model.image(,,incl,1,6)+1.0e-30);
  limits ;  logxy, 0, 0;
  plsys, 4;
  pli, log10(model.image(,,incl,1,7)+1.0e-30);
  limits;  logxy, 0, 0;
  */

}


func plot_pola(model,incl,&pola,&theta,win=,subsample=,scale=,color=,width=,superpose=,rt=,beam=) {
/* DOCUMENT : plot_pola(model,incl,win=,subsample=,scale=,color=,width=,superpose=,rt=)

   SEE ALSO:
 */
  if ( nameof(structof(model)) != "McfostImage") {
    print, "Error, McfostImage needed in func plot_pola";
    return [];
  }

  if (is_void(rt)) rt=0

  if (rt) {
    im=model.image_rt ;
  } else {
    im=model.image ;
  }

  if (!is_void(win)) {
    window, win;
  } else {
    win=1;
  }
  if (is_void(subsample)) subsample=1;
  if (is_void(scale)) scale=subsample;
  if (is_void(superpose)) superpose=0;
  if (is_void(width)) width=1;
  if (is_void(color)) color="black";

  nx = model.P.map.nx;
  ny = model.P.map.ny;

  x=(indgen(nx))(,-:1:ny) - 0.5;
  y=(indgen(ny))(-:1:nx,) -0.5;

  I0 = im(,,incl,1,1,1)+1.0e-30 ;


  if (!is_void(beam)) {
    I = convol2df(im(,,incl,1,1),beam) ;
    Q = convol2df(im(,,incl,1,2),beam) ;
    U = convol2df(im(,,incl,1,3),beam) ;
  } else {
    I = im(,,incl,1,1) ;
    Q = im(,,incl,1,2) ;
    U = im(,,incl,1,3) ;
  }

  Iu = I ;
  Qu = Q ;
  Uu = U ;

  Iu = undersample2d(I,subsample)+1.0e-30;
  Qu = undersample2d(Q,subsample);
  Uu = undersample2d(U,subsample);

  x=undersample2d(x,subsample);
  y=undersample2d(y,subsample);

  pola = sqrt((Qu/Iu)^2+(Uu/Iu)^2);
  theta = 0.5*atan(Uu,Qu) ;

  pola_x = - pola * sin(theta); // Ref is N (vertical axis) --> sin,  and Est is toward left --> -
  pola_y = pola * cos(theta);

  // Correction pola dans mcfost
  pola_x = pola * -cos(theta);
  pola_y = pola * sin(theta);

  write, "Max pola [%] =", max(pola)*100 ;

  if (superpose) pli, (I0+1.0e-30)^0.1,  0.5, 0.5, nx+0.5 , ny+0.5;
  pldj, x-0.5*scale*pola_x, y+0.5*scale*pola_y,x+0.5*scale*pola_x, y-0.5*scale*pola_y, color=color, width=width;

  /*
  window, 21 ; fma ;
  window, 21 ; gs_nm, 21, 2, 3, square=1 ;
  plsys, 1 ;
  pli, pola_x ;
  plsys, 2 ;
  pli, pola_y ;
  plsys, 3 ;
  pli, pola ;
  plsys, 4 ;
  pli, theta ;
  plsys, 5 ;
  */
  //pli, pola_x/pola_y ;

  //stat, pola_x/pola_y ;

  return pola ;
}

//------------

func plot_sed(model,incl,av=,rv=,contrib=,color=,width=,scale=,sed1=,w=,file=,type=,icontrib=,ou=,rt=) {
/* DOCUMENT plot_sed(model,incl,av=,rv=,contrib=,color=,width=,scale=,sed1=,rt=)
   rv = 3.1, 4.0 ou 5.5
   SEE ALSO:
 */
  if (is_void(icontrib)) icontrib = 1;
  if (is_void(rt)) rt=0;
  if (is_void(rv)) rv=3.1;
  if (is_void(av)) av=0.0;
  if (is_void(contrib)) contrib=0;
  if (is_void(color)) color="black";
  if (is_void(width)) width=1;
  if (is_void(scale)) scale=1;
  if (is_void(sed1)) sed1=0;
  if (is_void(w)) w=0;
  if (!is_void(file)) w=1 ;
  if (is_void(type)) type=1;


  srv = swrite(rv,format="%3.1f") ;

   if ( nameof(structof(model)) != "McfostSED") {
     print, "Error, McfostSED needed in func plot_sed";
     return [];
   }

   //if ((max(model.sed) < 1.0e-60) && (rt==0)) {
   if ((max(model.sed) < 1.0e-60) ) {
     lamb = model.lamb1 ;
     sed1 = 1;
     "Warning : sed2 = 0, plotting sed1";
   } else {
     lamb = model.lamb ;
   }
   if (is_void(ou)) ou = where(lamb > 0.);


   if (av > 0.) {
     extinction = OpenASCII(get_cwd()+"/SED/kext_albedo_WD_MW_"+srv+"_D03.all",prompt=0) ;

     ext_V = abs(extinction.lamb - 5.47e-01)(mnx) ;
     kext = extinction.kpa / (1.0-extinction.albedo) ;

     XX =  kext /  kext(ext_V) ;
     XXm = interp( XX , extinction.lamb ,  lamb) ;

     tau_V=0.4*log(10)*av ;
     correct_Av = exp(-tau_V * XXm);
   } else {
     correct_Av=1.0;
   }

   if (model.P.simu.lsepar_pola=="T") {
     istart=4 ;
   } else {
     istart=1 ;
   }

   if (sed1) {
     SED = model.sed1 ;
   } else if (rt) {
     SED = model.sed_rt ;
     print, "i=", model.P.map.RT_incl(incl) ;
   }  else {
     SED = model.sed ;
     istart = 4 ;
   }

   ou = ou(sort(lamb(ou))) ;

   plg, (SED(,incl,1,icontrib)*correct_Av*scale)(ou), lamb(ou), color=color,width=width, type=type;
   if (contrib) {
     order = sort(model.lamb) ;
     pla, (SED(order,incl,1,istart+[1,2,3,4])*correct_Av*scale), model.lamb(order), type=[2,3,4,5],  color=["magenta","blue","red","green"],width=width;
   }

   if (w) {
     if (!is_void(file)) {
       f = open(file,"w") ;
       write, f, "# lambda (micron), lambda.F_lambda (W.m^-2)";
       write, f, model.lamb, SED(ou,incl,1,icontrib)*correct_Av*scale;
       close, f ;
     } else {
       write, "# lambda (micron), lambda.F_lambda (W.m^-2)";
       write, model.lamb, SED(ou,incl,1,icontrib)*correct_Av*scale;
     }
   }


   return ;
}

func sum_mcfost(model) {
/* DOCUMENT sum_mcfost(model)

   SEE ALSO:
 */

  if (dimsof(model)(1) != 1) {
    print, "Error : can only sum 1D array of models";
    return [];
  } else {
    n_mod = dimsof(model)(2);
  }

  if ( nameof(structof(model)) == "McfostSED") {
    mod = model(1);
    mod.sed = model.sed(..,avg);
    mod.sed1 = model.sed1(..,avg);

    // Pour la temperature, on moyenne les energies
    mie, mod.P.dust_pop(1,1).amin, mod.P.dust_pop(1,1).amax, mod.P.dust_pop(1,1).aexp, mod.P.dust_pop(1,1).n_grains, mod.lamb1, "~/mcfost_utils/"+mod.P.dust_pop(1,1).file,1,kappa,albedo,S11,S12,g; // Calcul des kappa

    kappa_abs = kappa * (1.0-albedo) ;

    // Calcul du corps noir
    //cn = bb(model.T,mod.lamb1);
    // Integration (en log) de kappa*B + moyennage sur les modeles
    //E_abs = (kappa*mod.lamb1 * cn())(sum,..,avg);

    // Boucle a la main pour economiser mem
    E_abs = mod.T*0.; //allocation
    // Somme sur les modeles
    for (n=1 ; n<= n_mod ; n++) {
      write, "modele", n, "/", n_mod;
      somme=0.0;
      // Integration (en log) de kappa*B + moyennage sur les modeles
      for (l=1 ; l<=mod.P.wave.nlamb1; l++) {
        cn = bb(model(n).T,mod.lamb1(l))(1,);
        somme = somme + kappa_abs(l)*mod.lamb1(l)*cn;
      }
      E_abs = E_abs + somme;
    }
    E_abs = E_abs/n_mod;


    // Pretabulation des Energie en fonction des temperatures
    Temp = spanl(mod.P.gridT.T_min,mod.P.gridT.T_max,mod.P.gridT.n_T);
    cn = bb(Temp,mod.lamb1);
    E_em =  (kappa_abs*mod.lamb1 * cn())(sum,..);

    mod.T = interp(Temp, E_em, E_abs);
    //    mod.T = ((model.T^4)(..,avg))^0.25; // Pas propre
  } else if ( nameof(structof(model)) == "McfostImage") {
    mod=model(1);
    mod.image=model.image(..,avg);
  } else {
    print, "Error, Mcfost model needed in func sum_mcfost";
    return [];
  }

  return mod;
}



func open_sum_mcfost(dir,type) {

/* DOCUMENT open_sum_mcfost(dir,type) {
     Lit et somme progressivement une serie de modeles mcfost
     en economisant la memoire
     Ne marche que pour les images pour le moment
   SEE ALSO:
 */
  n = numberof(dir);

  1 ;
  mod = open_mcfost(dir(1),type);
  tmp_im = mod.image

  for (i=2 ; i<=n ; i++) {
    i ;
    tmp_mod = open_mcfost(dir(i),type);
    tmp_im = tmp_im + tmp_mod.image ;
  }
  mod.image = tmp_im / n ;
  return mod

}

func plot_az(mod,i) {

  I = 2.5 * log(mod.image(sum,sum,i,,1)/max(mod.image(sum,sum,i,,1))) ;
  P = sqrt(test.image(sum,sum,i,,2)^2+test.image(sum,sum,i,,3)^2)/test.image(sum,sum,i,,1)* 100.;

  fma;
  plg, I;
  plg, P, color="red";

  print, "attenuation =",  2.5 * log(max(mod.image(sum,sum,1,,1))/max(mod.image(sum,sum,i,,1)))
  print, "Delta I (mag) =", - min(I);
  print, "Pola", min(P), max(P) - min(P);
}

func plot_line(model,itrans,i,insert=,substract_cont=) {

  if (is_void(insert)) insert=1 ;
  if (is_void(substract_cont)) substract_cont = 0 ;

  if (is_void(i)) {
    "ERROR: argument missing" ;
    return[] ;
  }

  vmax = model.P.mol.vmax(1) ;
  v = span(-vmax,vmax,2*model.P.mol.n_speed(1)+1) ;
  dv = vmax / n_speed * SI.km ; // m/s

  inu = where(itrans==model.ifreq) ;

  if (!numberof(inu)) {
    "Error : transition not computed" ;
    return  ;
  } else {
    inu = inu(1) ;
  }


  nu = model.freq(inu) ;
  lambda = SI.c/nu ;
  dlambda = dv / c * lambda ;

  get_style, landscape,systems,legends,clegends ;
  systems = systems(-:1:2);
  systems.viewport(,2) = [0.2,0.35,0.7,0.85];
  set_style, landscape,systems,legends,clegends ;

  plsys, 1;

  Fline = model.Fline(,,,inu,i) ;

  F_lambda = Fline / lambda ;
  F_nu = Fline / nu * 1e26 ;

  //"Fline:"
  //info, Fline ;
  //write, "max Fline:", max(Fline(sum,sum,)) ;

  Fcont =  model.Fcont(,,inu,i)

  Fc_lambda = Fcont / lambda ;
  Fc_nu = Fcont / nu * 1e26 ;

  //"Fcont:"
  //info, Fcont ;
  //write, SI.c/model.freq(inu), "Fcont:", model.Fcont(sum,sum,inu,i), "W.m^-2" ;
  //sum(Fcont) ;

  if (substract_cont) {
    plg, F_nu(sum,sum,) - Fc_nu(sum,sum), v, color="red" ;
  } else {
    plg, F_nu(sum,sum,), v, color="red" ;
    plg, [Fc_nu(sum,sum), Fc_nu(sum,sum)], [-vmax,vmax], type=2;
  }

  write, "Fline = ",  sum([F_nu(sum,sum,) - Fc_nu(sum,sum)]) * dv / SI.km, "Jy.km/s"

  limits ; relimits ;
  l = limits() ;
  if (l(4) < 1.01 * l(3)) {
    l(3) = 0.99 * l(3) ;
    l(4) = 1.01 * l(4) ;
    limits, l ;
  }
  xytitles, "v (km/s)", "F_!n (Jy)" ;

  //pltr, name, 60, 90, height=12 ;
  Fl = sum(F_lambda(,sum,sum) - Fc_lambda(sum,sum)) * dlambda ;

  pltr, model.mol_name+"   !l = "+swrite(lambda * 1e6,format="%3.3f")+" !mm", 60, 90, height=12 ;
  pltr, "F_l_ = "+swrite(Fl,format="%3.2e")+" W.m^-2^", 60, 80, height=12 ;
  pltr, "F_c_ = "+swrite(Fc_nu(sum,sum),format="%3.2e")+" Jy", 60, 70, height=12 ;

  if (insert) {
    if ((dimsof(F_nu))(0) > 1) {
      plsys, 2 ; gridxy,0x200;
      plim, (F_nu(,,avg) - Fc_nu(,))  ;
      limits ;
    }
  }

  if (substract_cont) {
    return [F_nu(sum,sum,) - Fc_nu(sum,sum), v] ;
  } else {
    return [F_nu(sum,sum,), v] ;
  }

}


func L2R(L,Teff) {
  /* DOCUMENT  L2R(log10(Lstar/Lsun),Teff) renvoie Rstar/Rsun
   */
  return sqrt(10^L * (5777./Teff)^4) ;
}

func R2L(R,Teff) {
  /* DOCUMENT  R2L(Rstar/Rsun,Teff) renvoie Lstar/Lsun
   */
  return R^2 * (Teff/5777.)^4 ;
}


func mmVis(model,i,champ=,gauss=,bin=,start=) {

/* DOCUMENT
   calcule visibilite d'un modele mcfost suivant l'inclinaison dans des bins de largeur donnes
   moyenne azimuthale
   renvoie visibilite dans les bins en Jy
   SEE ALSO:
 */

  if ( nameof(structof(model)) != "McfostImage") {
    print, "Error, McfostImage needed in func mmVis" ;
    return [] ;
  }


  lambda = model.lamb ;


  if (is_void(color)) color="black";
  if (is_void(type)) type=1;
  if (is_void(width)) width=1;
  if (!is_void(dist)) model.P.map.dist = dist;
  if (is_void(gauss)) gauss=0;
  if (is_void(bin)) bin=10.;
  if (is_void(r_start)) r_start=0.;

  // taille pixel image en AU
  pix_size_AU = (model.P.map.size_neb*2.0) / (model.P.map.ny*model.P.map.zoom);
  // taille pixel image en sec
  pix_size = pix_size_AU / model.P.map.dist ;

  //  print, "pix_size", pix_size;

  ou = 1;
  if (!is_void(champ)) {
    size=model.P.map.ny;
    x=indgen(size)(,-:1:size) - (size/2+1);
    y=indgen(size)(-:1:size,) - (size/2+1);
    distance = abs(x,y);

    ou = (distance * pix_size < 0.5 * champ) ;

    if (gauss==1) {
      FWHM = champ / pix_size; // OK : le FWHM est equivalent a la largeur de la porte !!! Cool !
      sigma = FWHM / (2*sqrt(2*log(2))); // en sec
      ou = gauss_kernel(size, sigma) ;
    }
  }

  // images
  im = model.image(,,i,1,1) * ou;

  // fft
  fim=abs(roll(fft(im,1)));

  //window, 20 ; fma ; plim, fim;
  //pause, 3000;


  // lignes de base en m
  size = dimsof(im)(2);
  center = size/2+1;

  pix_size = pix_size/3600. * pi/180.;
  pix_fft = 1.0/pix_size;
  pix=lambda*1e-6*pix_fft;
  B=span(0,pix/2,size/2);

  pix = pix / (size/2) ;

  // Surechantillonage image


  // Distance des pixels a l'etoile
  tmp1 = indgen(nx)(,-:1:ny) - center;
  tmp2 = indgen(ny)(-:1:nx,) - center;
  pixel_map = [tmp1, tmp2];
  star_pos = [center,center];
  vector_map = pixel_map;
  distance = abs(vector_map(,,1),vector_map(,,2));
  distance = distance * pix;

  nbin = int(max(distance)/sqrt(2.)/bin);

  V = array(double,nbin);

  for (i=1 ; i<=nbin ; i++) {
    rmin = r_start + (i-1) * bin;
    rmax = r_start +  i   * bin;

    ou = where( (distance >= rmin) * (distance < rmax) ) ;

    //    plim, fim * (distance >= rmin) * (distance < rmax);
    //pause, 1000;

    if (numberof(ou)) {
      V(i) = avg(fim(ou)) ;
    } else {
      V(i) = 1e-20;
    }
  }

  // En Jy
  //V = V * lambda*1e-6 / c * 1.0e26;

  window, 0;
  plg, fim(center,center:center+size/2-1),  B;
  plg, fim(center:center+size/2-1,center),  B;
  plg, 0.5*(fim(center:center+size/2-1,center) +  fim(center,center:center+size/2-1)),  B;

  plp, V, (indgen(numberof(V)) - 0.5) * bin;

  return V;


}

func indice_spectral(mod,incl,lambda1,lambda2) {
/* DOCUMENT indice_spectral(mod,incl,lambda1,lambda2)

   SEE ALSO:
 */
  lFnu = log(mod.lamb*mod.sed(,incl,1)) ;
  llamb=log(mod.lamb) ;

  // C'est bizarre c'est pas du nu.Fnu
  f1 = interp(lFnu,llamb,log(lambda1));
  f2 = interp(lFnu,llamb,log(lambda2));

  return -((f1) - (f2) ) / (log(lambda1) - log(lambda2));
}


//**********************************************************

func integ_mass_size_distrib(amin,amax,p) {
  return(amax^(p+4) - amin^(p+4)) / (p+4)  ;
}

//**********************************************************

func surface_density(dens,grid,plot=,color=,smooth_length=) {
/* DOCUMENT
     Calcul et plot la densite de surface d'un modele mcfost
     a partir de density.fits.gz et grid.fits.gz
   SEE ALSO:
 */

  if (is_void(plot)) plot=0 ;
  if (is_void(color)) color="black" ;

  dz = (grid(,2,1,2) - grid(,1,1,2)) * SI.AU / SI.cm ; // en cm

  Sigma = 2 * dens(,sum) * dz ;
  R = grid(,1, 1) ;

  if (plot) {
    if (is_void(smooth_length)) {
      plg, Sigma, R, color=color ;
    } else {
      plg, smooth(Sigma,smooth_length), R, color=color ;
    }
    //plp, Sigma, R, color=color ;
    xytitles, "R [AU]", "!S_dust_ [g.cm^-2^]", [0.02,0.02] ;
    limits ;
    limits, , , 0.5 * min(Sigma(where(Sigma > 1e-30))), 2* max(Sigma) ;
  }

  //write, "R(AU)   Sigma (g.cm-2)" ;
  //write, R, Sigma ;

  //write, "M", sum( (Sigma * 2*pi * R)(zcen) * R(dif)) * (SI.AU/SI.cm)^2 * SI.g/ SI.Msun ;

  return [R,Sigma] ;
}

/* ************************************* */

func IsFile(name, restrict)
/* DOCUMENT yocoTypeIsFile(name, restrict)

   DESCRIPTION
     Return an array of int of the same dimension of name.

   PARAMETERS
     - name: name to be tested
     - restrict (optional):
       If not given (or 0) :
        - 1 if name is an existing file
        - 2 if name is an existing directory
        - 0 else
       If set to 1, speed up the search by avoiding a directory listing:
        - 1 if name is an existing file or directory
        - 0 else

   EXAMPLES
     > yocoTypeIsFile("~/");
     2
     > yocoTypeIsFile("~/",1);
     1
     > yocoTypeIsFile( ["./","test","/home","~/.emacs"] );
     [2,0,0,1]
*/
{
    local type, i, flag;

    type = array(int,dimsof(name));

    for(i=1;i<=numberof(name);i++) {
        /* test is file exist */
        if(open(name(i),"r",1)) type(i) += 1;
        /* add 1 if name is a directory */
        if(type(i)==1 & !restrict) {
            flag = lsdir(name(i));
            if(!(dimsof(flag)(1)==0 && flag==0))  type(i) += 1;
        }
    }
    return type;
}


//**********************************************************

/* // Myriam's version
func computeV2New(model,u,v,incli,PA,&CV)
{
  incl = incli;
  lambda = model.lamb;
  dist = model.P.map.dist;
  //write, "Field of view = ", (model.P.map.size_neb*2.0)  / (model.P.map.zoom * model.P.map.dist) * 1000, "mas";
  //write, "zoom =", model.P.map.zoom ;


  pix_size_AU = (model.P.map.size_neb*2.0) / (model.P.map.ny*model.P.map.zoom);
  pix_size = pix_size_AU / model.P.map.dist ;
  resx = resy = pix_size_AU / model.P.map.dist;
  //print, "pix size = "+pr1(resx);

  if (dimsof(model.image)(0) > 4) {
    im_scatt = model.image(,,incl,5,1) ;
    im_em = model.image(,,incl,6,1) ;
    im_em_scatt = model.image(,,incl,7,1) ;
    im_star = im - im_scatt - im_em - im_em_scatt;
    im_star(where(im_star < 1.0e-18)) = 0.0;
  }

  im = model.image_rt(,,incl,1,1) ;

  Pa = PA-90; //si mcfost a une convention de PA defini par 1/2 petit axe
  PAr = Pa*pi/180;
  U = u*sin(PAr) + v*cos(PAr);
  V = -u*cos(PAr) + v*sin(PAr);

  Nx = dimsof(im)(2);
  Ny = dimsof(im)(3);

  //res = lambda/B; Bmax = B/2;
  res= sqrt(resx^2 + resy^2); // in arcsec
  Bmax = (lambda/res*4.8)/2;

  resx_1surD = resx*4.8/lambda; //in m-1 freq in B/lambda units, inverse of arcsec
  resy_1surD = resy*4.8/lambda; //in m-1

  alpha_x = U*resx_1surD;
  alpha_y = V*resy_1surD;

  matx_tmp = indgen(1:Nx);
  maty_tmp = indgen(1:Ny);

  //create a 2-D array similar to:
  matx = maty = array(double,Nx,Ny);
  for (i=1;i<=Nx;i++){
    matx(i,) = matx_tmp;
  }
  for (i=1;i<=Ny;i++){
    maty(,i) = maty_tmp;
  }

  FT = im(,,-)*exp(-2*1i*pi*(maty(,,-)*alpha_x(-,-,) + matx(,,-)*alpha_y(-,-,)));
  CV = FT(sum,sum,)/sum(im);
  VIS2 = (CV.re)^2 + (CV.im)^2;

  return VIS2 ;

}

//**********************************************************

func computeT3(model,u1,v1,u2,v2,incli,PA)
{

  U1 = u1; V1 = v1;
  U2 = u2; V2 = v2;
  U3 = -u1-u2; V3 = -v1-v2;

  U123 = V123 = array(double,3,numberof(U1));

  U123(1,) = U1;  U123(2,) = U2;  U123(3,) = U3;
  V123(1,) = V1;  V123(2,) = V2;  V123(3,) = V3;

  VC = array(complex,3,numberof(U1));

  for (i=1;i<=3;i++){
    tmp = computeV2New(model,U123(i,),V123(i,),incli,PA,CV);
    VC(i,) = CV;
  }

  BS = VC(1,)*VC(2,)*VC(3,);
  CP = atan(BS.im, BS.re) * 180/pi;

  return CP;

}
*/

//**********************************************************

func computeV2New(model,u,v,incli,PA,&CV)
{
  incl = incli;
  lambda = model.lamb;
  dist = model.P.map.dist;


  pix_size_AU = (model.P.map.size_neb*2.0) / (model.P.map.ny*model.P.map.zoom);
  pix_size = pix_size_AU / model.P.map.dist ;
  resx = resy = pix_size_AU / model.P.map.dist;

  if (dimsof(model.image)(0) > 4) {
    im_scatt = model.image(,,incl,5,1) ;
    im_em = model.image(,,incl,6,1) ;
    im_em_scatt = model.image(,,incl,7,1) ;
    im_star = im - im_scatt - im_em - im_em_scatt;
    im_star(where(im_star < 1.0e-18)) = 0.0;
  }

  im = model.image_rt(,,incl,1,1) ;

  //*******************************

  fim=(roll(fft(roll(im),1)))/sum(im) ;

  pix_size_AU = (model.P.map.size_neb*2.0) / (model.P.map.ny*model.P.map.zoom);
  pix_size = pix_size_AU / model.P.map.dist ;      // taille pixel image en sec

  // lignes de base (en m)
  size = dimsof(im)(2);
  center = size/2+1;

  pix_size = pix_size/3600. * pi/180.;

  pix_fft = 1.0/pix_size;
  pix=lambda*1e-6*pix_fft;
  B = span(-1,1,size) * pix/2. ;

  Uplane = B(,-:1:size) ; // signe - car Est a gauche dans convention astro
  Vplane = B(-:1:size,) ;

   /*
  Pa = PA - 90; // mcfost a une convention de PA defini par 1/2 petit axe
  PAr = Pa*pi/180 ;
  U = u*sin(PAr) + v*cos(PAr); // pourqoi ?? ca revient a faire 2 fois 90 degrees
  V = -u*cos(PAr) + v*sin(PAr);
  */

  PAr = PA*pi/180 ;
  U = u*cos(PAr) + v*sin(PAr) ;
  V = v*cos(PAr) - u*sin(PAr) ;



  Vre = interp2(U,V,  fim.re,Uplane,Vplane) ;
  Vim = interp2(U,V,  fim.im,Uplane,Vplane) ;

  //Vre = spline2(fim.re, U/max(B) *size,V/max(B)*size) ;
  //Vim = spline2(fim.im, U/max(B) *size,V/max(B)*size) ;

  CV = Vre + 1i * Vim ;
  VIS2 = abs(CV)^2 ;

  return VIS2

}

//**********************************************************

func computeT3(model,u1,v1,u2,v2,incli,PA)
{

  U1 = u1; V1 = v1;
  U2 = u2; V2 = v2;
  U3 = -u1-u2; V3 = -v1-v2;

  U123 = V123 = array(double,3,numberof(U1));

  U123(1,) = U1;  U123(2,) = U2;  U123(3,) = U3;
  V123(1,) = V1;  V123(2,) = V2;  V123(3,) = V3;

  VC = array(complex,3,numberof(U1));

  for (i=1;i<=3;i++){
    tmp = computeV2New(model,U123(i,),V123(i,),incli,PA,CV);
    VC(i,) = CV;
  }

  BS = VC(1,)*VC(2,)*VC(3,);
  CP = atan(BS.im, BS.re) * 180/pi;

  return CP;

}


func partial_mass(model,R1,R2) {

/* DOCUMENT partial_mass(model,R1,R2)
     compute the dust mass [Msun] of a given mcfost model between the radii R1 and R2 [AU]
   SEE ALSO:
 */

  if (R2 < R1) {
    R_tmp = R2 ;
    R2 = R1 ;
    R1 = R_tmp ;
  }

  P = model.P ;
  M_tot = 0. ;
  for (i=1 ; i<= P.simu.n_zones ; i++) {
    mass = P.zones(i).mass ;
    rin = P.zones(i).rin ;
    rout = P.zones(i).rout ;
    p = P.zones(i).dens ;

    if ( ((R1 >= rin) && (R1 <= rout)) ||  ((R2 >= rin) & (R2 <= rout)) ) {
      Mzone =  mass * (min(rout,R2)^(p+1) - max(rin,R1)^(p+1)) / (rout^(p+1) - rin^(p+1)) ;
    } else {
      Mzone = 0. ;
    }

    write, "zone", i, "mass=", Mzone ;
    M_tot += Mzone ;
  }
  return M_tot ;

}

write, format=" MCFOST-Yorick %3.2f package loaded\n", MCFOST_version ;


func ProDiMo_star_2_mcfost_star(prodimo_file,mcfost_file) {

  star = OpenASCII(prodimo_file,noheader=1,valtype=float,prompt=0) ;

  // Flux W.m^2 pour 1 Rsun a 1pc
  Flux = 4*pi * star.X2 * SI.c/(star.X1*1e-9) * 1e-3 * (SI.Rsun)^2 / (SI.pc^2) ;
  wl = star.X1*1e-3 ;

  n = numberof(wl) ;

  S = array(float,n,3) ;
  S(,1) = wl ;
  S(,2) = Flux/wl ;

  cfitsioWrite, mcfost_file, S ;
}


func mcfost(P, options=,dir=,dev=,clean=) {
/* DOCUMENT  mcfost(P, options=,dir=,dev=)

   SEE ALSO:
 */

  if (is_void(clean)) clean=0 ;
  if (is_void(dev)) dev=0 ;

  if (dev==1) {
    write, "Using MCFOST development version" ;
    mcfost_binary = "~/mcfost/src/mcfost " ;
  } else {
    mcfost_binary = "~/Software/mcfost/bin/mcfost " ;
  }

  if (is_void(dir)) {
    if (clean) system, "rm -rf data_*" ;
    system, mcfost_binary+P+" "+options ;
  } else {
    dir0 = get_cwd() ;
    cd, dir ;
    if (clean) system, "rm -rf data_*" ;
    mcfost_binary+P+" "+options ;
    system, mcfost_binary+P+" "+options ;
    cd, dir0 ;
  }
}


func T2d_to_T3d(T,naz) {

  D = dimsof(T) ;
  nx = D(2) ; nz = D(3) ;

  T3D = array(float,nx,2*nz+1,naz) ;

  for (j=1 ; j<=nz ; j++) {
    T3D(:,nz+1+j,:) = T(:,j) ;
    T3D(:,nz+1-j,:) = T(:,j) ;
  }
  return T3D ;
}


func pressure(gas_density,temperature) {
/* DOCUMENT pressure(gas_density,temperature)
   gas_density en g/cm3
   Temperature en K
   can use mcfost fits files

   SEE ALSO:
 */

  masseH = 1./SI.Na ;
  masse_mol = 2.3 * masseH ; // masse molecule de gaz (en g)

  rho = gas_density / masse_mol / SI.cm^3 ; // density in molecules per m3

  P = rho * SI.kB * T ; // Pa

  return P ;
}
