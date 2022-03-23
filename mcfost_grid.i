require, "init/mcfost_struct.i" ;

func create_grid(name,reference_file,
                 diskmass=,H0=,Rin=,Rout=,edge=,beta=,surf=,amin=,amax=,aexp=, gas_dust=, viscosity=,
                 diskmass_zone=,H0_zone=,Rin_zone=,Rout_zone=,edge_zone=,beta_zone=,surf_zone=,amin_zone=,amax_zone=,aexp_zone=,
                 diskmass_2=,H0_2=,Rin_2=,Rout_2=,edge_2=,beta_2=,surf_2=,amin_2=,amax_2=,aexp_2=,
                 diskmass_zone_2=,H0_zone_2=,Rin_zone_2=,Rout_zone_2=,edge_zone_2=,beta_zone_2=,surf_zone_2=,amin_zone_2=,amax_zone_2=,aexp_zone_2=,
                 diskmass_3=,H0_3=,Rin_3=,Rout_3=,edge_3=,beta_3=,surf_3=,amin_3=,amax_3=,aexp_3=,
                 diskmass_zone_3=,H0_zone_3=,Rin_zone_3=,Rout_zone_3=,edge_zone_3=,beta_zone_3=,surf_zone_3=,amin_zone_3=,amax_zone_3=,aexp_zone_3=,
                 diskmass_4=,H0_4=,Rin_4=,Rout_4=,edge_4=,beta_4=,surf_4=,amin_4=,amax_4=,aexp_4=,
                 diskmass_zone_4=,H0_zone_4=,Rin_zone_4=,Rout_zone_4=,edge_zone_4=,beta_zone_4=,surf_zone_4=,amin_zone_4=,amax_zone_4=,aexp_zone_4=,
                 amin_specie=,amax_specie=,aexp_specie=,
                 amin_specie_2=,amax_specie_2=,aexp_specie_2=,
                 amin_specie_3=,amax_specie_3=,aexp_specie_3=,
                 amin_specie_4=,amax_specie_4=,aexp_specie_4=,
                 Rstar=,Tstar=,fUV=,slope_UV=,star=,strat=,R0Rin=,correct_Rsub=,istart=,zone=,gasp=,basename=,extension=,return_param=) {
  /* DOCUMENT create_grid(name,reference_file,diskmass=,H0=,Rin=,Rout=,edge=,beta=,surf=,Rstar=,Tstar=,star=,amin=,amax=,aexp=,strat=,R0Rin=,istart=,zone=,gasp=,basename=,extension=,return_param=,
   diskmass_zone=,H0_zone=,Rin_zone=,Rout_zone=,edge_zone=,beta_zone=,surf_zone=,amin_zone=,amax_zone=,aexp_zone=)

   name : name of the grid of models to calculate
   reference_file : mcfost reference parameter file

   Parameters that can be varied :
      - add  _#  for zone #
      - add extra keword   PARAMETER_zone(_#) for the zones it applies to (ie in case of correlated zones)
      - add extra keword   PARAMETER_specie(_#) in case of multiple species

   diskmass=
   H0=
   Rin=
   Rout=
   edge=
   beta=
   surf=
   correct_Rsub=

   amin=
   amax=
   aexp=
   strat=

   Tstar=
   Rstar=
   star=

   R0Rin=1 -->  H0 defines the ratio H0(Rin)/Rin instead of H0(R0)

   Command line options : default = none --> SED calculation
   cmd_opt=


   SEE ALSO:
   _create_grid
*/

  extern n_params_tot;

  //tic ;

  V = 2.20 ;

  if (is_void(istart)) istart = 0 ;
  if (is_void(zone)) zone = 1 ;
  if (is_void(R0Rin)) R0Rin = 0 ;
  if (is_void(gasp)) gasp = 0 ;
  if (is_void(extension)) extension = ".par" ;
  if (is_void(return_param)) return_param=0 ;


  if (is_void(diskmass_zone)) diskmass_zone=zone ;
  if (is_void(H0_zone)) H0_zone=zone ;
  if (is_void(Rin_zone)) Rin_zone=zone ;
  if (is_void(Rout_zone)) Rout_zone=zone ;
  if (is_void(edge_zone)) edge_zone=zone ;
  if (is_void(beta_zone)) beta_zone=zone ;
  if (is_void(surf_zone)) surf_zone=zone ;
  if (is_void(amin_zone)) amin_zone=zone ;
  if (is_void(amax_zone)) amax_zone=zone ;
  if (is_void(aexp_zone)) aexp_zone=zone ;
  if (is_void(amin_specie)) amin_specie=1 ;
  if (is_void(amax_specie)) amax_specie=1 ;
  if (is_void(aexp_specie)) aexp_specie=1 ;

  if (is_void(diskmass_zone_2)) diskmass_zone_2=2 ;
  if (is_void(H0_zone_2)) H0_zone_2=2 ;
  if (is_void(Rin_zone_2)) Rin_zone_2=2 ;
  if (is_void(Rout_zone_2)) Rout_zone_2=2 ;
  if (is_void(edge_zone_2)) edge_zone_2=2 ;
  if (is_void(beta_zone_2)) beta_zone_2=2 ;
  if (is_void(surf_zone_2)) surf_zone_2=2 ;
  if (is_void(amin_zone_2)) amin_zone_2=2 ;
  if (is_void(amax_zone_2)) amax_zone_2=2 ;
  if (is_void(aexp_zone_2)) aexp_zone_2=2 ;
  if (is_void(amin_specie_2)) amin_specie_2=1 ;
  if (is_void(amax_specie_2)) amax_specie_2=1 ;
  if (is_void(aexp_specie_2)) aexp_specie_2=1 ;


  if (is_void(diskmass_zone_3)) diskmass_zone_3=3 ;
  if (is_void(H0_zone_3)) H0_zone_3=3 ;
  if (is_void(Rin_zone_3)) Rin_zone_3=3 ;
  if (is_void(Rout_zone_3)) Rout_zone_3=3 ;
  if (is_void(edge_zone_3)) edge_zone_3=3 ;
  if (is_void(beta_zone_3)) beta_zone_3=3 ;
  if (is_void(surf_zone_3)) surf_zone_3=3 ;
  if (is_void(amin_zone_3)) amin_zone_3=3 ;
  if (is_void(amax_zone_3)) amax_zone_3=3 ;
  if (is_void(aexp_zone_3)) aexp_zone_3=3 ;
  if (is_void(amin_specie_3)) amin_specie_3=1 ;
  if (is_void(amax_specie_3)) amax_specie_3=1 ;
  if (is_void(aexp_specie_3)) aexp_specie_3=1 ;


  if (is_void(diskmass_zone_4)) diskmass_zone_4=4 ;
  if (is_void(H0_zone_4)) H0_zone_4=4 ;
  if (is_void(Rin_zone_4)) Rin_zone_4=4 ;
  if (is_void(Rout_zone_4)) Rout_zone_4=4 ;
  if (is_void(edge_zone_4)) edge_zone_4=4 ;
  if (is_void(beta_zone_4)) beta_zone_4=4 ;
  if (is_void(surf_zone_4)) surf_zone_4=4 ;
  if (is_void(amin_zone_4)) amin_zone_4=4 ;
  if (is_void(amax_zone_4)) amax_zone_4=4 ;
  if (is_void(aexp_zone_4)) aexp_zone_4=4 ;
  if (is_void(amin_specie_4)) amin_specie_4=1 ;
  if (is_void(amax_specie_4)) amax_specie_4=1 ;
  if (is_void(aexp_specie_4)) aexp_specie_4=1 ;


  // Lecture du fichier de ref mcfost
  ref = read_params(reference_file) ;

  // Mise a jour des structures
  include, struct_file ;

  if (abs(V -ref.simu.version) > 1.e-4) {
    write, "GRID: Version ", V, "needed. Exiting !!!";
    return;
  }


  // Structure de la grille
  mc_grid = array(MCFOST_grid) ;
  mc_grid.ref = ref ;

  // Definition grille de structure de modeles
  dims = array(1,n_params_tot);

  // Parametres des zones
  if (!is_void(diskmass)) dims(1)=dimsof(diskmass)(0); mc_grid.names(1) = "diskmass" ;
  if (!is_void(H0)) dims(2)=dimsof(H0)(0);  mc_grid.names(2) = "H0" ;
  if (!is_void(R0)) dims(3)=dimsof(R0)(0);  mc_grid.names(3) = "R0" ;
  if (!is_void(Rin)) dims(4)=dimsof(Rin)(0);  mc_grid.names(4) = "Rin" ;
  if (!is_void(Rout)) dims(5)=dimsof(Rout)(0);  mc_grid.names(5) = "Rout" ;
  if (!is_void(edge)) dims(6)=dimsof(edge)(0); mc_grid.names(6) = "edge" ;
  if (!is_void(beta)) dims(7)=dimsof(beta)(0); mc_grid.names(7) = "beta" ;
  if (!is_void(surf)) dims(8)=dimsof(surf)(0); mc_grid.names(8) = "surf" ;

  if (!is_void(amin)) dims(9)=dimsof(amin)(0); mc_grid.names(9) = "amin" ;
  if (!is_void(amax))dims(10)=dimsof(amax)(0); mc_grid.names(10) = "amax" ;
  if (!is_void(aexp)) dims(11)=dimsof(aexp)(0); mc_grid.names(11) = "aexp" ;

  I0 = 11 ;
  if (!is_void(diskmass_2)) dims(I0+1)=dimsof(diskmass_2)(0); mc_grid.names(I0+1) = "diskmass_2" ;
  if (!is_void(H0_2)) dims(I0+2)=dimsof(H0_2)(0);  mc_grid.names(I0+2) = "H0_2" ;
  if (!is_void(R0_2)) dims(I0+3)=dimsof(R0_2)(0);  mc_grid.names(I0+3) = "R0_2" ;
  if (!is_void(Rin_2)) dims(I0+4)=dimsof(Rin_2)(0);  mc_grid.names(I0+4) = "Rin_2" ;
  if (!is_void(Rout_2)) dims(I0+5)=dimsof(Rout_2)(0);  mc_grid.names(I0+5) = "Rout_2" ;
  if (!is_void(edge_2)) dims(I0+6)=dimsof(edge_2)(0); mc_grid.names(I0+6) = "edge_2" ;
  if (!is_void(beta_2)) dims(I0+7)=dimsof(beta_2)(0); mc_grid.names(I0+7) = "beta_2" ;
  if (!is_void(surf_2)) dims(I0+8)=dimsof(surf_2)(0); mc_grid.names(I0+8) = "surf_2" ;

  if (!is_void(amin_2)) dims(I0+9)=dimsof(amin_2)(0); mc_grid.names(I0+9) = "amin_2" ;
  if (!is_void(amax_2))dims(I0+10)=dimsof(amax_2)(0); mc_grid.names(I0+10) = "amax_2" ;
  if (!is_void(aexp_2)) dims(I0+11)=dimsof(aexp_2)(0); mc_grid.names(I0+11) = "aexp_2" ;


  if (!is_void(diskmass_3)) dims(2*I0+1)=dimsof(diskmass_3)(0); mc_grid.names(2*I0+1) = "diskmass_3" ;
  if (!is_void(H0_3)) dims(2*I0+2)=dimsof(H0_3)(0);  mc_grid.names(2*I0+2) = "H0_3" ;
  if (!is_void(R0_3)) dims(2*I0+3)=dimsof(R0_3)(0);  mc_grid.names(2*I0+3) = "R0_3" ;
  if (!is_void(Rin_3)) dims(2*I0+4)=dimsof(Rin_3)(0);  mc_grid.names(2*I0+4) = "Rin_3" ;
  if (!is_void(Rout_3)) dims(2*I0+5)=dimsof(Rout_3)(0);  mc_grid.names(2*I0+5) = "Rout_3" ;
  if (!is_void(edge_3)) dims(2*I0+6)=dimsof(edge_3)(0); mc_grid.names(2*I0+6) = "edge_3" ;
  if (!is_void(beta_3)) dims(2*I0+7)=dimsof(beta_3)(0); mc_grid.names(2*I0+7) = "beta_3" ;
  if (!is_void(surf_3)) dims(2*I0+8)=dimsof(surf_3)(0); mc_grid.names(2*I0+8) = "surf_3" ;

  if (!is_void(amin_3)) dims(2*I0+9)=dimsof(amin_3)(0); mc_grid.names(2*I0+9) = "amin_3" ;
  if (!is_void(amax_3))dims(2*I0+10)=dimsof(amax_3)(0); mc_grid.names(2*I0+10) = "amax_3" ;
  if (!is_void(aexp_3)) dims(2*I0+11)=dimsof(aexp_3)(0); mc_grid.names(2*I0+11) = "aexp_3" ;


  if (!is_void(diskmass_4)) dims(3*I0+1)=dimsof(diskmass_4)(0); mc_grid.names(3*I0+1) = "diskmass_4" ;
  if (!is_void(H0_4)) dims(3*I0+2)=dimsof(H0_4)(0);  mc_grid.names(3*I0+2) = "H0_4" ;
  if (!is_void(R0_4)) dims(3*I0+3)=dimsof(R0_4)(0);  mc_grid.names(3*I0+3) = "R0_4" ;
  if (!is_void(Rin_4)) dims(3*I0+4)=dimsof(Rin_4)(0);  mc_grid.names(3*I0+4) = "Rin_4" ;
  if (!is_void(Rout_4)) dims(3*I0+5)=dimsof(Rout_4)(0);  mc_grid.names(3*I0+5) = "Rout_4" ;
  if (!is_void(edge_4)) dims(3*I0+6)=dimsof(edge_4)(0); mc_grid.names(3*I0+6) = "edge_4" ;
  if (!is_void(beta_4)) dims(3*I0+7)=dimsof(beta_4)(0); mc_grid.names(3*I0+7) = "beta_4" ;
  if (!is_void(surf_4)) dims(3*I0+8)=dimsof(surf_4)(0); mc_grid.names(3*I0+8) = "surf_4" ;

  if (!is_void(amin_4)) dims(3*I0+9)=dimsof(amin_4)(0); mc_grid.names(3*I0+9) = "amin_4" ;
  if (!is_void(amax_4))dims(3*I0+10)=dimsof(amax_4)(0); mc_grid.names(3*I0+10) = "amax_4" ;
  if (!is_void(aexp_4)) dims(3*I0+11)=dimsof(aexp_4)(0); mc_grid.names(3*I0+11) = "aexp_4" ;

  // Parametres communs a toutes les zones
  if (!is_void(strat)) dims(4*I0+1)=dimsof(strat)(0); mc_grid.names(4*I0+1) = "strat" ;
  if (!is_void(Rstar)) dims(4*I0+2)=dimsof(Rstar)(0); mc_grid.names(4*I0+2) = "Rstar" ;
  if (!is_void(Tstar)) dims(4*I0+3)=dimsof(Tstar)(0); mc_grid.names(4*I0+3) = "Tstar" ;
  if (!is_void(fUV)) dims(4*I0+4)=dimsof(fUV)(0); mc_grid.names(4*I0+4) = "fUV" ;
  if (!is_void(slope_UV)) dims(4*I0+5)=dimsof(slope_UV)(0); mc_grid.names(4*I0+5) = "slope_UV" ;
  if (!is_void(star)) dims(4*I0+6)=dimsof(star)(0); mc_grid.names(4*I0+6) = "star" ;

  if (!is_void(correct_Rsub)) dims(4*I0+7)=dimsof(correct_Rsub)(0); mc_grid.names(4*I0+7) = "correct_Rsub" ;

  if (!is_void(viscosity)) dims(4*I0+8)=dimsof(viscosity)(0); mc_grid.names(4*I0+8) = "viscosity" ;

  if (!is_void(gas_dust)) dims(4*I0+9)=dimsof(gas_dust)(0); mc_grid.names(4*I0+9) = "gas_dust" ;

  mc_grid.dims = dims ;

  // Indices des parametres varies
  ou = where(dims > 1);
  n_params_explores = numberof(ou);
  dims = dims(ou);

  n_models = 1;
  for (i=1 ; i<=n_params_explores ; i++) {
    n_models *= dims(i);
  }

  parameters = ref(1);

  n_zones = parameters.simu.n_zones ;


  write, "Ordre des parametres:" ;

  // diskmass
  if (!is_void(diskmass)) {
    write, "Masse", numberof(diskmass);
    mc_grid.values(1,1:numberof(diskmass))= diskmass ;
    parameters = parameters(..,-:1:numberof(diskmass)) ;
    for (i=1 ; i<= numberof(diskmass) ; i++) {
      for (z=1 ; z<=numberof(diskmass_zone) ; z++) {
        parameters.zones(diskmass_zone(z),..).mass(..,i) = diskmass(i);
      }
    }
  }

  if (!is_void(diskmass_2)) {
    write, "Masse zone 2", numberof(diskmass_2);
    mc_grid.values(I0+1,1:numberof(diskmass_2))= diskmass_2 ;   // TODO : index !!! et pour tous les suivants
    parameters = parameters(..,-:1:numberof(diskmass_2)) ;
    for (i=1 ; i<= numberof(diskmass_2) ; i++) {
      for (z=1 ; z<=numberof(diskmass_zone_2) ; z++) {
        parameters.zones(diskmass_zone_2(z),..).mass(..,i) = diskmass_2(i);
      }
    }
  }

  if (!is_void(diskmass_3)) {
    write, "Masse zone 3", numberof(diskmass_3);
    mc_grid.values(2*I0+1,1:numberof(diskmass_3))= diskmass_3 ;
    parameters = parameters(..,-:1:numberof(diskmass_3)) ;
    for (i=1 ; i<= numberof(diskmass_3) ; i++) {
      for (z=1 ; z<=numberof(diskmass_zone_3) ; z++) {
        parameters.zones(diskmass_zone_3(z),..).mass(..,i) = diskmass_3(i);
      }
    }
  }

  if (!is_void(diskmass_4)) {
    write, "Masse zone 4", numberof(diskmass_4);
    mc_grid.values(3*I0+1,1:numberof(diskmass_4))= diskmass_4 ;
    parameters = parameters(..,-:1:numberof(diskmass_4)) ;
    for (i=1 ; i<= numberof(diskmass_4) ; i++) {
      for (z=1 ; z<=numberof(diskmass_zone_4) ; z++) {
        parameters.zones(diskmass_zone_4(z),..).mass(..,i) = diskmass_4(i);
      }
    }
  }


  // H0
  if (!is_void(H0)) {
    write, "H0", numberof(H0);
    mc_grid.values(2,1:numberof(H0))= H0 ;
    parameters = parameters(..,-:1:numberof(H0)) ;
    for (i=1 ; i<= numberof(H0) ; i++) {
      for (z=1 ; z<=numberof(H0_zone) ; z++) {
        parameters.zones(H0_zone(z),..).Ho(..,i) = H0(i);
      }
    }
  }

  if (!is_void(H0_2)) {
    write, "H0 zone 2", numberof(H0_2);
    mc_grid.values(I0+2,1:numberof(H0_2))= H0_2 ;
    parameters = parameters(..,-:1:numberof(H0_2)) ;
    for (i=1 ; i<= numberof(H0_2) ; i++) {
      for (z=1 ; z<=numberof(H0_zone_2) ; z++) {
        parameters.zones(H0_zone_2(z),..).Ho(..,i) = H0_2(i);
      }
    }
  }

  if (!is_void(H0_3)) {
    write, "H0 zone 3", numberof(H0_3);
    mc_grid.values(2*I0+2,1:numberof(H0_3))= H0_3 ;
    parameters = parameters(..,-:1:numberof(H0_3)) ;
    for (i=1 ; i<= numberof(H0_3) ; i++) {
      for (z=1 ; z<=numberof(H0_zone_3) ; z++) {
        parameters.zones(H0_zone_3(z),..).Ho(..,i) = H0_3(i);
      }
    }
  }

  if (!is_void(H0_4)) {
    write, "H0 zone 4", numberof(H0_4);
    mc_grid.values(3*I0+2,1:numberof(H0_4))= H0_4 ;
    parameters = parameters(..,-:1:numberof(H0_4)) ;
    for (i=1 ; i<= numberof(H0_4) ; i++) {
      for (z=1 ; z<=numberof(H0_zone_4) ; z++) {
        parameters.zones(H0_zone_4(z),..).Ho(..,i) = H0_4(i);
      }
    }
  }




  // Rin
  if (!is_void(Rin)) {
    write, "Rin", numberof(Rin) ;
    mc_grid.values(4,1:numberof(Rin))= Rin ;
    parameters = parameters(..,-:1:numberof(Rin)) ;
    for (i=1 ; i<= numberof(Rin) ; i++) {
      for (z=1 ; z<=numberof(Rin_zone) ; z++) {
        parameters.zones(Rin_zone(z),..).rin(..,i) = Rin(i);
      }

      if (R0Rin) { // echelle de hauteur definie en h/R a Rin
        for (z=1 ; z<=numberof(Rin_zone) ; z++) {
          parameters.zones(Rin_zone(z),..).Ro(..,i) = Rin(i);
          parameters.zones(Rin_zone(z),..).Ho(..,i) = parameters.zones(Rin_zone,..).Ho(..,i) * Rin(i) ;
        }
      }
    }
  }


  if (!is_void(Rin_2)) {
    write, "Rin zone 2", numberof(Rin_2) ;
    mc_grid.values(I0+4,1:numberof(Rin_2))= Rin_2 ;
    parameters = parameters(..,-:1:numberof(Rin_2)) ;
    for (i=1 ; i<= numberof(Rin_2) ; i++) {
      for (z=1 ; z<=numberof(Rin_zone_2) ; z++) {
        parameters.zones(Rin_zone_2(z),..).rin(..,i) = Rin_2(i);
      }

      if (R0Rin) { // echelle de hauteur definie en h/R a Rin
        for (z=1 ; z<=numberof(Rin_zone_2) ; z++) {
          parameters.zones(Rin_zone_2(z),..).Ro(..,i) = Rin_2(i);
          parameters.zones(Rin_zone_2(z),..).Ho(..,i) = parameters.zones(Rin_zone_2,..).Ho(..,i) * Rin_2(i) ;
        }
      }
    }
  }


  if (!is_void(Rin_3)) {
    write, "Rin zone 3", numberof(Rin_3) ;
    mc_grid.values(2*I0+4,1:numberof(Rin_3))= Rin_3 ;
    parameters = parameters(..,-:1:numberof(Rin_3)) ;
    for (i=1 ; i<= numberof(Rin_3) ; i++) {
      for (z=1 ; z<=numberof(Rin_zone_3) ; z++) {
        parameters.zones(Rin_zone_3(z),..).rin(..,i) = Rin_3(i);
      }

      if (R0Rin) { // echelle de hauteur definie en h/R a Rin
        for (z=1 ; z<=numberof(Rin_zone_3) ; z++) {
          parameters.zones(Rin_zone_3(z),..).Ro(..,i) = Rin_3(i);
          parameters.zones(Rin_zone_3(z),..).Ho(..,i) = parameters.zones(Rin_zone_3,..).Ho(..,i) * Rin_3(i) ;
        }
      }
    }
  }

  if (!is_void(Rin_4)) {
    write, "Rin zone 4", numberof(Rin_4) ;
    mc_grid.values(3*I0+4,1:numberof(Rin_4))= Rin_4 ;
    parameters = parameters(..,-:1:numberof(Rin_4)) ;
    for (i=1 ; i<= numberof(Rin_4) ; i++) {
      for (z=1 ; z<=numberof(Rin_zone_4) ; z++) {
        parameters.zones(Rin_zone_4(z),..).rin(..,i) = Rin_4(i);
      }

      if (R0Rin) { // echelle de hauteur definie en h/R a Rin
        for (z=1 ; z<=numberof(Rin_zone_4) ; z++) {
          parameters.zones(Rin_zone_4(z),..).Ro(..,i) = Rin_4(i);
          parameters.zones(Rin_zone_4(z),..).Ho(..,i) = parameters.zones(Rin_zone_4,..).Ho(..,i) * Rin_4(i) ;
        }
      }
    }
  }

  // Rout
  if (!is_void(Rout)) {
    write, "Rout", numberof(Rout) ;
    mc_grid.values(5,1:numberof(Rout))= Rout ;
    parameters = parameters(..,-:1:numberof(Rout)) ;
    for (i=1 ; i<= numberof(Rout) ; i++) {
      for (z=1 ; z<=numberof(Rout_zone) ; z++) {
        parameters.zones(Rout_zone(z),..).rout(..,i) = Rout(i);
      }
    }
  }

  if (!is_void(Rout_2)) {
    write, "Rout zone 2", numberof(Rout_2) ;
    mc_grid.values(I0+5,1:numberof(Rout_2))= Rout_2 ;
    parameters = parameters(..,-:1:numberof(Rout_2)) ;
    for (i=1 ; i<= numberof(Rout_2) ; i++) {
      for (z=1 ; z<=numberof(Rout_zone_2) ; z++) {
        parameters.zones(Rout_zone_2(z),..).rout(..,i) = Rout_2(i);
      }
    }
  }

  if (!is_void(Rout_3)) {
    write, "Rout zone 3", numberof(Rout_3) ;
    mc_grid.values(2*I0+5,1:numberof(Rout_3))= Rout_3 ;
    parameters = parameters(..,-:1:numberof(Rout_3)) ;
    for (i=1 ; i<= numberof(Rout_3) ; i++) {
      for (z=1 ; z<=numberof(Rout_zone_3) ; z++) {
        parameters.zones(Rout_zone_3(z),..).rout(..,i) = Rout_3(i);
      }
    }
  }

  if (!is_void(Rout_4)) {
    write, "Rout zone 4", numberof(Rout_4) ;
    mc_grid.values(3*I0+5,1:numberof(Rout_4))= Rout_4 ;
    parameters = parameters(..,-:1:numberof(Rout_4)) ;
    for (i=1 ; i<= numberof(Rout_4) ; i++) {
      for (z=1 ; z<=numberof(Rout_zone_4) ; z++) {
        parameters.zones(Rout_zone_4(z),..).rout(..,i) = Rout_4(i);
      }
    }
  }

  //edge
  if (!is_void(edge)) {
    write, "edge", numberof(edge) ;
    mc_grid.values(6,1:numberof(edge))= edge ;
    parameters = parameters(..,-:1:numberof(edge)) ;
    for (i=1 ; i<= numberof(edge) ; i++) {
      for (z=1 ; z<=numberof(edge_zone) ; z++) {
        parameters.zones(edge_zone(z),..).edge(..,i) = edge(i);
      }
    }
  }

  if (!is_void(edge_2)) {
    write, "edge zone 2", numberof(edge_2) ;
    mc_grid.values(I0+6,1:numberof(edge_2))= edge_2 ;
    parameters = parameters(..,-:1:numberof(edge_2)) ;
    for (i=1 ; i<= numberof(edge_2) ; i++) {
      for (z=1 ; z<=numberof(edge_zone_2) ; z++) {
        parameters.zones(edge_zone_2(z),..).edge(..,i) = edge_2(i);
      }
    }
  }

  if (!is_void(edge_3)) {
    write, "edge zone 3", numberof(edge_3) ;
    mc_grid.values(2*I0+6,1:numberof(edge_3))= edge_3 ;
    parameters = parameters(..,-:1:numberof(edge_3)) ;
    for (i=1 ; i<= numberof(edge_3) ; i++) {
      for (z=1 ; z<=numberof(edge_zone_3) ; z++) {
        parameters.zones(edge_zone_3(z),..).edge(..,i) = edge_3(i);
      }
    }
  }

  if (!is_void(edge_4)) {
    write, "edge zone 4", numberof(edge_4) ;
    mc_grid.values(3*I0+6,1:numberof(edge_4))= edge_4 ;
    parameters = parameters(..,-:1:numberof(edge_4)) ;
    for (i=1 ; i<= numberof(edge_4) ; i++) {
      for (z=1 ; z<=numberof(edge_zone_4) ; z++) {
        parameters.zones(edge_zone_4(z),..).edge(..,i) = edge_4(i);
      }
    }
  }


  // beta
  if (!is_void(beta)) {
    write, "beta", numberof(beta) ;
    mc_grid.values(7,1:numberof(beta))= beta ;
    parameters = parameters(..,-:1:numberof(beta)) ;
    for (i=1 ; i<= numberof(beta) ; i++) {
      for (z=1 ; z<=numberof(beta_zone) ; z++) {
        parameters.zones(beta_zone(z),..).beta(..,i) = beta(i);
      }
    }
  }

  if (!is_void(beta_2)) {
    write, "beta zone 2", numberof(beta_2) ;
    mc_grid.values(I0+7,1:numberof(beta_2))= beta_2 ;
    parameters = parameters(..,-:1:numberof(beta_2)) ;
    for (i=1 ; i<= numberof(beta_2) ; i++) {
      for (z=1 ; z<=numberof(beta_zone_2) ; z++) {
        parameters.zones(beta_zone_2(z),..).beta(..,i) = beta_2(i);
      }
    }
  }

  if (!is_void(beta_3)) {
    write, "beta zone 3", numberof(beta_3) ;
    mc_grid.values(2*I0+7,1:numberof(beta_3))= beta_3 ;
    parameters = parameters(..,-:1:numberof(beta_3)) ;
    for (i=1 ; i<= numberof(beta_3) ; i++) {
      for (z=1 ; z<=numberof(beta_zone_3) ; z++) {
        parameters.zones(beta_zone_3(z),..).beta(..,i) = beta_3(i);
      }
    }
  }

  if (!is_void(beta_4)) {
    write, "beta zone 4", numberof(beta_4) ;
    mc_grid.values(3*I0+7,1:numberof(beta_4))= beta_4 ;
    parameters = parameters(..,-:1:numberof(beta_4)) ;
    for (i=1 ; i<= numberof(beta_4) ; i++) {
      for (z=1 ; z<=numberof(beta_zone_4) ; z++) {
        parameters.zones(beta_zone_4(z),..).beta(..,i) = beta_4(i);
      }
    }
  }

  // Surface density
  if (!is_void(surf)) {
    write, "surf", numberof(surf) ;
    mc_grid.values(8,1:numberof(surf))= surf ;
    parameters = parameters(..,-:1:numberof(surf)) ;
    for (i=1 ; i<= numberof(surf) ; i++) {
      for (z=1 ; z<=numberof(surf_zone) ; z++) {
        parameters.zones(surf_zone(z),..).dens(..,i) = surf(i);
      }
    }
  }

  if (!is_void(surf_2)) {
    write, "surf zone 2", numberof(surf_2) ;
    mc_grid.values(I0+8,1:numberof(surf_2))= surf_2 ;
    parameters = parameters(..,-:1:numberof(surf_2)) ;
    for (i=1 ; i<= numberof(surf_2) ; i++) {
      for (z=1 ; z<=numberof(surf_zone_2) ; z++) {
        parameters.zones(surf_zone_2(z),..).dens(..,i) = surf_2(i);
      }
    }
  }

  if (!is_void(surf_3)) {
    write, "surf zone 3", numberof(surf_3) ;
    mc_grid.values(2*I0+8,1:numberof(surf_3))= surf_3 ;
    parameters = parameters(..,-:1:numberof(surf_3)) ;
    for (i=1 ; i<= numberof(surf_3) ; i++) {
      for (z=1 ; z<=numberof(surf_zone_3) ; z++) {
        parameters.zones(surf_zone_3(z),..).dens(..,i) = surf_3(i);
      }
    }
  }

  if (!is_void(surf_4)) {
    write, "surf zone 4", numberof(surf_4) ;
    mc_grid.values(3*I0+8,1:numberof(surf_4))= surf_4 ;
    parameters = parameters(..,-:1:numberof(surf_4)) ;
    for (i=1 ; i<= numberof(surf_4) ; i++) {
      for (z=1 ; z<=numberof(surf_zone_4) ; z++) {
        parameters.zones(surf_zone_4(z),..).dens(..,i) = surf_4(i);
      }
    }
  }


  // amin
  if (!is_void(amin)) {
    write, "amin", numberof(amin) ;
    mc_grid.values(9,1:numberof(amin))= amin ;
    parameters = parameters(..,-:1:numberof(amin)) ;
    for (i=1 ; i<= numberof(amin) ; i++) {
      for (z=1 ; z<=numberof(amin_zone) ; z++) {
        for (s=1 ; s<=numberof(amin_specie) ; s++) {
          parameters.dust_pop.amin(amin_zone(z),amin_specie(s),..,i) = amin(i);
        }
      }
    }
  }


  if (!is_void(amin_2)) {
    write, "amin zone 2", numberof(amin_2) ;
    mc_grid.values(I0+9,1:numberof(amin_2))= amin_2 ;
    parameters = parameters(..,-:1:numberof(amin_2)) ;
    for (i=1 ; i<= numberof(amin_2) ; i++) {
      for (z=1 ; z<=numberof(amin_zone_2) ; z++) {
        for (s=1 ; s<=numberof(amin_specie_2) ; s++) {
          parameters.dust_pop.amin(amin_zone_2(z),amin_specie_2(s),..,i) = amin_2(i);
        }
      }
    }
  }

  if (!is_void(amin_3)) {
    write, "amin zone 3", numberof(amin_3) ;
    mc_grid.values(2*I0+9,1:numberof(amin_3))= amin_3 ;
    parameters = parameters(..,-:1:numberof(amin_3)) ;
    for (i=1 ; i<= numberof(amin_3) ; i++) {
      for (z=1 ; z<=numberof(amin_zone_3) ; z++) {
        for (s=1 ; s<=numberof(amin_specie_3) ; s++) {
          parameters.dust_pop.amin(amin_zone_3(z),amin_specie_3(s),..,i) = amin_3(i);
        }
      }
    }
  }

  if (!is_void(amin_4)) {
    write, "amin zone 4", numberof(amin_4) ;
    mc_grid.values(3*I0+9,1:numberof(amin_4))= amin_4 ;
    parameters = parameters(..,-:1:numberof(amin_4)) ;
    for (i=1 ; i<= numberof(amin_4) ; i++) {
      for (z=1 ; z<=numberof(amin_zone_4) ; z++) {
        for (s=1 ; s<=numberof(amin_specie_4) ; s++) {
          parameters.dust_pop.amin(amin_zone_4(z),amin_specie_4(s),..,i) = amin_4(i);
        }
      }
    }
  }


  // amax
  if (!is_void(amax)) {
    write, "amax", numberof(amax) ;
    mc_grid.values(10,1:numberof(amax))= amax ;
    parameters = parameters(..,-:1:numberof(amax)) ;
    for (i=1 ; i<= numberof(amax) ; i++) {
      for (z=1 ; z<=numberof(amax_zone) ; z++) {
          for (s=1 ; s<=numberof(amax_specie) ; s++) {
            parameters.dust_pop.amax(amax_zone(z),amax_specie(s),..,i) = amax(i);
          }
      }
    }
  }

  if (!is_void(amax_2)) {
    write, "amax zone 2", numberof(amax_2) ;
    mc_grid.values(I0+10,1:numberof(amax_2))= amax_2 ;
    parameters = parameters(..,-:1:numberof(amax_2)) ;
    for (i=1 ; i<= numberof(amax_2) ; i++) {
      for (z=1 ; z<=numberof(amax_zone_2) ; z++) {
          for (s=1 ; s<=numberof(amax_specie_2) ; s++) {
            parameters.dust_pop.amax(amax_zone_2(z),amax_specie_2(s),..,i) = amax_2(i);
          }
      }
    }
  }


  if (!is_void(amax_3)) {
    write, "amax zone 3", numberof(amax_3) ;
    mc_grid.values(2*I0+10,1:numberof(amax_3))= amax_3 ;
    parameters = parameters(..,-:1:numberof(amax_3)) ;
    for (i=1 ; i<= numberof(amax_3) ; i++) {
      for (z=1 ; z<=numberof(amax_zone_3) ; z++) {
        for (s=1 ; s<=numberof(amax_specie_3) ; s++) {
          parameters.dust_pop.amax(amax_zone_3(z),amax_specie_3(s),..,i) = amax_3(i);
        }
      }
    }
  }

  if (!is_void(amax_4)) {
    write, "amax zone 4", numberof(amax_4) ;
    mc_grid.values(3*I0+10,1:numberof(amax_4))= amax_4 ;
    parameters = parameters(..,-:1:numberof(amax_4)) ;
    for (i=1 ; i<= numberof(amax_4) ; i++) {
      for (z=1 ; z<=numberof(amax_zone_4) ; z++) {
        for (s=1 ; s<=numberof(amax_specie_4) ; s++) {
          parameters.dust_pop.amax(amax_zone_4(z),amax_specie_4(s),..,i) = amax_4(i);
        }
      }
    }
  }

  // aexp
  if (!is_void(aexp)) {
    write, "aexp", numberof(aexp) ;
    mc_grid.values(11,1:numberof(aexp))= aexp ;
    parameters = parameters(..,-:1:numberof(aexp)) ;
    for (i=1 ; i<= numberof(aexp) ; i++) {
      for (z=1 ; z<=numberof(aexp_zone) ; z++) {
        for (s=1 ; s<=numberof(aexp_specie) ; s++) {
          parameters.dust_pop.aexp(aexp_zone(z),aexp_specie(s),..,i) = aexp(i);
        }
      }
    }
  }

  if (!is_void(aexp_2)) {
    write, "aexp zone 2", numberof(aexp_2) ;
    mc_grid.values(I0+11,1:numberof(aexp_2))= aexp_2 ;
    parameters = parameters(..,-:1:numberof(aexp_2)) ;
    for (i=1 ; i<= numberof(aexp_2) ; i++) {
      for (z=1 ; z<=numberof(aexp_zone_2) ; z++) {
        for (s=1 ; s<=numberof(aexp_specie_2) ; s++) {
          parameters.dust_pop.aexp(aexp_zone_2(z),aexp_specie_2(s),..,i) = aexp_2(i);
        }
      }
    }
  }


  if (!is_void(aexp_3)) {
    write, "aexp zone 3", numberof(aexp_3) ;
    mc_grid.values(2*I0+11,1:numberof(aexp_3))= aexp_3 ;
    parameters = parameters(..,-:1:numberof(aexp_3)) ;
    for (i=1 ; i<= numberof(aexp_3) ; i++) {
      for (z=1 ; z<=numberof(aexp_zone_3) ; z++) {
        for (s=1 ; s<=numberof(aexp_specie_3) ; s++) {
          parameters.dust_pop.aexp(aexp_zone_3(z),aexp_specie_3(s),..,i) = aexp_3(i);
        }
      }
    }
  }

  if (!is_void(aexp_4)) {
    write, "aexp zone 4", numberof(aexp_4) ;
    mc_grid.values(3*I0+11,1:numberof(aexp_4))= aexp_4 ;
    parameters = parameters(..,-:1:numberof(aexp_4)) ;
    for (i=1 ; i<= numberof(aexp_4) ; i++) {
      for (z=1 ; z<=numberof(aexp_zone_4) ; z++) {
        for (s=1 ; s<=numberof(aexp_specie_4) ; s++) {
          parameters.dust_pop.aexp(aexp_zone_4(z),aexp_specie_4(s),..,i) = aexp_4(i);
        }
      }
    }
  }

  if (!is_void(strat)) {
    write, "strat", numberof(strat) ;
    mc_grid.values(4*I0+1,1:numberof(strat))= strat ;
    parameters = parameters(..,-:1:numberof(strat)) ;
    for (i=1 ; i<= numberof(strat) ; i++) {
      parameters.dust.strat(..,i) = strat(i);
    }
  }

  if (!is_void(Rstar)) {
    write, "Rstar", numberof(Rstar) ;
    mc_grid.values(4*I0+2,1:numberof(Rstar))= Rstar ;
    parameters = parameters(..,-:1:numberof(Rstar)) ;
    for (i=1 ; i<= numberof(Rstar) ; i++) {
      parameters.stars.R(..,i) = Rstar(i);
    }
  }

  if (!is_void(Tstar)) {
    write, "Tstar", numberof(Tstar) ;
    mc_grid.values(4*I0+3,1:numberof(Tstar))= Tstar ;
    parameters = parameters(..,-:1:numberof(Tstar)) ;
    for (i=1 ; i<= numberof(Tstar) ; i++) {
      parameters.stars.T(..,i) = Tstar(i);
    }
  }

  if (!is_void(fUV)) {
    write, "fUV", numberof(fUV) ;
    mc_grid.values(4*I0+4,1:numberof(fUV))= fUV ;
    parameters = parameters(..,-:1:numberof(fUV)) ;
    for (i=1 ; i<= numberof(fUV) ; i++) {
      parameters.stars.fUV(..,i) = fUV(i);
    }
  }

  if (!is_void(slope_UV)) {
    write, "slope_UV", numberof(slope_UV) ;
    mc_grid.values(4*I0+4,1:numberof(slope_UV))= slope_UV ;
    parameters = parameters(..,-:1:numberof(slope_UV)) ;
    for (i=1 ; i<= numberof(slope_UV) ; i++) {
      parameters.stars.slope_UV(..,i) = slope_UV(i);
    }
  }


  if (!is_void(star)) {
    write, "star", numberof(star) ;
    mc_grid.values(4*I0+6,1:numberof(star))= indgen(numberof(star)) ;
    mc_grid.stars(1:numberof(star))= star ;
    parameters = parameters(..,-:1:numberof(star)) ;
    for (i=1 ; i<= numberof(star) ; i++) {
      parameters.stars(..,i) = star(i);
    }
  }

  if (!is_void(correct_Rsub)) {
    write, "correct_Rsub", numberof(correct_Rsub) ;
    mc_grid.values(4*I0+7,1:numberof(correct_Rsub))= correct_Rsub ;
    parameters = parameters(..,-:1:numberof(correct_Rsub)) ;
    for (i=1 ; i<= numberof(correct_Rsub) ; i++) {
      parameters.simu.correct_Rsub(..,i) = correct_Rsub(i);
    }
  }

  if (!is_void(viscosity)) {
    write, "viscosity", numberof(viscosity) ;
    mc_grid.values(4*I0+8,1:numberof(viscosity))= viscosity ;
    parameters = parameters(..,-:1:numberof(viscosity)) ;
    for (i=1 ; i<= numberof(viscosity) ; i++) {
      parameters.simu.viscosity(..,i) = viscosity(i);
    }
  }

    if (!is_void(gas_dust)) {
    write, "gas_dust", numberof(gas_dust) ;
    mc_grid.values(4*I0+8,1:numberof(gas_dust))= gas_dust ;
    parameters = parameters(..,-:1:numberof(gas_dust)) ;
    for (i=1 ; i<= numberof(gas_dust) ; i++) {
      parameters.dust.gas_dust(..,i) = gas_dust(i);
    }
  }


  info, parameters;
  write, numberof(parameters), "models";

  // Creation du tar avec les fichiers de parametres
  if (gasp==1) {
    _create_grid_gasp, name,parameters,basename,extension ;
  } else {
    "Creating parameter files ..." ;
    _create_grid, name,parameters,istart,extension ;
    "Making tar file ..." ;
    system, "tar --directory "+name+" -cf "+get_cwd()+"/"+name+".tar . " ;// --remove-files" ;
    system, "gzip -9 "+name+".tar" ;
    system, "mv "+name+".tar.gz "+name+"_P.tgz" ;
    //system, "rm -rf "+name ;

  }

  // Fichier binaire avec l'info sur la grille
  f = createb(name+".grid");
  save, f, mc_grid ;
  close, f;

  //time = tac();
  //write, "Grid generated in ", time, "sec";

  if (return_param) return parameters ;

  }

//***************************************************

func _create_grid(name,parameters,istart,extension,last_model=) {

  // Creation des repertoires
  system, "mkdir "+name;

  // Ecriture des fichiers de parametres
  n_models = numberof(parameters) ;
  if (is_void(last_model)) last_model = istart + n_models ;
  n_chiffres = floor(log10(last_model) + 1);
  format = "%"+swrite(format="%1i",int(n_chiffres))+"i";

  for (i=1 ; i<=n_models ; i++) {
    number=swrite(format=format,i+istart);
    // Ajout des 0
    for (j=1 ; j<=n_chiffres -1 ; j++) {
      number = streplace(number,strgrep(" ",number),"0");
    }
    file = name+"/"+number+extension;

    write_params, parameters(i),file ;

  }


}

//***************************************************

func _create_grid_gasp(name,parameters,basename,extension) {

  // Creation des repertoires
  system, "mkdir "+name;
  system, "mkdir "+name+"/parameter_files";


  // Ecriture des fichiers de parametres
  dims = dimsof(parameters) ;

  n1 = dims(2) ;
  n2 = dims(3) ;
  n3 = dims(4) ;
  n4 = dims(5) ;
  n5 = dims(6) ;
  n6 = dims(7) ;
  n7 = dims(8) ;

  for (i1=1 ; i1<=n1 ; i1++) {
    for (i2=1 ; i2<=n2 ; i2++) {
      for (i3=1 ; i3<=n3 ; i3++) {
        for (i4=1 ; i4<=n4 ; i4++) {
          for (i5=1 ; i5<=n5 ; i5++) {
            for (i6=1 ; i6<=n6 ; i6++) {
              for (i7=1 ; i7<=n7 ; i7++) {
                number= swrite(format="%1i",i1) + swrite(format="%1i",i2) + swrite(format="%1i",i3) + swrite(format="%1i",i4) + swrite(format="%1i",i5) + swrite(format="%1i",i6) + swrite(format="%1i",i7) ;
                file = name+"/parameter_files/"+basename+"_"+number+extension;
                write_params, parameters(i1,i2,i3,i4,i5,i6,i7), file ;
              }
            }
          }
        }
      }
    }
  }

}


//***************************************************

func write_params(param,file) {

  f=open(file,"w");

  //write, f, param.simu.version, "mcfost version" ;
  write, f, " 3.0       mcfost version" ;
  write, f, "";

  //photon;
  write, f, "#Number of photon packages";
  write, f, param.phot.n2, "nbr_photons_eq_th  : T computation ";
  write, f, param.phot.nlambda, "nbr_photons_lambda : SED computation";
  write, f, param.phot.nimage, "nbr_photons_image : images computation";

  write, f, "";

  //Wavelength
  write, f, "#Wavelength";
  write, f, param.wave.nlamb1, param.wave.lambda_min, param.wave.lambda_max, "n_lambda, lambda_min, lambda_max";
  write, f, param.simu.ltemp , param.simu.lsed ,param.simu.lsed_complete, " compute temperature?, compute sed?, use default wavelength grid ?";
  write, f, param.wave.file,  "wavelength file";
  write, f, param.simu.lsepar,  param.simu.lsepar_pola, "separation of different contributions, stokes parameters";

  write, f, "";

  //Grid
  write, f, "#Grid geometry and size";
  write, f, param.grid.geometry, "1 = cylindrical, 2 = spherical";
  write, f, param.grid.nrad, param.grid.nz, param.grid.n_az, param.grid.nrad_in, "n_rad, nz (or n_theta), n_az, n_rad_in";

  write, f, "";


  //map;
  write, f, "#Maps";
  write, f, param.map.nx, param.map.ny, param.map.map_size, "Nx Ny  map_size";
  //write, f, param.map.nthet, param.map.nphi, param.map.delta, "MC : N_bin_incl, N_bin_az, # of bin where MC is converged ";
  write, f, param.map.RT_imin, param.map.RT_imax, param.map.RT_ntheta, param.map.lRT_i_centered, "RT: imin, imax, n_theta, centered?";
  write, f, param.map.RT_az_min, param.map.RT_az_max, max(param.map.RT_naz,1), "RT: az_min, az_max, n_az angles" ;
  write, f, param.map.dist, "distance (pc)";
  write, f, param.map.PA, "disk PA";

  write, f, "";

  // Scattering method
  write, f, "#Scattering method";
  write, f, param.simu.scatt_method, "0=auto, 1=grain prop, 2=cell prop";
  write, f, param.simu.aniso_method, "1=Mie, 2=hg (2 means loosing pola information)";

  write, f, "";

  // Symetries
  write, f, "#Symetries";
  write, f, param.simu.lsym_ima, "image symetry";
  write, f, param.simu.lsym_centrale, "central symetry";
  write, f, param.simu.lsym_axiale, "axial symmetry (important only if N_phi > 1)" ;

  write, f, "";

  // Dust global properties
  write, f, "#Dust global properties";
  write, f, param.dust.strat_type, param.dust.strat, param.dust.a_strat, "dust settling, exp_strat, a_strat";
  write, f, param.dust.lmigration,  " dust migration?";
  write, f, param.simu.ldust_sublimation,  " sublimate dust?";
  write, f, param.simu.lhydro_eq,  " hydrotstatic equilibrium?";
  write, f, param.simu.lchauff_int, param.simu.viscosity, "viscous heating, viscosity";

  write, f, "";

  //Nbre de zones
  write, f, "#Number of zones : 1 zone = 1 density structure + corresponding grain properties"
  write,f, param.simu.n_zones;

  write, f, "";

  //Density structure
  write, f, "#Density structure";
  for (i=1 ; i <=  param.simu.n_zones ; i++) {
    write, f, param.zones(i).geometry, "zone type : 1 = disk, 2 = tappered-edge disk, 3 = envelope";
    write, f, param.zones(i).mass, param.zones(i).gas_to_dust, "zone dust mass, gas-to-dust mass ratio";
    write, f, param.zones(i).Ho, param.zones(i).Ro, param.zones(i).vert_exp, " scale height, reference radius (AU), unused for envelope, vertical profile exponent (only for debris disk)";
    write, f, param.zones(i).rin,  param.zones(i).edge, param.zones(i).rout, param.zones(i).Rc, "Rin, edge, Rout, Rc (AU) Rc is only used for tappered-edge disks (Rout set to 8*Rc if Rout==0)";
    write, f, param.zones(i).beta, "flaring exponent, unused for envelope";
    write, f, param.zones(i).dens, param.zones(i).gamma_exp, "surface density exponent (or -gamma for tappered-edge disk), usually < 0, -gamma_exp";
  }

  //write, f, "";
  //write, f, "#Cavity";
  //write, f, param.cavity.lcavity ;
  //write, f, param.cavity.H0,  param.cavity.R0 ;
  //write, f, param.cavity.beta ;

  write, f, "";

  //Grain
  write, f, "#Grain properties";
  for (i=1 ; i <=  param.simu.n_zones ; i++) {
    write, f, param.zones(i).n_especes, "Number of species";
    for (j=1 ; j<= param.zones(i).n_especes ; j++) {
      write, f, param.dust_pop(i,j).type, param.dust_pop(i,j).n_components, param.dust_pop(i,j).mixing_rule, param.dust_pop(i,j).porosity, param.dust_pop(i,j).mass_fraction,param.dust_pop(i,j).dhs_maxf, " Grain type (Mie or DHS), N_components, mixing rule (1 = EMT or 2 = coating),  porosity, mass fraction, Vmax (for DHS)"  ;

      for (k=1 ; k<=param.dust_pop(i,j).n_components ; k++) {
        write, f, param.dust_pop(i,j).file(k), param.dust_pop(i,j).component_volume_fraction(k), "optical indices file, volume fraction";
      }
      write, f, param.dust_pop(i,j).type_grain, "heating method : 1 = RE + LTE, 2 = RE + NLTE, 3 = NRE";
      if (param.dust_pop(i,j).n_grains == 1) {
        AMAX = param.dust_pop(i,j).amin ;
      } else {
        AMAX = param.dust_pop(i,j).amax ;
      }
      write, f, param.dust_pop(i,j).amin,  AMAX,  param.dust_pop(i,j).aexp, param.dust_pop(i,j).n_grains, "amin, amax, aexp, nbr_grains";
    }
  }
  write, f, "";

  // mol
  write, f, "#Molecular RT settings";
  write, f, param.mol.lpop(1),  param.mol.lpop_precise(1), param.mol.lLTE(1), param.mol.profile_width(1), "lpop, lprecise_pop, LTE, level_max, profile width";

  write, f, param.mol.vturb(1), "v_turb";
  write, f, param.simu.nmol(1), "nmol";

  for (i=1 ; i<=param.simu.nmol ; i++) {
    write, f, param.mol(i).molecule_file, param.mol(i).level_max, "molecular data filename, level_max" ;
    write, f, param.mol(i).vmax , param.mol(i).n_speed, "vmax (m.s-1), v_turb, n_speed";
    write, f, param.mol(i).cst_abundance, param.mol(i).abundance, param.mol(i).abundance_file, "cst molecule abundance ?, abundance, abundance file";
    write, f, param.mol(i).ray_tracing, param.mol(i).n_lines_ray_tracing, "ray tracing ?,  number of lines in ray-tracing";
    N = param.mol(i).n_lines_ray_tracing ;
    write, f, param.mol(i).transition_numbers(1:N);//, "transition numbers";
  }
  write, f, "";

  // Star
  write, f, "#Star properties"
    write, f, param.simu.n_stars, "Number of stars";
  for (i=1 ; i <= param.simu.n_stars ; i++) {
    write, f, param.stars(i).T, param.stars(i).R, param.stars(i).Mass, param.stars(i).x, param.stars(i).y, param.stars(i).z, param.stars(i).is_bb, "Temp, radius (solar radius),M (solar mass),x,y,z (AU), is a blackbody ?";
    write, f, param.stars(i).spectre ;
    write, f, param.stars(i).fUV, param.stars(i).slope_fUV, "fUV, slope_fUV"
  }

  close, f;
}


//***************************************************

func sedfit(name,dir,Av,obs,lect=,win=,winP=,color=,clear=,erreur=,noX=,progressive_plot=,output_pdf=,rt=,proba_filename=,extension=,correct_dist=,sed_write=,return_chi2=,force_fit=) {
/* DOCUMENT sedfit(name,dir,Av,obs,lect=,win=,winP=,color=,clear=,erreur=,noX=,progressive_plot=,output_pdf=,rt=,proba_filename=,extension=,distance=,sed_write=)
   name must be the same as the one used to create the grid
   dir is the directory with the results of the grid calculation
   SEE ALSO: get_chi2, analyse_bayesienne

   lect if you want to read the models and recompute the chi2
 */

  extern n_params_tot ;

  if (is_void(sed_write)) sed_write=0 ;
  if (is_void(return_chi2)) return_chi2=0;
  if (is_void(force_fit)) force_fit=0 ;
  if (is_void(rt)) rt=0 ;

  if (is_void(extension)) {
    extension="";
  } else {
    extension = "_"+extension ;
  }

  if (is_void(correct_dist)) {
    correct_dist = [1.0] ;
    correct_flux_dist = 1.0 ;
    log_correct_flux_dist = 0.0
    n_D = 1 ;
  } else {
    correct_flux_dist = 1.0/correct_dist^2 ;
    log_correct_flux_dist = log(correct_flux_dist) ;
    n_D = numberof(correct_dist) ;
  }

  f = openb(name+".grid");
  restore, f;
  close, f;


  if (rt) {
    n_incl = mc_grid.ref.map.RT_ntheta ;
  } else {
    n_incl = mc_grid.ref.map.nthet ;
  }
  n_Av = numberof(Av);

  if (is_void(lect)) lect=1 ;
  if (is_void(win)) win=1 ;
  if (is_void(winP)) winP=win+1 ;
  if (is_void(clear)) clear=0 ;
  if (is_void(color)) color="blue" ;
  if (is_void(rt)) rt=0 ;
  if (is_void(proba_filename)) proba_filename=0 ;
  if (is_void(progressive_plot)) progressive_plot=0 ;
  if (is_void(noX)) noX=0 ;
  if (is_void(output_pdf)) output_pdf=0 ;
  if (noX) {
    progressive_plot=0 ;
    output_pdf = 1 ;
  }


  OBS = obs ;
  ou = where(OBS.flux > 0.) ; // donnees
  ou2 = where(OBS.flux == 0.) ;  // upper limit
  if (!is_void(erreur)) {
    OBS(ou).erreur = erreur * OBS(ou).flux ;
  }


  // Fit en log sur les points de donnees
  // en lineaire sur les limites superieures
  OBS(ou).erreur = OBS(ou).erreur / OBS(ou).flux ;
  OBS(ou).flux = log(OBS(ou).flux) ;


  // Dimensions de l'espace des parametres
  D = mc_grid.dims(where(mc_grid.dims > 1)) ;
  n_Parametres = numberof(D) + 1 ; // + 1 pour inclinaison
  if (n_Av > 1) n_Parametres++ ;
  if (n_D > 1) n_Parametres++ ;

  n_models = 1 ; for (i=1 ; i<=numberof(D) ; i++) n_models *= D(i) ;

  N_max = 31998 ;
  n_Parts = n_models / N_max + 1 ;
  Parts = array(string,n_Parts) ;

  // Definition des sous-dossiers
  if (n_Parts > 1) {
    n_chiffres = floor(log10(n_Parts) + 1);
    format = "%"+swrite(format="%1i",int(n_chiffres))+"i";
    for (i=1 ; i<=n_Parts ; i++) {
      Parts(i) = "Part"+swrite(i,format=format)+"/" ;
    }
  }


  // Lecture du premier modele
  n_chiffres = floor(log10(n_models) + 1);
  format = "%"+swrite(format="%1i",int(n_chiffres))+"i";
  number1=swrite(format=format,1);
  for (j=1 ; j<=n_chiffres -1 ; j++) {
    number1 = streplace(number1,strgrep(" ",number1),"0");
  }

  if (n_Parts > 1) {
    dir1 = dir+Parts(1)+number1
  } else {
    dir1 = dir+number1 ;
  }

  mod = open_mcfost(dir1,"th");
  mod_lamb = mod.lamb ;


  // Correction modeles par A_V
  extinction = OpenASCII(get_cwd()+"/SED/extinction_law.dat",prompt=0);

  ext_V = abs(extinction.lamb - 5.47e-01)(mnx);
  kext = extinction.kpa / (1.0-extinction.albedo) ;

  XX =  kext /  kext(ext_V);
  //XX =  extinction.kpa /  extinction.kpa(ext_V);
  XXm = interp( XX , extinction.lamb ,  mod_lamb);

  tau_V=0.4*log(10)*Av;
  correct_Av = exp(-tau_V(,-:1:numberof(mod_lamb)) * XXm(-:1:n_Av,));

  if (numberof(OBS) - n_Parametres <= 0) {
    "Error: too many parameters !!!";
    write, "N obs=", numberof(OBS), "N Param =", n_Parametres ;
    if (!force_fit) {
      return ;
    } else {
      "Forcing fit in 5 sec" ;
      pause, 5000  ;
    }
  }


  n_chiffres = floor(log10(n_models) + 1) ;
  format = "%"+swrite(format="%1i",int(n_chiffres))+"i";

  // Lecture eventuelle
  if (lect) {
    chi2=[];
    chi2=array(double, n_incl, n_models, n_Av, n_D);

    minchi2 = 1.e300 ;

    // Boucle sur les modeles
    for (k=1 ; k<=n_models ; k++) {

      if (!(k%1000)) write, "Model", k, "/", n_models ;

      sk=swrite(format=format,k); for (j=1 ; j<=n_chiffres -1 ; j++) sk = streplace(sk,strgrep(" ",sk),"0"); // Ajout des 0
      ipart = k / N_max + 1 ;


      if (rt) {
        filename = dir+"/"+Parts(ipart)+sk+"/data_th/sed_rt.fits.gz" ;
      } else {
        filename = dir+"/"+Parts(ipart)+sk+"/data_th/sed_mc.fits.gz" ;
      }

      if (IsFile(filename)) {
        sed = cfitsRead(filename) ;
        for (i=1 ; i <= n_incl ; i++) {
          for (j=1 ; j <= n_Av ; j++) {
            //m = exp(interp(log(mod.sed(,i,1,1) * correct_Av(j,)), log(mod.lamb), log(OBS.lamb))) ;
            if (rt) {
              m0 = interp(sed_rt(,i,1,1) * correct_Av(j,), mod_lamb, OBS.lamb) ;
            } else {
              m0 = interp(sed(,i,1,1) * correct_Av(j,), mod_lamb, OBS.lamb) ;
            }
            // m = interp(sed(,i,1,1) * correct_Av(j,), mod_lamb, OBS.lamb) ;

            // Fit en log sur les points de donnees
            m0(ou) = log(max(m0(ou),1e-30)) ;
            m = m0 ;

            for (d=1 ; d<= n_D ; d++) {

              m(ou) = m0(ou) + log(correct_flux_dist(d)) ;
              if (numberof(ou2)) m(ou2) = m0(ou2) * correct_flux_dist(d) ;

              C2 = ( (( m - OBS.flux)/(OBS.erreur))^2)(sum) / (numberof(OBS) - n_Parametres);

              chi2(i,k,j,d) = C2 ;
              if (C2 < minchi2) {
                minchi2 = C2 ;
                //write, "New best model =", filename, "chi2=", C2 ;
                if (progressive_plot) {
                  window, win ; fma ;    lxy, 1, 1 ;
                  plp, obs.flux(ou), obs.lamb(ou), dy=obs.erreur(ou) ;
                  // plot des limites sup : 3 sigma
                  if (numberof(ou2)) plp, 3*obs.erreur(ou2) , obs.lamb(ou2), symbol=7 ;

                  if (clear) {limits ; relimits_log;}

                  order = sort(mod_lamb) ;

                  if (rt) {
                    plg, (sed_rt(,i,1,1) * correct_Av(j,) * correct_flux_dist(d))(order) , mod_lamb(order), color=color, marks=0;
                  } else {
                    plg, (sed(,i,1,1) * correct_Av(j,) * correct_flux_dist(d))(order), mod_lamb(order), color=color, marks=0;
                  }
                }
              }

              if (ieee_test(chi2(i,k,j,d))) chi2(i,k,j,d) = 1.0e9;
            } //d
          } //j

          // TMP
          //if (i==1) write, "Model", k, "/", n_models, "chi2=",  chi2(i,k,:,:)
          //if (i > 1)  chi2(i,k,:,:) = 1e9 ;
          // Fin TMP


        } // i
        //        print, filename, "chi2=", C2, "minchi2=", minchi2 ;


      } else {
        chi2(:,k,:,:) = 1.0e9 ;
      }
    }



    if (rt) {
      f = createb(name+extension+"_rt_chi2.dat");
    } else {
      f = createb(name+extension+"_chi2.dat");
    }
    save, f, chi2;
    close, f;
  } else {
    if (rt) {
      f = openb(name+extension+"_rt_chi2.dat");
    } else {
      f = openb(name+extension+"_chi2.dat");
    }
    restore, f;
    close, f;
  }


  // Analyse Bayesienne
  D2 = [n_Parametres,n_incl] ;
  grow, D2, D ;
  if (n_Av > 1) grow, D2, [n_Av] ;
  if (n_D > 1) grow, D2, [n_D] ;

  D3 = numberof(D) ;
  grow, D3, D ;

  P = array(double,D2) ; // Le nbre de dim est n_Parametres
  liste2 = array(string,D3) ;
  P(*) = exp(-chi2(*)/2.) ;

  if (sum(P) > 0) P = P / sum(P) ;

  // Meilleur modele
  best0 = chi2(*)(mnx);
  best = indexof(chi2,best0);
  best2 = indexof(P,best0);

  //print, "BEST MODEL: ", best , "chi2=", min(chi2);
  write, "Best model: chi2=", min(chi2);

  sbest=swrite(format=format,best(2)); for (j=1 ; j<=n_chiffres -1 ; j++) sbest = streplace(sbest,strgrep(" ",sbest),"0"); // Ajout des 0
  ipart = best(2) / N_max + 1 ;
  filename = dir+"/"+Parts(ipart)+sbest ;

  write, "Model="+filename ; //, "i=", best(1) ; //, "Av=", Av(best(3)), "correct_dist=", correct_dist(best(4)) ;

  /**** Best model values */
  j = 0 ;
  if (n_incl > 1) {
    j++ ;
    x_name = "i" ;
    if (rt) {
      cos_imin= cos(mc_grid.ref.map.RT_imin/180. * pi) ;
      cos_imax= cos(mc_grid.ref.map.RT_imax/180. * pi) ;
      x = acos(span(cos_imin,cos_imax,n_incl+1)(zcen)) * 180./pi ;
    } else {
      x = acos(span(1,0,n_incl+1)(zcen)) * 180./pi ;
    }
    write, x_name, x(best2(j)), best2(j) ;
  }

  for(i=1; i<= numberof(mc_grid.dims); i++) {
    if (mc_grid.dims(i) > 1) {
      j++ ;
      x_name = mc_grid.names(i) ;
      x = mc_grid.values(i,1:mc_grid.dims(i)) ;
      write, x_name, x(best2(j)) ;
    }
  }

  if (numberof(Av) > 1) {
    j++ ;
    x_name = "A_V_" ;
    x = Av ;
    write, x_name, x(best2(j)) ;
  }

  if (n_D > 1) {
    j++ ;
    x_name = "Dist factor" ;
    x = correct_dist ;
    write, x_name, x(best2(j)) ;
  }
  /*****/

  mod = open_mcfost(filename,"th");



  if (noX) {
    window, win, style="landscape.gs", display="", hcp="hcp" ;
  } else {
    window, win, style="landscape.gs" ;
  }
  if (clear) fma;
  plp, obs.flux(ou), obs.lamb(ou), dy=obs.erreur(ou) ;
  if (numberof(ou2)) plp, 3*obs.erreur(ou2), obs.lamb(ou2), symbol=7 ;


  if (clear) {limits ; relimits_log;}
  if (rt) {
    plg, mod.sed_rt(,best(1),1,1) * correct_Av(best(3),) * correct_flux_dist(best(4)), mod_lamb, color=color;
    if (sed_write)  {
      write, "# lambda (micron) lambda.F_lambda (W.m^-2)" ;
      write, mod.lamb,  mod.sed_rt(,best(1),1,1) * correct_Av(best(3),) * correct_flux_dist(best(4)) ;
    }
  } else {
    plg, mod.sed(,best(1),1,1) * correct_Av(best(3),) * correct_flux_dist(best(4)), mod_lamb, color=color;
    if (sed_write)  {
      write, "# lambda (micron) lambda.F_lambda (W.m^-2)" ;
      write, mod.lamb,  mod.sed(,best(1),1,1) * correct_Av(best(3),) * correct_flux_dist(best(4)) ;
    }
  }

  lxy, 1, 1 ;
  if (output_pdf) pdf, name+extension+"_SED.pdf" ;

  P = _analyse_bayesienne(name,mc_grid,chi2, Av=Av,win=winP,noX=noX,rt=rt,clear=clear,proba_filename=proba_filename,output_pdf=output_pdf,extension=extension) ;

  // most probable
  //liste2(*) = liste(*) ;
  //liste2(m) ;
  //write, "most probable", m ;

  if (return_chi2) {
    C2 = array(double,D2) ; // Le nbre de dim est n_Parametres
    C2(*) = chi2(*) ;
    return C2 ;
  } else {
    return P ;
  }
}

//***************************************************

func _analyse_bayesienne(name, mc_grid,chi2, Av=,win=,noX=,rt=,clear=,proba_filename=,output_pdf=,extension=) {

  // Dimensions de l'espace des parametres
  D = mc_grid.dims(where(mc_grid.dims > 1)) ;
  n_Parametres = numberof(D) + 1 ; // + 1 pour inclinaison

  n_Av = numberof(Av) ; n_D = 0 ;
  if (n_Av > 1) n_Parametres++ ;
  if (n_D > 1) n_Parametres++ ;
  n_models = 1 ; for (i=1 ; i<=numberof(D) ; i++) n_models *= D(i) ;

  if (rt) {
    n_incl = mc_grid.ref.map.RT_ntheta ;
  } else {
    n_incl = mc_grid.ref.map.nthet ;
  }

  D2 = [n_Parametres,n_incl] ;
  grow, D2, D ;
  if (n_Av > 1) grow, D2, [n_Av] ;
  if (n_D > 1) grow, D2, [n_D] ;

  // Tableau de proba
  P = array(double,D2) ; // Le nbre de dim est n_Parametres
  liste2 = array(string,D3) ;
  P(*) = exp(-chi2(*)/2.) ;

  if (sum(P) > 0) P = P / sum(P) ;


  if (noX) {
    window, win, display="", hcp="hcp" ;
  } else {
    window, win;
  }
  if (clear) fma;
  gs_nm, win, 3, 4, dx=0.075, dy=0.075, square=0 ;

 j = 0 ;
 if (n_incl > 1) {
   plsys, j+1 ;
   x_name = "i" ;
    if (rt) {
      cos_imin= cos(mc_grid.ref.map.RT_imin/180. * pi) ;
      cos_imax= cos(mc_grid.ref.map.RT_imax/180. * pi) ;
      x = acos(span(cos_imin,cos_imax,n_incl+1)(zcen)) * 180./pi ;
    } else {
      x = acos(span(1,0,n_incl+1)(zcen)) * 180./pi ;
    }
    p = P(,*)(,sum) ;
    plh, p, x, color=color ;
    if (proba_filename) {
      write, "writing ", proba_filename ;
      f = open(proba_filename,"w") ;
      write, f, "# "+x_name;
      write, f, x, p ;
    }

    xytitles, x_name, , [0.02,0.02];
    limits ;
    limits, 0, 90, 0., 1.1*limits()(4) ;

    gridxy, degrees=1;
  }

  m = array(int,numberof(D)) ;

  for(i=1; i<= numberof(mc_grid.dims); i++) {
    if (mc_grid.dims(i) > 1) {
      j++ ;
      plsys, j+1 ;
      x_name = mc_grid.names(i) ;
      x = mc_grid.values(i,1:mc_grid.dims(i)) ;

      X = x(sort(x)) ;
      if ( abs(X(0)-X(-1)) > 2*abs(X(2)- X(1)) ) {
        lxy, 1, 0 ;
      } else {
        lxy, 0, 0 ;
      }

      if (j==1) {
        p = P(sum,,*)(,sum) ;
        m(1) = p(mxx) ;
      } else if (j==2) {
        p = P(sum,sum,,*)(,sum) ;
        m(2) = p(mxx) ;
      } else if (j==3) {
        p = P(sum,sum,sum,,*)(,sum) ;
        m(3) = p(mxx) ;
      } else if (j==4) {
        p = P(sum,sum,sum,sum,,*)(,sum) ;
        m(4) = p(mxx) ;
      } else if (j==5) {
        p = P(sum,sum,sum,sum,sum,,*)(,sum) ;
        m(5) = p(mxx) ;
      } else if (j==6) {
        p = P(sum,sum,sum,sum,sum,sum,,*)(,sum) ;
        m(6) = p(mxx) ;
      } else if (j==7) {
        p = P(sum,sum,sum,sum,sum,sum,sum,,*)(,sum) ;
        m(7) = p(mxx) ;
      } else if (j==8) {
        p = P(sum,sum,sum,sum,sum,sum,sum,sum,,*)(,sum) ;
        m(8) = p(mxx) ;
      } else if (j==9) {
        p = P(sum,sum,sum,sum,sum,sum,sum,sum,sum,,*)(,sum) ;
        m(9) = p(mxx) ;
      } else if (j==10) {
        p = P(sum,sum,sum,sum,sum,sum,sum,sum,sum,sum,,*)(,sum) ;
        m(10) = p(mxx) ;
      } else if (j==11) {
        p = P(sum,sum,sum,sum,sum,sum,sum,sum,sum,sum,sum,,*)(,sum) ;
        m(11) = p(mxx) ;
      } else if (j==12) {
        p = P(sum,sum,sum,sum,sum,sum,sum,sum,sum,sum,sum,sum,,*)(,sum) ;
        m(12) = p(mxx) ;
      } else if (j > 12) {
        error, "Too many dims " ;
      }

      plh, p, x, color=color ;

      xytitles, x_name, , [0.02,0.02];
      limits ;
      limits, , , 0., 1.1*limits()(4) ;
      if (proba_filename) {
        write,f, " " ;
        write, f, "# "+x_name;
        write, f, x, p ;
      }

    }
  }


  if (numberof(Av) > 1) {
    plsys, j+2;
    x_name = "A_V_" ;
    x = Av ;
    if (n_D > 1) {
      p = P(..,sum) ;
      p = p(*,)(sum,) ;
    } else {
      p = P(*,)(sum,) ;
    }
    plh, p, x, color=color ;
    xytitles, x_name, , [0.02,0.02];



    limits ;
    limits, , , 0., 1.1*limits()(4) ;

    if (proba_filename) {
      write, f, " " ;
      write, f, "# "+x_name;
      write, f, x, p ;
    }
  }

  if (n_D > 1) {
    plsys, j+3;
    x_name = "Dist factor" ;
    x = correct_dist ;
    p = P(*,)(sum,) ;
    plh, p, x, color=color ;
    xytitles, x_name, , [0.02,0.02];

    limits ;
    limits, , , 0., 1.1*limits()(4) ;

    if (proba_filename) {
      write,f, " " ;
      write, f, "# "+x_name;
      write, f,x, p ;
    }
  }

  if (proba_filename) close, f ;
  if (output_pdf) pdf, name+extension+"_Proba.pdf" ;


  //*************

  best0 = chi2(*)(mnx);
  best = indexof(chi2,best0);
  best2 = indexof(P,best0);

  //print, "BEST MODEL: ", best , "chi2=", min(chi2);
  write, "Best model: chi2=", min(chi2);

  sbest=swrite(format=format,best(2)); for (j=1 ; j<=n_chiffres -1 ; j++) sbest = streplace(sbest,strgrep(" ",sbest),"0"); // Ajout des 0
  ipart = best(2) / N_max + 1 ;
  filename = dir+"/"+Parts(ipart)+sbest ;

  write, "Model="+filename ; //, "i=", best(1) ; //, "Av=", Av(best(3)), "correct_dist=", correct_dist(best(4)) ;

  /**** Best model values */
  j = 0 ;
  if (n_incl > 1) {
    j++ ;
    x_name = "i" ;
    if (rt) {
      cos_imin= cos(mc_grid.ref.map.RT_imin/180. * pi) ;
      cos_imax= cos(mc_grid.ref.map.RT_imax/180. * pi) ;
      x = acos(span(cos_imin,cos_imax,n_incl+1)(zcen)) * 180./pi ;
    } else {
      x = acos(span(1,0,n_incl+1)(zcen)) * 180./pi ;
    }
    write, x_name, x(best2(j)), best2(j) ;
  }

  for(i=1; i<= numberof(mc_grid.dims); i++) {
    if (mc_grid.dims(i) > 1) {
      j++ ;
      x_name = mc_grid.names(i) ;
      x = mc_grid.values(i,1:mc_grid.dims(i)) ;
      write, x_name, x(best2(j)) ;
    }
  }

  if (numberof(Av) > 1) {
    j++ ;
    x_name = "A_V_" ;
    x = Av ;
    write, x_name, x(best2(j)) ;
  }

  if (n_D > 1) {
    j++ ;
    x_name = "Dist factor" ;
    x = correct_dist ;
    write, x_name, x(best2(j)) ;
  }

  //****

  return P ;

}

//***************************************************

func sedfit2(name,dir,Av,obs,lect=,win=,color=,clear=,erreur=,noX=,progressive_plot=,output_pdf=,rt=,output_proba=,extension=,correct_dist=,dist_source=,sed_write=,return_chi2=,force_fit=) {
/* DOCUMENT sedfit2(name,dir,Av,obs,lect=,win=color=,clear=,erreur=,noX=,progressive_plot=,output_pdf=,rt=,output_proba=,extension=,distance=,sed_write=)
   name must be the same as the one used to create the grid
   dir is the directory with the results of the grid calculation

   Version multi-source, pas sur qu'il faille garder
   SEE ALSO: get_chi2, analyse_bayesienne
 */

  extern n_params_tot ;

  if (is_void(sed_write)) sed_write=0 ;
  if (is_void(return_chi2)) return_chi2=0;
  if (is_void(force_fit)) force_fit=0 ;

  if (is_void(extension)) {
    extension="";
  } else {
    extension = "_"+extension ;
  }

  if (is_void(correct_dist)) {
    correct_dist = [1.0] ;
    correct_flux_dist = 1.0 ;
    log_correct_flux_dist = 0.0
    n_D = 1 ;
  } else {
    correct_flux_dist = 1.0/correct_dist^2 ;
    log_correct_flux_dist = log(correct_flux_dist) ;
    n_D = numberof(correct_dist) ;
  }


  f = openb(name+".grid");
  restore, f;
  close, f;


  if (rt) {
    n_incl = mc_grid.ref.map.RT_ntheta ;
  } else {
    n_incl = mc_grid.ref.map.nthet ;
  }
  n_Av = numberof(Av);

  if (is_void(lect)) lect=1 ;
  if (is_void(win)) win=1 ;
  if (is_void(clear)) clear=0 ;
  if (is_void(color)) color="blue" ;
  if (is_void(rt)) rt=0 ;
  if (is_void(output_proba)) output_proba=0 ;
  if (is_void(progressive_plot)) progressive_plot=0 ;
  if (is_void(noX)) noX=0 ;
  if (is_void(output_pdf)) output_pdf=0 ;
  if (noX) {
    progressive_plot=0 ;
    output_pdf = 1 ;
  }


  // obs est une structure multi-source
  // TODO : ajouter un flag pour la presence des lambda

  OBS = obs ;

  if (!is_void(erreur)) {
    OBS(ou).erreur = erreur * OBS(ou).flux ;
  }


  // Dimensions de l'espace des parametres
  D = mc_grid.dims(where(mc_grid.dims > 1)) ;
  n_Parametres = numberof(D) + 1 ; // + 1 pour inclinaison
  n_models = 1 ; for (i=1 ; i<=numberof(D) ; i++) n_models *= D(i) ;


  N_max = 31998 ;
  n_Parts = n_models / N_max + 1 ;
  Parts = array(string,n_Parts) ;

  // Definition des sous-dossiers
  if (n_Parts > 1) {
    n_chiffres = floor(log10(n_Parts) + 1);
    format = "%"+swrite(format="%1i",int(n_chiffres))+"i";
    for (i=1 ; i<=n_Parts ; i++) {
      Parts(i) = "Part"+swrite(i,format=format)+"/" ;
    }
  }


  // Lecture du premier modele
  n_chiffres = floor(log10(n_models) + 1);
  format = "%"+swrite(format="%1i",int(n_chiffres))+"i";
  number1=swrite(format=format,1);
  for (j=1 ; j<=n_chiffres -1 ; j++) {
    number1 = streplace(number1,strgrep(" ",number1),"0");
  }

  if (n_Parts > 1) {
    dir1 = dir+Parts(1)+number1
  } else {
    dir1 = dir+number1 ;
  }

  mod = open_mcfost(dir1,"th");
  mod_lamb = mod.lamb ;


  if (is_void(dist_source)) {
    "ERROR : dist_source needed"
      error ;
  } else {
    CorrectFluxSource = (mod.P.map.dist/dist_source)^2 ;
    log_CorrectFluxSource = log(CorrectFluxSource) ;
  }


  // Correction modeles par A_V
  extinction = OpenASCII(get_cwd()+"/SED/extinction_law.dat",prompt=0);

  ext_V = abs(extinction.lamb - 5.47e-01)(mnx);
  kext = extinction.kpa / (1.0-extinction.albedo) ;

  XX =  kext /  kext(ext_V);
  //XX =  extinction.kpa /  extinction.kpa(ext_V);
  XXm = interp( XX , extinction.lamb ,  mod_lamb);

  tau_V=0.4*log(10)*Av;
  correct_Av = exp(-tau_V(,-:1:numberof(mod_lamb)) * XXm(-:1:n_Av,));

  if (n_Av > 1) n_Parametres++ ;
  if (n_D > 1) n_Parametres++ ;

  Dims_obs = dimsof(OBS) ;
  if (Dims_obs(1) == 2) {
    // sources multiples
    n_sources = Dims_obs(2) ;
    obs_lamb = OBS(1,).lamb ;
  } else {
    // 1 seule source : on ajoute 1 dim
    n_sources = 1 ;
    obs_lamb = OBS.lamb ;
    OBS = OBS(-:1:1,) ;
  }

  // TODO : faire le test pour chaque source
  for (s=1; s<=n_sources ; s++) {
    if (numberof(where(OBS(s,).flag)) - n_Parametres <= 0) {
      write, "Error: too many parameters for source",s,"!!!";
      write, "N obs=",numberof(where(OBS(s,).flag)) , "N Param =", n_Parametres ;
      if (!force_fit) {
        return ;
      } else {
        "Forcing fit in 5 sec" ;
        pause, 5000  ;
      }
    }
  }

  // Lecture eventuelle
  if (lect) {
    chi2=[];
    write, "Size of chi2 array = ", n_sources * n_incl * n_models * n_Av * n_D * 4/1024.^3, "GB" ;

    chi2=array(float, n_sources, n_incl, n_models, n_Av, n_D);


    minchi2 = 1.e300 ;
    //for (k=1 ; k<=n_models ; k++) {
    for (k=1 ; k<=10000 ; k++) {
      // TODO : nom du modele
      //Parts + string avec les 0

      n_chiffres = floor(log10(n_models) + 1);
      format = "%"+swrite(format="%1i",int(n_chiffres))+"i";
      number=swrite(format=format,k);
      // Ajout des 0
      for (j=1 ; j<=n_chiffres -1 ; j++) number = streplace(number,strgrep(" ",number),"0");
      ipart = k / N_max + 1 ;


      if (rt) {
        filename = dir+"/"+Parts(ipart)+number+"/data_th/sed_rt.fits.gz" ;
      } else {
        filename = dir+"/"+Parts(ipart)+number+"/data_th/sed_mc.fits.gz" ;
      }
      write, filename ;

      if (IsFile(filename)) {
        k ;
        sed = cfitsRead(filename) ;
        for (i=1 ; i <= n_incl ; i++) {
          for (j=1 ; j <= n_Av ; j++) {
            m0 = interp(sed(,i,1,1) * correct_Av(j,), mod_lamb, obs_lamb) ;
            // Fit en log sur les points de donnees
            ml0 = log(max(m0,1e-30)) ;

            m = m0 ; ml = ml0 ;

            for (d=1 ; d<= n_D ; d++) {

              // Les corrections possibles sur les sources
              ml1 = ml0 + log_correct_flux_dist(d) ;
              m1 = m0 * correct_flux_dist(d) ;

              for (s=1 ; s<=n_sources ; s++) {
                ou_data = where(OBS(s,).flag & !OBS(s,).upper) ;
                ou_upper = where(OBS(s,).flag & OBS(s,).upper) ;

                // Correction specifique a la source
                ml = ml1 + log_CorrectFluxSource(s) ;
                m = m1 * CorrectFluxSource(s) ;

                //C2 = ( (( m - OBS.flux)/(OBS.erreur))^2)(sum) / (numberof(OBS) - n_Parametres);
                C1 = 0. ; C2 = 0.
                if (numberof(ou_data))  C1 = ( (( ml(ou_data) - OBS(s,ou_data).flux)/(OBS(s,ou_data).erreur))^2 )(sum) ;
                if (numberof(ou_upper)) C2 = ( (( m(ou_upper) - OBS(s,ou_upper).flux)/(OBS(s,ou_upper).erreur))^2 )(sum) ;
                C2 = (C1+C2) / (numberof(OBS) - n_Parametres);

                write, s, C2 ;
                chi2(s,i,k,j,d) = C2 ;

                if (ieee_test(chi2(s,i,k,j,d))) chi2(s,i,k,j,d) = 1.0e9;
              } // s
            } //d
          } //j
        } // i
      } else {   // fichier manquant
        chi2(s,:,k,:,:) = 1.0e9 ;
      }
    }

    /*
    if (rt) {
      f = createb(name+extension+"_rt_chi2.dat");
    } else {
      f = createb(name+extension+"_chi2.dat");
    }
    save, f, chi2;
    close, f;
     */
  } else {
    if (rt) {
      f = openb(name+extension+"_rt_chi2.dat");
    } else {
      f = openb(name+extension+"_chi2.dat");
    }
    restore, f;
    close, f;
  }

  }

//***************************************************

func read_obs(liste_sources,no_upper=) {
  /* DOCUMENT
     Lit une liste de fichiers d'obs et les placent dans un structure commune
     flag=1 indique que la longueur d'onde est presente pour la source
   SEE ALSO:
 */

  n_sources = numberof(liste_sources) ;
  if (is_void(no_upper)) no_upper = 0 ;

  obs = OpenASCII(liste_sources(1),prompt=0) ;
  lambda = obs.lamb ;

  // Detection de toutes les longeurs d'onde
  for (i=2 ; i<=n_sources ; i++) {
    liste_sources(i)
    obs = OpenASCII(liste_sources(i),prompt=0) ;

    for (l=1 ; l<= numberof(obs) ; l++) {
      ou = where(abs(lambda -obs.lamb(l)) < 1e-6) ;
      if (!numberof(ou)) grow, lambda, obs.lamb(l) ;
    }
  }
  lambda = lambda(sort(lambda)) ;

  write, numberof(lambda), "wavelengths found" ;

  // Creation de la structure
  OBS = array(obs_struct,n_sources,numberof(lambda)) ;
  OBS(,).lamb = lambda(-,) ;

  OBS.flag = 0 ;
  OBS.upper = 0 ;
  for (i=1 ; i<=n_sources ; i++) {
    obs = OpenASCII(liste_sources(i),prompt=0) ;

    for (l=1 ; l<= numberof(obs) ; l++) {
      ou = where(abs(lambda -obs.lamb(l)) < 1e-6) ;
      OBS(i,ou).flux = obs(l).flux ;
      OBS(i,ou).erreur = obs(l).erreur ;
      if (!no_upper) OBS(i,ou).upper = obs(l).upper ;
      OBS(i,ou).flag = 1 ;
    }
  }

  return OBS ;
}

//***************************************************

struct obs_struct {
  float lamb, flux, erreur ;
  int upper, flag ;
}

//***************************************************

func get_chi2(name,extension=,rt=) {

/* DOCUMENT get_chi2(name,extension=,rt=)

   SEE ALSO: sedfit
 */
  if (is_void(rt)) rt=0 ;

  if (is_void(extension)) {
    extension="";
  } else {
    extension = "_"+extension ;
  }

  if (rt) {
    f = openb(name+extension+"_rt_chi2.dat");
  } else {
    f = openb(name+extension+"_chi2.dat");
  }
  restore, f ;
  close, f ;

  D2 = [n_Parametres,n_incl] ;
  grow, D2, D ;
  if (n_Av > 1) grow, D2, [n_Av] ;
  if (n_D > 1) grow, D2, [n_D] ;

  D3 = numberof(D) ;
  grow, D3, D ;

  C2 = array(double,D2) ; // Le nbre de dim est n_Parametres
  C2(*) = chi2(*) ;

  return chi2 ;
}

//***************************************************

func analyse_bayesienne(name,dir,Av,P,win=,color=,clear=,rt=,extension=,correct_dist=,sed_write=,proba_filename=,output_pdf=,noX=) {

  if (is_void(sed_write)) sed_write=0 ;
  if (is_void(proba_filename)) proba_filename=0;
  if (is_void(lect)) lect=1 ;
  if (is_void(win)) win=1 ;
  if (is_void(clear)) clear=0 ;
  if (is_void(color)) color="blue" ;
  if (is_void(rt)) rt=0 ;
  if (is_void(extension)) {
    extension="";
  } else {
    extension = "_"+extension ;
  }
  if (is_void(output_pdf)) output_pdf=0 ;
  if (is_void(noX)) noX = 0 ;
  if (noX) output_pdf = 1 ;


  if (is_void(correct_dist)) {
    correct_dist = [1.0] ;
    correct_flux_dist = 1.0 ;
    log_correct_flux_dist = 0.0
    n_D = 1 ;
  } else {
    correct_flux_dist = 1.0/correct_dist^2 ;
    log_correct_flux_dist = log(correct_flux_dist) ;
    n_D = numberof(correct_dist) ;
  }

  f = openb(name+".grid");
  restore, f;
  close, f;


  // Dimensions de l'espace des parametres
  D = mc_grid.dims(where(mc_grid.dims > 1)) ;
  n_Parametres = numberof(D) + 1 ; // + 1 pour inclinaison
  n_models = 1 ; for (i=1 ; i<=numberof(D) ; i++) n_models *= D(i) ;


  N_max = 31998 ;
  n_Parts = n_models / N_max + 1 ;
  Parts = array(string,n_Parts) ;

  // Definition des sous-dossiers
  if (n_Parts > 1) {
    n_chiffres = floor(log10(n_Parts) + 1);
    format = "%"+swrite(format="%1i",int(n_chiffres))+"i";
    for (i=1 ; i<=n_Parts ; i++) {
      Parts(i) = "Part"+swrite(i,format=format)+"/" ;
    }
  }

  // Meilleur modele
  if (numberof(Av) > 1) {
    if (numberof(correct_dist) > 1) {
      A = P(,*,,) ;
    } else {
      A = P(,*,) ;
    }
  } else {
    if (numberof(correct_dist) > 1) {
      A = P(,*,) ;
    } else {
      A = P(,*) ;
    }
  }


  best = A(*)(mxx); best;
  best = indexof(A,best);

  //liste = dir+lsdir(dir);
  //liste = liste(sort(liste));
  //  n_models = numberof(liste);
  if (rt) {
    n_incl = mc_grid.ref.map.RT_ntheta ;
  } else {
    n_incl = mc_grid.ref.map.nthet ;
  }
  n_Av = numberof(Av) ;


  //print, "BEST MODEL: ", best ;//, "chi2=", min(chi2);



  // Lecture du premier


  n_chiffres = floor(log10(n_models) + 1);

  format = "%"+swrite(format="%1i",int(n_chiffres))+"i";
  sbest=swrite(format=format,best(2)); for (j=1 ; j<=n_chiffres -1 ; j++) sbest = streplace(sbest,strgrep(" ",sbest),"0"); // Ajout des 0
  ipart = best(2) / N_max + 1 ;
  filename = dir+"/"+Parts(ipart)+sbest ;

  if (numberof(best) > 2) {
    write, "Model="+filename, "i=", best(1), "Av=", Av(best(3)) ; //, "correct_dist=", correct_dist(best(4)) ;
  } else {
    write, "Model="+filename, "i=", best(1) ;
  }
  mod = open_mcfost(filename,"th");
  mod_lamb = mod.lamb ;

  // Correction modeles par A_V
  extinction = OpenASCII(get_cwd()+"/SED/extinction_law.dat",prompt=0);

  ext_V = abs(extinction.lamb - 5.47e-01)(mnx);
  XX =  extinction.kpa /  extinction.kpa(ext_V);
  XXm = interp( XX , extinction.lamb ,  mod_lamb);

  tau_V=0.4*log(10)*Av;
  correct_Av = exp(-tau_V(,-:1:numberof(mod_lamb)) * XXm(-:1:n_Av,));

  if (noX) {
    window, win, style="landscape.gs", display="", hcp="hcp" ;
  } else {
    window, win, style="landscape.gs" ;
  }
  if (clear) fma;
  if (numberof(ou2)) plp, 3*obs.erreur(ou2), obs.lamb(ou2), symbol=7 ;

  if (clear) {limits ; } //relimits_log;}

  /*
  if (rt == 1) {
    plg, mod.sed_rt(,best(1),1,1) * correct_Av(best(3),) * correct_flux_dist(best(4)), mod_lamb, color=color;
    if (sed_write)  {
      write, "# lambda (micron) lambda.F_lambda (W.m^-2)" ;
      if (numberof(correct_dist) > 1) {
        //write, mod.lamb,  mod.sed_rt(,best(1),1,1) * correct_Av(best(3),) * correct_flux_dist(best(4)) ;
      } else {
        //write, mod.lamb,  mod.sed_rt(,best(1),1,1) * correct_Av(best(3),) ;
      }
    }
  } else {
    plg, mod.sed(,best(1),1,1) * correct_Av(best(3),) * correct_flux_dist(best(4)), mod_lamb, color=color;
    if (sed_write)  {
      write, "# lambda (micron) lambda.F_lambda (W.m^-2)" ;
      write, mod.lamb,  mod.sed(,best(1),1,1) * correct_Av(best(3),) * correct_flux_dist(best(4)) ;
    }
  }
  */

  lxy, 1, 1 ;

  if (clear) {
    if (noX) {
      window, win+1, style="landscape.gs", display="", hcp="hcp" ;
    } else {
      window, win+1, style="landscape.gs" ;
    }
    gs_nm, win+1, 3, 4, dx=0.075, dy=0.075, square=0 ;
  } else {
    window, win+1 ;
  }

  plsys, 1 ;
  x_name = "i" ;
  if (rt) {
    cos_imin= cos(mc_grid.ref.map.RT_imin/180. * pi) ;
    cos_imax= cos(mc_grid.ref.map.RT_imax/180. * pi) ;
    x = acos(span(cos_imin,cos_imax,n_incl+1)(zcen)) * 180./pi ;
  } else {
    x = acos(span(1,0,n_incl+1)(zcen)) * 180./pi ;
  }
  p = P(,*)(,sum) ;
  plh, p, x, color=color ;
  if (proba_filename) {
    write, "writing "+proba_filename ;
    f = open(proba_filename,"w") ;
    write, f, "# "+x_name;
    write, f, x, p ;
  }

  xytitles, x_name, , [0.02,0.02];
  limits ;
  limits, 0, 90, 0., 1.1*limits()(4) ;


  gridxy, degrees=1;

  D = mc_grid.dims(where(mc_grid.dims > 1)) ;
  m = array(int,numberof(D)) ;


  j = 0 ;
  for(i=1; i<= min(n_params_tot,numberof(mc_grid.dims)) ; i++) {
    if (mc_grid.dims(i) > 1) {
      j++ ;
      plsys, j+1 ;
      x_name = mc_grid.names(i) ;
      x = mc_grid.values(i,1:mc_grid.dims(i)) ;
      if (j==1) {
        p = P(sum,,*)(,sum) ;
        m(1) = p(mxx) ;
      } else if (j==2) {
        p = P(sum,sum,,*)(,sum) ;
        m(2) = p(mxx) ;
      } else if (j==3) {
        p = P(sum,sum,sum,,*)(,sum) ;
        m(3) = p(mxx) ;
      } else if (j==4) {
        p = P(sum,sum,sum,sum,,*)(,sum) ;
        m(4) = p(mxx) ;
      } else if (j==5) {
        p = P(sum,sum,sum,sum,sum,,*)(,sum) ;
        m(5) = p(mxx) ;
      } else if (j==6) {
        p = P(sum,sum,sum,sum,sum,sum,,*)(,sum) ;
        m(6) = p(mxx) ;
      } else if (j==7) {
        p = P(sum,sum,sum,sum,sum,sum,sum,,*)(,sum) ;
        m(7) = p(mxx) ;
      } else if (j==8) {
        p = P(sum,sum,sum,sum,sum,sum,sum,sum,,*)(,sum) ;
        m(8) = p(mxx) ;
      }

      plh, p, x, color=color ;
      xytitles, x_name, , [0.02,0.02];
      limits ;
      limits, , , 0., 1.1*limits()(4) ;

      X = x(sort(x)) ;
      if ( abs(X(0)-X(-1)) > 2*abs(X(2)- X(1)) ) {
        lxy, 1, 0 ;
      } else {
        lxy, 0, 0 ;
      }

      if (proba_filename) {
        write, f, " " ;
        write, f, "# "+x_name;
        write, f, x, p ;
      }

    }
  }


  if (numberof(Av) > 1) {
    plsys, j+2;
    x_name = "A_V_" ;
    x = Av ;
    if (n_D > 1) {
      p = P(..,sum) ;
      p = p(*,)(sum,) ;
    } else {
      p = P(*,)(sum,) ;
    }
    plh, p, x, color=color ;
    xytitles, x_name, , [0.02,0.02];



    limits ;
    limits, , , 0., 1.1*limits()(4) ;

    if (proba_filename) {
      write, f, " " ;
      write, f, "# "+x_name;
      write, f, x, p ;
    }
  }

  if (n_D > 1) {
    plsys, j+3;
    x_name = "Dist factor" ;
    x = correct_dist ;
    p = P(*,)(sum,) ;
    plh, p, x, color=color ;
    xytitles, x_name, , [0.02,0.02];

    limits ;
    limits, , , 0., 1.1*limits()(4) ;

    if (proba_filename) {
      write, f, " " ;
      write, f, "# "+x_name;
      write, f, x, p ;
    }
  }

  if (proba_filename) close, f ;
  if (output_pdf) pdf, name+extension+"_Proba_tot.pdf" ;

  /**** Best model values */
  j = 0 ;
  if (n_incl > 1) {
    j++ ;
    x_name = "i" ;
    if (rt) {
      cos_imin= cos(mc_grid.ref.map.RT_imin/180. * pi) ;
      cos_imax= cos(mc_grid.ref.map.RT_imax/180. * pi) ;
      x = acos(span(cos_imin,cos_imax,n_incl+1)(zcen)) * 180./pi ;
    } else {
      x = acos(span(1,0,n_incl+1)(zcen)) * 180./pi ;
    }
    write, x_name, x(best2(j)), best2(j) ;
  }

  for(i=1; i<= numberof(mc_grid.dims); i++) {
    if (mc_grid.dims(i) > 1) {
      j++ ;
      x_name = mc_grid.names(i) ;
      x = mc_grid.values(i,1:mc_grid.dims(i)) ;
      write, x_name, x(best2(j)) ;
    }
  }

  if (numberof(Av) > 1) {
    j++ ;
    x_name = "A_V_" ;
    x = Av ;
    write, x_name, x(best2(j)) ;
  }

  if (n_D > 1) {
    j++ ;
    x_name = "Dist factor" ;
    x = correct_dist ;
    write, x_name, x(best2(j)) ;
  }
  /*****/

}

//***************************************************

func add_grids(grid1,grid2) {


}

//***************************************************

func start_grid(grid_name,cluster=,omp_num_threads=,n_nodes=,walltime=,root_dir=,options=) {
/* DOCUMENT start_grid(grid_name,cluster=,omp_num_threads=,n_nodes=,walltime=,root_dir=,options=)

   SEE ALSO:
 */

  write, "Starting grid "+grid_name+"..." ;

  if (is_void(cluster)) cluster="fostino" ;
  if (is_void(omp_num_threads)) omp_num_threads=8 ;
  if (is_void(options)) options = " " ;
  if (is_void(root_dir)) root_dir = "." ;


  if (is_void(n_nodes)) {
    write, "n_nodes is required in start_grid" ;
    return ;
  }

  if (is_void(walltime)) walltime=240 ;

  write_OAR_script, omp_num_threads, options=options ;

  dir = root_dir+"/"+grid_name ;

  // To be updated depending on the cluster
  system, "ssh "+cluster+" 'mkdir "+dir+"'" ;
  system, "scp "+grid_name+"_P.tgz distrib.sh "+cluster+":"+dir ;
  write, "Untaring files and submitting ..." ;

  system, "ssh "+cluster+" 'cd "+dir+" && chmod +x distrib.sh && tar xzf "+grid_name+"_P.tgz && "+
    "/bin/tcsh -c \"repeat "+swrite(n_nodes,format="%i")+" oarsub -n "+grid_name+" -lnodes=1,walltime="+swrite(walltime,format="%i")+" ./distrib.sh\"' " ;
  write, "Done" ;


  write, "Grid started" ;
 }

//***************************************************

func write_OAR_script(omp_num_threads,options=,post_process=,dir=) {

  if (is_void(omp_num_threads)) omp_num_threads = 8 ;
  if (is_void(options)) options = " " ;
  if (is_void(dir)) dir = "." ;

  n_options = numberof(options) ;

  f = open(dir+"/distrib.sh","w") ;
  write, f, "#!/bin/bash" ;
  write, f, "export mcfost=/home/ciment/cpinte/mcfost_cigri/mcfost" ;   // ${HOME} --> cpinte
  write, f, "export MCFOST_UTILS=/home/ciment/cpinte/mcfost/utils/" ;
  write, f, "export MCFOST_NO_DISCLAIMER=1" ;
  write, f, "export OMP_NUM_THREADS="+swrite(omp_num_threads,format="%i") ;
  write, f, " " ;
  for (i=1 ; i<=n_options ; i++) {
    if (options(i) == "-prodimo") {
      write, f, "export prodimo=/home/ciment/cpinte/ProDiMo/src_develop/prodimo" ;
      write, f, "export ProDiMo_datapath=/home/ciment/cpinte/ProDiMo/data" ;
      write, f, " " ;
    }
  }
  write, f, "export node_dir=\"job_$OAR_JOBID\"" ;
  write, f, "mkdir $node_dir" ;
  write, f, "cp *.lambda ../*.lambda $node_dir" ;
  write, f, " ";
  write, f, "for param in *.par;do" ;
  write, f, "if [ -f $param ]; then" ;
  write, f, "      mv \"$param\" \"$node_dir\"" ;
  write, f, "      cd \"$node_dir\"" ;
  write, f, "      new_dir=\"`ls \"$param\" | sed s/.par//`\"" ;
  write, f, "      mkdir \"$new_dir\"" ;
  write, f, "      rm -rf data_*";
  for (i=1 ; i<=n_options ; i++) {
    if (i == 1) {
      write, f, "      $mcfost $param "+options(i)+" > mcfost.log" ;
    } else {
      write, f, "      $mcfost $param "+options(i)+" >> mcfost.log" ;
    }
    if (options(i) == "-prodimo") {
      //write, f, "      mkdir data_ProDiMo ; cd data_ProDiMo ; ln -s ../data_th/forProDiMo.fits.gz . ; cp ~/mcfost/src/forProDiMo/*.in ." ;
      write, f, "      cd data_ProDiMo" ;
      write, f, "      $prodimo > prodimo.log" ;
      write, f, "      cd .." ;
    }
  } // i
  //write, f, "      mv mcfost.log data_th" ;
  write, f, "      mv data_* mcfost.log \"$new_dir\"" ;
  //write, f, "      rm -rf \"$param\"" ;
  // Commande de post-processing
  if (!is_void(post_process)) {
    write, f, "      cd \"$new_dir\"" ;
    for (i=1 ; i<=numberof(post_process) ; i++) {
      write, f, "      "+post_process(i) ;
    } // i
    write, f, "      cd .." ;
  } // post_process
  write, f, "      mv \"$new_dir\" .." ;
  write, f, "      cd .." ;
  write, f, "    fi" ;
  write, f, "done" ;
  write, f, "rm -rf \"$node_dir\"" ;
  close, f ;

  system, "chmod +x "+dir+"/distrib.sh" ;
}

//***************************************************


func stat_grid(grid_name,cluster=,file=) {

    if (is_void(cluster)) cluster="fostino" ;
    if (is_void(file)) file="sed_mc.fits.gz" ;

    write, "Searching for "+file ;
    write, "I found :"

    cmd = "ssh "+cluster+" 'oarstat | grep "+strpart(grid_name,1:14)+" ; "+
      "find "+grid_name+"/ -name "+file+" | wc -l '" ;

    //cmd  ;
    system, cmd ;

}

//***************************************************

func tar_grid(grid_name,cluster=,root_dir=) {

  if (is_void(cluster)) cluster="fostino" ;
  if (is_void(root_dir)) root_dir = "." ;
  dir = root_dir+"/"+grid_name ;

  system, "ssh "+cluster+" 'cd "+dir+"&& ~/Software/bin/merge_job.sh "+
    "&& cd .. && tar czf "+grid_name+".tgz "+grid_name+"'" ;
}

//***************************************************


func visfit(name,dir,slambda,obs,lect=,win=,winP=,color=,clear=,erreur=,noX=,progressive_plot=,output_pdf=,rt=,proba_file=,extension=,correct_dist=,sed_write=,return_chi2=,force_fit=,klambda=,Jy=) {
/* DOCUMENT visfit(name,dir,slambda,obs,lect=,win=,winP=,color=,clear=,erreur=,noX=,progressive_plot=,output_pdf=,rt=,proba_file=,extension=,distance=,sed_write=,klambda=,Jy=)
   name must be the same as the one used to create the grid
   dir is the directory with the results of the grid calculation

   obs.baseline must be in m or klambda
   obs.real  must be in W/m^2 or Jy

   SEE ALSO: get_chi2, analyse_bayesienne
 */

  extern n_params_tot ;

  if (is_void(return_chi2)) return_chi2=0;
  if (is_void(force_fit)) force_fit=0 ;

  if (is_void(extension)) {
    extension="";
  } else {
    extension = "_"+extension ;
  }

  if (is_void(correct_dist)) {
    correct_dist = [1.0] ;
    correct_flux_dist = 1.0 ;
    log_correct_flux_dist = 0.0
    n_D = 1 ;
  } else {
    correct_flux_dist = 1.0/correct_dist^2 ;
    log_correct_flux_dist = log(correct_flux_dist) ;
    n_D = numberof(correct_dist) ;
  }


  f = openb(name+".grid");
  restore, f;
  close, f;


  if (rt) {
    n_incl = mc_grid.ref.map.RT_ntheta ;
  } else {
    n_incl = mc_grid.ref.map.nthet ;
  }

  if (is_void(lect)) lect=1 ;
  if (is_void(win)) win=1 ;
  if (is_void(winP)) winP=win+1 ;
  if (is_void(clear)) clear=0 ;
  if (is_void(color)) color="blue" ;
  if (is_void(rt)) rt=0 ;
  if (is_void(proba_filename)) proba_filename=0 ;
  if (is_void(progressive_plot)) progressive_plot=0 ;
  if (is_void(noX)) noX=0 ;
  if (is_void(output_pdf)) output_pdf=0 ;
  if (noX) {
    progressive_plot=0 ;
    output_pdf = 1 ;
  }
  if (is_void(klambda)) klambda = 0 ;
  if (is_void(Jy)) Jy = 0 ;


  // Dimensions de l'espace des parametres
  D = mc_grid.dims(where(mc_grid.dims > 1)) ;
  n_Parametres = numberof(D) + 1 ; // + 1 pour inclinaison
  n_models = 1 ; for (i=1 ; i<=numberof(D) ; i++) n_models *= D(i) ;
  D2 = [n_Parametres,n_incl] ;
  grow, D2, D ;
  if (n_D > 1) grow, D2, [n_D] ;


  N_max = 31998 ;
  n_Parts = n_models / N_max + 1 ;
  Parts = array(string,n_Parts) ;

  // Definition des sous-dossiers
  if (n_Parts > 1) {
    n_chiffres = floor(log10(n_Parts) + 1);
    format = "%"+swrite(format="%1i",int(n_chiffres))+"i";
    for (i=1 ; i<=n_Parts ; i++) {
      Parts(i) = "Part"+swrite(i,format=format)+"/" ;
    }
  }

  // Lecture du premier modele
  n_chiffres = floor(log10(n_models) + 1);
  format = "%"+swrite(format="%1i",int(n_chiffres))+"i";
  number1=swrite(format=format,1);
  for (j=1 ; j<=n_chiffres -1 ; j++) {
    number1 = streplace(number1,strgrep(" ",number1),"0");
  }

  if (n_Parts > 1) {
    dir1 = dir+Parts(1)+number1
  } else {
    dir1 = dir+number1 ;
  }

  model = open_mcfost(dir1,slambda);
  if (rt) {
    im = model.image_rt ;
  } else {
    im = model.image ;
  }

  // taille pixel image en AU
  pix_size = (model.P.map.size_neb*2.0) / (model.P.map.ny*model.P.map.zoom);

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


  // speed up Fourier transform by selecting optimized size ;
  n = dimsof(im)(2);
  n_fft = fft_best_size(n) ; n_extra = (n_fft -n)/2 ;

  write, "Padding with",n_extra," extra pixels for fft" ;

  // lignes de base
  size = n_fft ;
  center = size/2+1;

  pix_size = pix_size/3600. * pi/180.;
  pix_fft = 1.0/pix_size;
  pix=model.lamb*1e-6*pix_fft;
  baseline=span(0,pix/2,size/2);
  if (klambda==1) baseline = baseline / (model.lamb * 1e-3) ; // baseline en m ou klambda


  if (n_D > 1) n_Parametres++ ;

  if (numberof(obs) - n_Parametres <= 0) {
    "Error: too many parameters !!!";
    write, "N obs=", numberof(obs), "N Param =", n_Parametres ;
    if (!force_fit) {
      return ;
    } else {
      "Forcing fit in 5 sec" ;
      pause, 5000  ;
    }
  }


  n_chiffres = floor(log10(n_models) + 1) ;
  format = "%"+swrite(format="%1i",int(n_chiffres))+"i";

  // Lecture eventuelle
  if (lect) {
    chi2=[];
    chi2=array(double, n_incl, n_models, n_D);

    minchi2 = 1.e300 ;

    // Boucle sur les modeles
    for (k=1 ; k<=n_models ; k++) {

      if (!(k%1000)) write, "Model", k ;

      sk=swrite(format=format,k); for (j=1 ; j<=n_chiffres -1 ; j++) sk = streplace(sk,strgrep(" ",sk),"0"); // Ajout des 0
      ipart = k / N_max + 1 ;


      if (rt) {
        filename = dir+"/"+Parts(ipart)+sk+"/data_"+slambda+"/RT.fits.gz" ;
      } else {
        filename = dir+"/"+Parts(ipart)+sk+"/data_"+slambda+"/MC.fits.gz" ;
      }

       if (IsFile(filename)) {
        image = cfitsRead(filename) ;
        for (i=1 ; i <= n_incl ; i++) {
          im = image(,,i) ;

          // Padding to speed up fft
          im2 = array(float,n_fft,n_fft) ;
          im2(n_extra+1:-n_extra,n_extra+1:-n_extra) = im ;
          im = im2 ; im2 = [] ;

          // fft
          fim=(fft(roll(im),1));

          // visibilités
          if (hor==1) {
            vis = fim(1:size/2,1);
          } else {
            vis = fim(1,1:size/2);
          }
          if (average==1) vis = 0.5 * (fim(1:size/2,1) + fim(1,1:size/2));
          if (Jy==1) vis = vis / (SI.c * 1e-26 / (model.lamb * 1e-6)) ;


          // Interpolation aux lignes de bases observees
          V = interp(vis.re,baseline,obs.baseline) ;

          for (d=1 ; d<= n_D ; d++) {

              C2 = ( (( V - obs.real)/(obs.rerr))^2)(sum) / (numberof(obs) - n_Parametres);

              chi2(i,k,d) = C2 ;
              if (C2 < minchi2) {
                minchi2 = C2 ;
                write, "New best model =", filename, "chi2=", C2 ;
                if (progressive_plot) {
                  window, win ; fma ;    lxy, 1, 0 ;
                  plp, obs.real, obs.baseline, dy=obs.rerr ;

                  if (clear) {limits ; relimits_log;}
                  plg, vis.re, baseline, color=color ;
                }
              }

              if (ieee_test(chi2(i,k,d))) chi2(i,k,d) = 1.0e9;
            } //d
        } // i
      } else {
        chi2(:,k,:) = 1.0e9 ;
      }
    }

    if (rt) {
      f = createb(name+extension+"_rt_chi2_vis"+slambda+".dat");
    } else {
      f = createb(name+extension+"_chi2_vis"+slambda+".dat");
    }
    save, f, chi2;
    close, f;
  } else {
    if (rt) {
      f = openb(name+extension+"_rt_chi2_vis"+slambda+".dat");
    } else {
      f = openb(name+extension+"_chi2_vis"+slambda+".dat");
    }
    restore, f;
    close, f;
  }



  // Meilleur modele
  best = chi2(*)(mnx);
  best = indexof(chi2,best);

  print, "BEST MODEL: ", best , "chi2=", min(chi2);

  sbest=swrite(format=format,best(2)); for (j=1 ; j<=n_chiffres -1 ; j++) sbest = streplace(sbest,strgrep(" ",sbest),"0"); // Ajout des 0
  ipart = best(2) / N_max + 1 ;
  filename = dir+"/"+Parts(ipart)+sbest ;

  write, "Model="+filename, "i=", best(1) ;
  mod = open_mcfost(filename,slambda);


  if (noX) {
    window, win, style="landscape.gs", display="", hcp="hcp" ;
  } else {
    window, win, style="landscape.gs" ;
  }
  if (clear) fma;

  if (rt) {
    im = mod.image_rt(,,best(1),1) ;
  } else {
    im = mod.image(,,best(1),1) ;
  }
  // Padding to speed up fft
  im2 = array(float,n_fft,n_fft) ;
  im2(n_extra+1:-n_extra,n_extra+1:-n_extra) = im ;
  im = im2 ; im2 = [] ;

  fim=(fft(roll(im),1));

  // visibilités
  if (hor==1) {
    vis = fim(1:size/2,1);
  } else {
    vis = fim(1,1:size/2);
  }
  if (average==1) vis = 0.5 * (fim(1:size/2,1) + fim(1,1:size/2));
  if (Jy==1) vis = vis / (SI.c * 1e-26 / (model.lamb * 1e-6)) ;


  fma ; lxy, 1, 0 ;
  plp, obs.real, obs.baseline, dy=obs.rerr ;
  limits ; relimits, 0, 8 ; relimits_log, 85, 100;
  plg, vis.re, baseline, color=color ;

  if (output_pdf) pdf, name+extension+"_Vis.pdf" ;

  Av = 0 ;
  P = _analyse_bayesienne(name,mc_grid,chi2, Av=Av,win=winP,noX=noX,rt=rt,clear=clear,proba_filename=proba_filename,output_pdf=output_pdf,extension=extension) ;


  if (return_chi2) {
    C2 = array(double,D2) ; // Le nbre de dim est n_Parametres
    C2(*) = chi2(*) ;
    return C2 ;
  } else {
    return P ;
  }
}

//***************************************************

func mcfost_grid2ascii(name,Av=,rt=){
/* DOCUMENT mcfost_grid2ascii(name,Av=,rt=)
   write an ASCII file with the parameter and chi2 values
   from the XX.grid and XX_chi2.dat files

   The Av array must be the same as the one used to compute the chi2

   SEE ALSO:
 */

  if (is_void(Av)) Av = 0. ;
  if (is_void(rt)) rt = 0 ;

  n_Av = numberof(Av);
  if (rt) {
    n_incl = mc_grid.ref.map.RT_ntheta ;
  } else {
    n_incl = mc_grid.ref.map.nthet ;
  }


  f = openb(name+".grid");
  restore, f;
  close, f;


  // Dimensions de l'espace des parametres
  D = mc_grid.dims(where(mc_grid.dims > 1)) ;
  n_Parametres = numberof(D) + 1 ; // + 1 pour inclinaison
  if (n_Av > 1) n_Parametres++ ;

  n_models = 1 ; for (i=1 ; i<=numberof(D) ; i++) n_models *= D(i) ;


  f = openb(name+"_chi2.dat");
  restore, f ;
  close, f ;

  D2 = [n_Parametres,n_incl] ;
  grow, D2, D ;
  if (n_Av > 1) grow, D2, [n_Av] ;

  C2 = array(double,D2) ; // Le nbre de dim est n_Parametres

  if (numberof(C2) != numberof(chi2)) {
    write, "ERROR in dimensions, make sure you used the correct number of Av points" ;
    return [] ;
  }
  C2(*) = chi2(*) ;

  j = 0 ;
  X = [] ;
  x_name = "# " ;
  if (n_incl > 1) {
    grow, x_name,  "i" ;
    if (rt) {
      cos_imin= cos(mc_grid.ref.map.RT_imin/180. * pi) ;
      cos_imax= cos(mc_grid.ref.map.RT_imax/180. * pi) ;
      x0 = acos(span(cos_imin,cos_imax,n_incl+1)(zcen)) * 180./pi ;
    } else {
      x0 = acos(span(1,0,n_incl+1)(zcen)) * 180./pi ;
    }
    X = x0 ;
  }
  //p = P(,*)(,sum) ;

  for(i=1; i<= numberof(mc_grid.dims); i++) {
    if (mc_grid.dims(i) > 1) {
      j++ ;
      grow, x_name, mc_grid.names(i) ;
      if (j==1) {
        x1 = mc_grid.values(i,1:mc_grid.dims(i)) ;
      } else if (j==2) {
        x2 = mc_grid.values(i,1:mc_grid.dims(i)) ;
      } else if (j==3) {
        x3 = mc_grid.values(i,1:mc_grid.dims(i)) ;
      } else if (j==4) {
        x4 = mc_grid.values(i,1:mc_grid.dims(i)) ;
      } else if (j==5) {
        x5 = mc_grid.values(i,1:mc_grid.dims(i)) ;
      } else if (j==6) {
        x6 = mc_grid.values(i,1:mc_grid.dims(i)) ;
      } else if (j==7) {
        x7 = mc_grid.values(i,1:mc_grid.dims(i)) ;
      } else if (j==8) {
        x8 = mc_grid.values(i,1:mc_grid.dims(i)) ;
      }
    }
  }

  // Av
  if (numberof(Av) > 1) {
    grow, x_name,  "A_V" ;
    if (n_Parametres == 2) x1 = Av ;
    if (n_Parametres == 3) x2 = Av ;
    if (n_Parametres == 4) x3 = Av ;
    if (n_Parametres == 5) x4 = Av ;
    if (n_Parametres == 6) x5 = Av ;
    if (n_Parametres == 7) x6 = Av ;
    if (n_Parametres == 8) x7 = Av ;
    if (n_Parametres == 9) x8 = Av ;
  }


  // Writing ascii file
  f = open(name+".txt","w") ;
  grow, x_name, "chi2" ;
  write, f, x_name ;
  if (n_Parametres == 10) {
    for (i0=1 ; i0<=numberof(x0) ; i0++) {
      for (i1=1 ; i1<=numberof(x1) ; i1++) {
        for (i2=1 ; i2<=numberof(x2) ; i2++) {
          for (i3=1 ; i3<=numberof(x3) ; i3++) {
            for (i4=1 ; i4<=numberof(x4) ; i4++) {
              for (i5=1 ; i5<=numberof(x5) ; i5++) {
                for (i6=1 ; i6<=numberof(x6) ; i6++) {
                  for (i7=1 ; i7<=numberof(x7) ; i7++) {
                    for (i8=1 ; i8<=numberof(x8) ; i8++) {
                      for (i9=1 ; i9<=numberof(x9) ; i9++) {
                        write, f, x0(i0), x1(i1), x2(i2), x3(i3), x4(i4), x5(i5), x6(i6), x7(i7), x8(i8), x9(i9), C2(i0,i1,i2,i3,i4,i5,i6,i7,i8,i9) ;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  } // n_Parametres

  if (n_Parametres == 9) {
    for (i0=1 ; i0<=numberof(x0) ; i0++) {
      for (i1=1 ; i1<=numberof(x1) ; i1++) {
        for (i2=1 ; i2<=numberof(x2) ; i2++) {
          for (i3=1 ; i3<=numberof(x3) ; i3++) {
            for (i4=1 ; i4<=numberof(x4) ; i4++) {
              for (i5=1 ; i5<=numberof(x5) ; i5++) {
                for (i6=1 ; i6<=numberof(x6) ; i6++) {
                  for (i7=1 ; i7<=numberof(x7) ; i7++) {
                    for (i8=1 ; i8<=numberof(x8) ; i8++) {
                      write, f, x0(i0), x1(i1), x2(i2), x3(i3), x4(i4), x5(i5), x6(i6), x7(i7), x8(i8), C2(i0,i1,i2,i3,i4,i5,i6,i7,i8) ;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  } // n_Parametres

  if (n_Parametres == 8) {
    for (i0=1 ; i0<=numberof(x0) ; i0++) {
      for (i1=1 ; i1<=numberof(x1) ; i1++) {
        for (i2=1 ; i2<=numberof(x2) ; i2++) {
          for (i3=1 ; i3<=numberof(x3) ; i3++) {
            for (i4=1 ; i4<=numberof(x4) ; i4++) {
              for (i5=1 ; i5<=numberof(x5) ; i5++) {
                for (i6=1 ; i6<=numberof(x6) ; i6++) {
                  for (i7=1 ; i7<=numberof(x7) ; i7++) {
                      write, f, x0(i0), x1(i1), x2(i2), x3(i3), x4(i4), x5(i5), x6(i6), x7(i7), C2(i0,i1,i2,i3,i4,i5,i6,i7) ;
                  }
                }
              }
            }
          }
        }
      }
    }
  } // n_Parametres

  if (n_Parametres == 7) {
    for (i0=1 ; i0<=numberof(x0) ; i0++) {
      for (i1=1 ; i1<=numberof(x1) ; i1++) {
        for (i2=1 ; i2<=numberof(x2) ; i2++) {
          for (i3=1 ; i3<=numberof(x3) ; i3++) {
            for (i4=1 ; i4<=numberof(x4) ; i4++) {
              for (i5=1 ; i5<=numberof(x5) ; i5++) {
                for (i6=1 ; i6<=numberof(x6) ; i6++) {
                      write, f, x0(i0), x1(i1), x2(i2), x3(i3), x4(i4), x5(i5), x6(i6), C2(i0,i1,i2,i3,i4,i5,i6) ;
                }
              }
            }
          }
        }
      }
    }
  } // n_Parametres

    if (n_Parametres == 6) {
    for (i0=1 ; i0<=numberof(x0) ; i0++) {
      for (i1=1 ; i1<=numberof(x1) ; i1++) {
        for (i2=1 ; i2<=numberof(x2) ; i2++) {
          for (i3=1 ; i3<=numberof(x3) ; i3++) {
            for (i4=1 ; i4<=numberof(x4) ; i4++) {
              for (i5=1 ; i5<=numberof(x5) ; i5++) {
                write, f, x0(i0), x1(i1), x2(i2), x3(i3), x4(i4), x5(i5), C2(i0,i1,i2,i3,i4,i5) ;
              }
            }
          }
        }
      }
    }
  } // n_Parametres


  close, f ;
}
