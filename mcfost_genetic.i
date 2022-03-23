//"Loading genetic algorithm" ;

/*
http://www.talkorigins.org/faqs/genalg/genalg.html#what:selection

Hierarchical selection: Individuals go through multiple rounds of selection each generation. Lower-level evaluations are faster and less discriminating, while those that survive to higher levels are evaluated more rigorously. The advantage of this method is that it reduces overall computation time by using faster, less selective evaluation to weed out the majority of individuals that show little or no promise, and only subjecting those who survive this initial test to more rigorous and more computationally expensive fitness evaluation.

*/
Genetic_file = "~/yorick/init/mcfost_genetic.i" ;
Genetic_N_explored_param  = (Genetic_N_explored_param ? Genetic_N_explored_param : 1) ; // will be updated

default_cluster = "gofree" ;


struct GeneticModel {
  int generation ;
  int id ;
  float chi2 ;
  float Proba ;
  float Parameters(Genetic_N_explored_param) ;
  int parents(2) ;
  int parameter_muted ;
  string dir ;
  int out ;
} ;

struct ClusterConfig {
  int n_nodes ;
  int omp_num_threads ;
  int walltime ;
  int waiting_time ;
  string running_dir ;
  string grid_name ;
  int n_options ;
  string options(10) ;
} ;

struct GeneticRun {
  ClusterConfig FostiConfig ;
  float pmut ;
  float pmut0 ;
  float pmut_max ;
  float pmut_min ;
  int generation_max ;
  int generation_reached ;
  int elitism ;
  int n_models ;
  string simu_config_file ;
  float seed ;
}

//************************************************************************

  func mcfost_genetic(grid_name,parameter_file,fitting_routine, simu_config_file=, n_nodes=,omp_num_threads=,walltime=,waiting_time=,root_dir=, options=,n_models=,pmut=,pmut_max=,pmut_min=,generation_max=,elitism=,no_mcfost_compute=,progressive_plot=, rt=) {
/* DOCUMENT mcfost_genetic(grid_name,parameter_file,fitting_routine,

   simu_config_file=,n_nodes=,omp_num_threads=,walltime=,root_dir=,

   options=,n_models=,pmut=,pmut_max=,pmut_min=,generation_max=,elitism=)

   Default values :
   - simu_config_file : root_dir+"/"+grid_name+"/config.genetic"
   - n_nodes : 20
   - omp_num_threads : 8
   - walltime (per generation) : 24h
   - waiting_time (to check if generation if finished) : 30s
   - root_dir : ~/Simus/Genetic/
   - options : none
   - n_models : 100
   - pmut : 0.2
   - pmut_max : 0.5
   - pmut_min : 1e-2
   - generation_max : 100
   - elitism : 1  (the value indicates the number of best parents kept in the next parents'pool)
   - no_mcfost_compute = 0 : skip the mcfost calculations, for fitting routines that do not need mcfost
   - progressive_plot = 1

   SEE ALSO:
 */

  extern Genetic_N_explored_param, Genetic_file ;

  log_file = open("Genetic_"+grid_name+".log","w") ;

  if (is_void(rt)) rt=0 ;
  Genetic_rt = rt ;

  // GA config
  if (is_void(pmut)) pmut=0.2 ;
  if (is_void(pmut_max)) pmut_max=0.5 ;
  if (is_void(pmut_min)) pmut_min=1e-2 ;
  if (is_void(generation_max)) generation_max = 100 ;
  if (is_void(elitism)) elitism = 1 ;
  if (is_void(n_models)) n_models=100 ;
  if (is_void(no_mcfost_compute)) no_mcfost_compute=0 ;
  if (is_void(progressive_plot))  progressive_plot=1 ;


  // Cluster configuration
  if (is_void(omp_num_threads)) omp_num_threads = 8 ;
  if (is_void(walltime)) walltime = 24 ;
  if (is_void(waiting_time)) waiting_time = 30 ;
  if (is_void(n_nodes)) n_nodes = 20 ;
  if (is_void(root_dir)) {
    running_dir = "~/Simus/Genetic/"+grid_name ;
  } else {
    running_dir = root_dir+"/"+grid_name ;
  }
  if (is_void(simu_config_file)) simu_config_file = running_dir+"/config.genetic" ;
  n_options = numberof(options) ;


  // Read configuration file
  SimuConfig = OpenASCII(simu_config_file,prompt=0) ;
  Genetic_N_explored_param = numberof(SimuConfig) ;
  // Recharge la structure pour redefinir GenPop
  include, Genetic_file, 1 ;

  // GA structure
  GA = array(GeneticRun) ;
  GA.pmut = pmut ;
  GA.pmut0 = pmut ;
  GA.pmut_max = pmut_max ;
  GA.pmut_min = pmut_min ;
  GA.generation_max = generation_max ;
  GA.generation_reached = 0 ;
  GA.elitism = elitism ;
  GA.n_models = n_models ;
  GA.simu_config_file = simu_config_file ;

  // Cluster config
  FostiConfig = array(ClusterConfig) ;
  FostiConfig.grid_name = grid_name ;
  FostiConfig.n_nodes = n_nodes ;
  FostiConfig.omp_num_threads = omp_num_threads ;
  FostiConfig.walltime = walltime ;
  FostiConfig.waiting_time = waiting_time ;
  FostiConfig.running_dir = running_dir ;
  n_options = numberof(options) ;
  FostiConfig.n_options = n_options ;
  for (i=1 ; i<=n_options ; i++) {
    FostiConfig.options(i) = options(i) ;
  }
  GA.FostiConfig = FostiConfig ;


  // Read refence parameter file
  V = 2.20 ;
  ref = read_params(parameter_file) ;
  if (abs(V - ref.simu.version) > 1.e-4) {
    write, "Genetic : Version ", V, "needed. Exiting !!!";
    return;
  }

  // Mise a jour des structures MCFOST
  include, struct_file ;

  parents_pop = [] ;

  BEST = [] ;
  MEDIAN = [] ;
  PMUT = [] ;

  pldefault, style="landscape.gs", marks=0, width=1, palette="heat.gp", legends=0, dpi=75;


  // Routine de fit
  if (typeof(fitting_routine) == "string") {
    post_process = fitting_routine ;      // commande excute sur chaque noeud
    fitting_routine = read_genetic_chi2 ; // on ne lit que le chi2
  } else {
    post_process = [] ;
  }


  for (generation=1 ; generation <= generation_max ; generation++) {

    write, "------------------------------------------------------" ;
    write, "Generation", generation ;
    write, log_file, "------------------------------------------------------" ;
    write, log_file, "Generation", generation ;

    // Creation du dossier pour la generation
    dir=running_dir+"/"+"gen"+swrite(format="%4i",generation) ;
    for (j=1 ; j<=3 ; j++) dir = streplace(dir,strgrep(" ",dir),"0");
    system, "rm -rf "+dir ; // removing eventual previous run
    mkdirp, dir ;

    // Create the new population
    write, "Creating the new population ..." ;
    write, log_file, "Creating the new population ..." ;
    if (generation == 1) {
      write, "Initializing parameter files for first generation" ;
      write, log_file, "Initializing parameter files for first generation" ;
      GA_seed = randomize() ;  GA.seed = GA_seed ;
      write, "Random seed =", GA_seed ;
      write, log_file, "Random seed =", GA_seed ;
    }
    pop = make_population(SimuConfig, parents_pop, generation, n_models, pmut, dir) ;


    if (!no_mcfost_compute) {
      // Create the corresponding parameter files
      write, "Creating the parameter files ..." ;
      write, log_file, "Creating the parameter files ..." ;
      make_parameter_files, SimuConfig, ref, pop, generation ;  // todo : a merger avec l'etape precedente ???

      // Run the models --> system call
      write, "Computing the MCFOST models ..." ;
      write, log_file, "Computing the MCFOST models ..." ;
      compute_mcfost_models, generation, FostiConfig, options, post_process, dir ;
    }

    // Compute the chi2
    compute_fits, fitting_routine, pop, generation, dir ;

    //write, pop.chi2 ;

    // enregistrement de la pop en fichier binaire yorick
    f = createb(dir+".bin") ;
    save, f, pop ;
    close, f ;

    // Select the models that will be parents for the next generation of models
    // this population is ranked
    parents_pop = make_parents_pop(parents_pop, pop, elitism=elitism) ;

    // Log of results
    write, "Best and median chi2", parents_pop(1).chi2, parents_pop(n_models/2).chi2 ;
    write, log_file, "Best and median chi2", parents_pop(1).chi2, parents_pop(n_models/2).chi2 ;
    write, "Best model from generation:", parents_pop(1).generation, "id=", parents_pop(1).id  ;
    write, log_file, "Best model from generation:", parents_pop(1).generation, "id=", parents_pop(1).id  ;

    // adjust mutation rate
    pmut = adjust_mutation_rate(parents_pop, pmut, pmut_max, pmut_min) ;
    write, "Adjusting mutation rate to:", pmut ;
    write, log_file, "Adjusting mutation rate to:", pmut ;


    // Log of results
    write, "Current best parameters :" ;
    write, log_file, "Current best parameters :"
    write, SimuConfig.parameter, parents_pop(1).Parameters ;
    write, log_file, SimuConfig.parameter, parents_pop(1).Parameters ;


    // Saving the state of the GA
    GA.generation_reached = generation ;
    GA.pmut = pmut ;
    f = createb(running_dir+"/GA.bin") ;
    save, f, GA ;
    close, f ;

    // Plots
    grow, BEST, parents_pop(1) ;
    grow, MEDIAN, [parents_pop(1:n_models/2)] ;
    grow, PMUT, pmut ;

    if ( (generation == generation_max) || progressive_plot )  {
      write, "Generating plots ..." ;
      write, log_file, "Generating plots ..." ;
      // Chi2 en fct generation
      window, 1, display="", hcp="hcp", style="landscape.gs" ;
      fma ; limits ; lxy, 0, 1 ;
      plg, BEST.chi2, color="red" ;
      plg, MEDIAN(0,).chi2 ;
      relimits_log, 100, 50.0 ;
      xytitles, "generation", "!c^2^" ;

      pltr, "  Best !c^2 = "+swrite(parents_pop(1).chi2), 70, 90 ;
      pltr, "Median !c^2 = "+swrite(parents_pop(n_models/2).chi2), 70, 80 ;

      pdf, "chi2" ; system, "mv chi2.pdf "+running_dir ;

      // Taux de mutation en fct generation
      window, 2, display="", hcp="hcp", style="landscape.gs" ;
      fma ; limits ; lxy, 0, 1 ;
      plg, PMUT ;
      relimits_log, 100, 50.0 ;
      xytitles, "generation", "Mutation rate" ;
      pdf, "mutation_rate" ; system, "mv mutation_rate.pdf "+running_dir ;

      // Distribution des Chi2 et analyse bayesienne
      //window, 3, display="", hcp="hcp", style="landscape.gs";
      //GA_stat, grid_name,running_dir=running_dir, win=3 ;
      // pdf bug pour le moment ????
      //pdf, "chi2_distrib" ; system, "mv chi2_distrib.pdf "+running_dir ;
      window, 4, display="", hcp="hcp", style="landscape.gs" ;
      //pdf, "bayesian_analysis" ; system, "mv bayesian_analysis.pdf "+running_dir ;
      write, "Done." ;
      write, log_file, "Done." ;
    }

  } // generation
  close, f ;

  // Population finale
  pop = parents_pop ;

  return pop ;
}

//************************************************************************

func make_population(SimuConfig,parents_pop,generation,n_models,pmut, dir) {
  /*
    Creates a table of models to compute for the generation specified.

    If this is the generation 0, then the parameters are sampled in an
    unbiased way within the ranges specified by the user. Otherwise, this
    method uses results from previous generations to determine which
    models to run.
  */

  extern Genetic_N_explored_param ;

  mutation_type = 2 ; // 1 : uniform sur 1 parametre, 2 : gaussien sur tous les parametres

  kfrac = 0.2 ;
  p_tournament = 0.9 ;

  pop = array(GeneticModel,n_models) ;
  pop.out = 0 ; // pour eMCMC

  if (generation == 1) {
    pop.generation() = 1 ;
    pop.id() = indgen(n_models) ;

    // Creating the array of parameter values
    for (p=1 ; p<=Genetic_N_explored_param ; p++) {
      CP = SimuConfig(p) ;
      if (CP.mode == "linear") {
        values =  CP.vmin + random(n_models) *  (CP.vmax - CP.vmin) ;
      } else if (CP.mode == "log") {
        values =  10^( log10(CP.vmin) + random(n_models) *  (log10(CP.vmax) - log10(CP.vmin)) ) ;
      } else {
        write, "Unknown mode:", CP.mode ;
        write, "Exiting" ;
        error ;
      }
      pop.Parameters(p,) = values ;
    }

  } else { // generation > 1
    n_children = n_models ;
    // Creation new population
    for (i = 1 ; i<= n_children ; i++) {
      pop(i).generation = generation ;
      pop(i).id = i ;

      //father = tournament_select(pop,kfrac,p_tournament) ;
      father = roulette_wheel_select(numberof(pop)); // queue bcp plus longue : mieux pour les proba bayesienne
      do {
        //mother = tournament_select(pop,kfrac,p_tournament) ;
        mother = roulette_wheel_select(numberof(pop)) ;
      } while (mother == father) ;


      pop(i).parents = [father,mother] ;

      //pop(i).parents

      //write, "father", parents_pop(father).Parameters
      //write, "mother", parents_pop(mother).Parameters

      // Crossover all the parameters from the parents
      for (p=1 ; p<=Genetic_N_explored_param ; p++) {
        f = random() ;

        if (SimuConfig(p).mode == "linear") {
          value = parents_pop(father).Parameters(p) * f + parents_pop(mother).Parameters(p) * (1.-f) ;
        } else {
          value = exp( log(parents_pop(father).Parameters(p)) * f + log(parents_pop(mother).Parameters(p)) * (1.-f) ) ;
        }
        pop(i).Parameters(p) = value ;
      } // p loop

      //write, "kid", pop(i).Parameters

      // Apply eventual mutation
      pop(i).parameter_muted = 0 ;

      if (random() < pmut) {
        if (mutation_type == 1) {
          // Select a gene for mutation & apply mutation
          p = int_random(Genetic_N_explored_param) ; CP = SimuConfig(p) ;
          pop(i).parameter_muted = p ;
          f = random() ;
          if (CP.mode == "linear") {
            value =  CP.vmin + f *  (CP.vmax - CP.vmin) ;
          } else {
            value =  10^( log10(CP.vmin) + f *  (log10(CP.vmax) - log10(CP.vmin)) ) ;
          }
          pop(i).Parameters(p) = value ;
        } else { // mutation_type=2
          // Gaussian mutation of all parameters
          r = Gauss_random(Genetic_N_explored_param) ;
          for (p=1 ; p<=Genetic_N_explored_param ; p++) {
            CP = SimuConfig(p) ;
            if (CP.mode == "linear") {
              sigma = 0.1 * (CP.vmax - CP.vmin)  ;
              pop(i).Parameters(p) += r(p) * sigma ;
            } else {
              sigma =  0.1 *(log10(CP.vmax) - log10(CP.vmin)) ;
              pop(i).Parameters(p) =  10^(log10(pop(i).Parameters(p)) + r(p) * sigma)  ;
            }
            // Make sure we don;t go above the limits
            if (pop(i).Parameters(p) > CP.vmax) pop(i).Parameters(p) = CP.vmax ;
            if (pop(i).Parameters(p) < CP.vmin) pop(i).Parameters(p) = CP.vmin

          }
        }

        //write, "mutation", pop(i).Parameters
      }
    } // boucle t

  } // test 1er generation

  // Dossier dans lequel va etre la simu
  n_chiffres = floor(log10(n_models) + 1);
  format = "%"+swrite(format="%i",int(n_chiffres))+"i";
  for (i=1 ; i<=n_models ; i++) {
    number=swrite(format=format,i);
    // Ajout des 0
    for (j=1 ; j<=n_chiffres -1 ; j++) {
      number = streplace(number,strgrep(" ",number),"0");
    }
    pop(i).dir = dir+"/"+number ;
  }



  return pop ;

}

//************************************************************************

func set_mcfost_parameters(P,Pname,Pvalues) {
  // Set 1 parameter to different values for a bunch of parameter files

  // test type de P

  // test dimension P et Pvalues

  //diskmass=,H0=,Rin=,Rout=,edge=,beta=,surf=,amin=,amax=,aexp=,


  // On fait les test sur Pname
  if (Pname == "amin") {
    P.dust_pop.amin(1,1,) = Pvalues ;
  } else if (Pname == "amax") {
    P.dust_pop.amax(1,1,) = Pvalues ;
  } else if (Pname == "aexp") {
    P.dust_pop.aexp(1,1,) = Pvalues ;
  } else if (Pname == "porosity") {
    P.dust_pop.porosity(1,1,) = Pvalues ;

    // grains dans la 2eme zone
  } else if (Pname == "amin_2") {
    P.dust_pop.amin(2,1,) = Pvalues ;
  } else if (Pname == "amax_2") {
    P.dust_pop.amax(2,1,) = Pvalues ;
  } else if (Pname == "aexp_2") {
    P.dust_pop.aexp(2,1,) = Pvalues ;
  } else if (Pname == "porosity_2") {
    P.dust_pop.porosity(2,1,) = Pvalues ;

  } else if (Pname == "strat") {
    P.dust().strat = Pvalues ;

  } else if (Pname == "viscosity") {
    P.simu.viscosity = Pvalues ;

  } else if (Pname == "diskmass") {
    P.zones(1,).mass = Pvalues ;
  } else if (Pname == "gas_to_dust") {
    P.zones(1,).gas_to_dust = Pvalues ;
  } else if (Pname == "H0") {
    P.zones(1,).Ho = Pvalues ;
  } else if (Pname == "R0") {
    P.zones(1,).Ro = Pvalues ;
  } else if (Pname == "Rin") {
    P.zones(1,).rin = Pvalues ;
  } else if (Pname == "Rout") {
    P.zones(1,).rout = Pvalues ;
  } else if (Pname == "Rc") {
    P.zones(1,).Rc = Pvalues ;
  } else if (Pname == "edge") {
    P.zones(1,).edge = Pvalues ;
  } else if (Pname == "r_edge") {
    P.zones(1,).edge = Pvalues *  P.zones(1,).rin;
  } else if (Pname == "beta") {
    P.zones(1,).beta = Pvalues ;
  } else if ( (Pname == "surf") || (Pname == "alpha") ) {
    P.zones(1,).dens = Pvalues ;
  } else if (Pname == "gamma_exp") {
    P.zones(1,).gamma_exp = Pvalues ;

  } else if (Pname == "diskmass_2") {
    P.zones(2,).mass = Pvalues ;
  } else if (Pname == "H0_2") {
    P.zones(2,).Ho = Pvalues ;
  } else if (Pname == "R0_2") {
    P.zones(2,).Ro = Pvalues ;
  } else if (Pname == "Rin_2") {
    P.zones(2,).rin = Pvalues ;
  } else if (Pname == "Rout_2") {
    P.zones(2,).rout = Pvalues ;
  } else if (Pname == "edge_2") {
    P.zones(2,).edge = Pvalues ;
  } else if (Pname == "r_edge_2") {
    P.zones(2,).edge = Pvalues *  P.zones(2,).rin;
  } else if (Pname == "beta_2") {
    P.zones(2,).beta = Pvalues ;
  } else if ( (Pname == "surf_2") || (Pname == "alpha_2") ) {
    P.zones(2,).dens = Pvalues ;


  } else if ( (Pname == "component_volume_fraction2") ) {
    P.dust_pop.component_volume_fraction(2,1,1,) = Pvalues ;
    P.dust_pop.component_volume_fraction(1,1,1,) = 1.0 - Pvalues ;

  } else if ( (Pname == "specie_mass_fraction2") ) {
    P.dust_pop.mass_fraction(1,2,) = Pvalues ;
    P.dust_pop.mass_fraction(1,1,) = 1.0 - Pvalues ;

  } else if ( (Pname == "component_volume_fraction3") ) {
    P.dust_pop.component_volume_fraction(3,1,1,) = Pvalues ;
    P.dust_pop.component_volume_fraction(2,1,1,) *= (1.0 - Pvalues) ;
    P.dust_pop.component_volume_fraction(1,1,1,) *= (1.0 - Pvalues) ;

  } else if ( (Pname == "specie_mass_fraction3") ) {
    P.dust_pop.mass_fraction(1,3,) = Pvalues ;
    P.dust_pop.mass_fraction(1,2,) *= (1.0 - Pvalues) ;
    P.dust_pop.mass_fraction(1,1,) *= (1.0 - Pvalues) ;

  } else if (Pname == "vturb") {
    for (imol=1 ; imol <= P(1).simu.nmol ; imol++) {
      P.mol.vturb(imol,) = Pvalues ;
    }

  } else if (Pname == "abundance1") {
    P.mol(1,).abundance = Pvalues ;
  } else if (Pname == "abundance2") {
    P.mol(2,).abundance = Pvalues ;
  } else if (Pname == "abundance3") {
    P.mol(3,).abundance = Pvalues ;


  } else if (Pname == "Rstar") {
    P.stars(1,).R = Pvalues ;
  } else if (Pname == "Tstar") {
    P.stars(1,).T = Pvalues ;
  } else if (Pname == "fUV") {
    P.stars(1,).fUV =  Pvalues ;

    // TODO pour tous les parametres !!!!

  } else {
    write, "Unknown parameter name:", Pname ;
    write, "Exiting" ;
    return [] ;
  }

}

//************************************************************************

func make_parents_pop(parents_pop,pop,elitism=) {

  if (is_void(elitism)) elitism=0 ;
  n = numberof(pop) ;   POP = pop ;

  // Elitism : we keep the N best elements of the previous parent pop with N = elitism ;
  if ( (elitism > 0) & !is_void(parents_pop)) grow, POP, parents_pop(1:elitism) ;

  // Ordering the population
  POP = POP(sort(POP.chi2)) ;

  // Elitism : we delete the worst elements
  if (elitism) POP = POP(1:n) ;

  return POP ;
}

//************************************************************************

func tournament_select(pop, k_frac, p) {
  /*
  Tournament selection routine
  Good parameters for getting ~10% are k_frac=0.2 and p=0.9

  TODO : a finir
  */

  n_models = numberof(pop) ;

  k = int(n_models * k_frac) ;
  if (k <= 0) {
    write, "k_frac is too small" ;
    return ;
  }

  model_id = indgen(numberof(pop)) ;

  prob = p * (1-p)^indgen(k) ;
  prob = prob / sum(prob) ;

  // Tournament draw : pick k random elements in the population
  players = int_random(n_models,n=k) ;
  // sort the players according to fitness ranking
  order = sort(pop(players).chi2) ;
  players = players(order) ;

  somme = 0. ;
  r = random() ;
  for (i=1 ; i<= k ; i++) {
    somme += prob(i) ;
    if (somme > r) {
      return players(i) ;
    }
  }


}




//************************************************************************

func roulette_wheel_select(n_models) {

  // proba inversement proportionelle au rang
  // models must be ranked here
  proba = 1./indgen(n_models) ;
  proba = proba / sum(proba) ;

  somme = 0. ;
  r = random() ;
  for (i=1 ; i<= n_models ; i++) {
    somme += proba(i) ;
    if (somme > r) {
      return i ;
    }
  }

}

//************************************************************************

func adjust_mutation_rate(pop,pmut,pmut_max,pmut_min) {
/* DOCUMENT adjust_mutation_rate(pop,pmut)
     pop must be ordered here !!!
   SEE ALSO:
 */

  // Parameters
  rdif_low  = 0.05 * 2;
  rdif_high = 0.25 * 2;
  delta = 1.5 ;

  n_models = numberof(pop) ;

  // Ajustement based on fitness differential
  denominator = (pop.chi2(n_models/2) + pop.chi2(1)) ;
  if (denominator == 0.) {
    write, "*********************************************" ;
    write, "You have the perfect model !!! (really ?????)" ;
    write, "Stopping the GA" ;
    write, "*********************************************" ;
    error ;
  }

  rdif = abs(pop.chi2(n_models/2) - pop.chi2(1)) / denominator ;

  if (rdif <= rdif_low) { // increase mutation rate
    pmut = min(pmut_max,pmut * delta) ;
  } else if (rdif >= rdif_high) { // decrease mutation rate
    pmut = max(pmut_min,pmut/delta) ;
  }

  return pmut ;
}

//************************************************************************

func make_parameter_files(SimuConfig, ref, pop, generation) {

  // make_parameter_files, SimuConfig, ref, pop, generation, dir ;  // todo : a merger avec l'etape precedente ???

  // Creation du tableau de structure de parametres
  n_Para = numberof(SimuConfig)
  n_models = numberof(pop) ;
  Params = ref(-:1:n_models) ;
  for (p=1 ; p<=n_Para ; p++) {
    set_mcfost_parameters, Params, SimuConfig(p).parameter, pop.Parameters(p,) ;
  }

  // Creation des fichiers de parametres
  for (i=1 ; i<=n_models ; i++) {
    filename = pop(i).dir+".par" ;
    write_params, Params(i), filename ;
  }

  return ;
}


//************************************************************************

func compute_mcfost_models(generation, FostiConfig, options, post_process, dir) {

  // Creation du script OAR
  write_OAR_script, omp_num_threads, options=options, post_process=post_process, dir=dir ;

  // Lancement de la grille
  write, "Submitting jobs ..." ;
  //ret = yocoSystem("repeat "+swrite(FostiConfig.n_nodes,format="%i")+" oarsub -n \""+FostiConfig.grid_name+"\" -lnodes=1,walltime="+swrite(FostiConfig.walltime,format="%i")+" distrib.sh") ;

  yodir = get_cwd() ;

  cd, dir ;
  ret = yocoSystem("for i in $(seq 1 "+swrite(FostiConfig.n_nodes,format="%i")+"); do oarsub -n \""+FostiConfig.grid_name+"\" -lnodes=1,walltime="+swrite(FostiConfig.walltime,format="%i")+" ./distrib.sh ; done") ;
  cd, yodir ;

  if (ret != 0) {
    yocoError("Submission failed !!!!", "", 2)
  }

  // Test si la grille est terminee, on ne fait rien tant que tous les jobs ne sont pas termines
  jobs_running = n_nodes ;
  buf = int() ;
  while (jobs_running) {
    write, "Gen", generation,":",jobs_running, "jobs running, waiting ", FostiConfig.waiting_time, " sec ..." ;
    pause, FostiConfig.waiting_time * 1000 ; // 30 sec par defaut
    f= popen("oarstat | grep "+strpart(FostiConfig.grid_name,1:14)+"| wc -l",0) ;
    read, f, buf ; if (!is_void(buf)) jobs_running = buf  ;
    close, f ;
  }

  write, "Model computation done for generation ",swrite(generation,format="%i") ;
}

//************************************************************************

func compute_fits(fitting_routine, pop, generation, dir) {

  N = numberof(pop) ;
  fit = array(float,N,2) ;

  write, (1.0*numberof(where(pop.out)))/N, "models out" ;

  for (i=1 ; i<=N ; i++) {
    if (pop(i).out) {
      fit(i,1) = 1e9 ;
      fit(i,2) = 0. ;
    } else {
      // chi2 et proba marginale sur parametres explores par routine de fit
      fit(i,) = fitting_routine(pop(i)) ;
      if (fit(i,1) == 0.) {
        fit(i,1) = 1e9 ;
        fit(i,2) = 0. ;
      }
    }
  }
  pop.chi2 = fit(,1) ;
  pop.Proba = fit(,2) ; // proba relative non normalisee
}

//************************************************************************

func GA_stat(grid_name,&best_model,running_dir=,color=,last_color=,win=,generation_max=,nx=,cluster=,clear=,symbol=,size=,fill=) {

  extern default_cluster ;

  if (is_void(running_dir)) running_dir = "~/Simus/Genetic/"+grid_name ;

  if (is_void(color)) {
    color="black" ;
    last_color="red" ;
  }
  if (is_void(last_color)) last_color=color
  if (is_void(win)) win = 0;
  if (is_void(nx)) nx = 3 ;
  if (is_void(cluster)) cluster=default_cluster ;
  if (is_void(clear)) clear=0 ;
  pltitle_height = 12 ;

  if (is_void(symbol)) symbol=6 ;
  if (is_void(size)) size=1 ;
  if (is_void(fill)) fill=0 ;

  Config = OpenASCII(running_dir+"/config.genetic",prompt=0) ;
  N_params = numberof(Config) ; // TODO : parametres virtuels : i et av par exemple
  window, win ;
  if (clear) {
    fma ;
  }
  gs_nm, win, nx, int(ceil((N_params*1.0)/nx)), landscape=1, dx=0.02, dy = 0.07, square=0 ;

  f = openb(running_dir+"/GA.bin") ;
  restore, f, GA ;
  close, f ;

  if (is_void(generation_max))   generation_max = GA.generation_reached ;
  write, generation_max, "generations found  --> ", generation_max * GA.n_models, "models" ;


  XX = [] ;
  P = [] ;
  chi2 = [] ;

  minchi2 = 1e300 ;
  for (generation=1 ; generation <= generation_max ; generation++) {
    pop = [] ;
    dir=running_dir+"/"+"gen"+swrite(format="%4i",generation) ;
    for (j=1 ; j<=3 ; j++) dir = streplace(dir,strgrep(" ",dir),"0");

    write, "------------------------------------" ;
    write, "reading "+dir+".bin" ;
    f = openb(dir+".bin") ;
    restore, f, pop ;
    close, f ;


    ou = pop.chi2(mnx) ; // best model of this generation
    if (pop.chi2(ou) < minchi2) {
      minchi2 = pop.chi2(ou) ;
      BEST = pop(ou) ;
    }
    write, "Current best model :" ;
    BEST ;


    //stat, pop.Parameters(1,)

    if (generation == generation_max) color=last_color ;
    for (iParameter=1 ; iParameter<=N_params ; iParameter++) {
      plsys, iParameter ;
      plp, pop.chi2, pop.Parameters(iParameter,), color=color, symbol=symbol, size=size, fill=fill  ;
    }

    /*
    window, 11 ;
    plp, pop.Parameters(1,) ;
    window, win
    */

    // TODO : avoid grow
    grow, chi2, pop.chi2 ;
    grow, P, pop.Proba ;
    grow, XX, pop.Parameters ;
  }

  // Plot axes
  for (iParameter=1 ; iParameter<=N_params ; iParameter++) {
    plsys, iParameter ;
    if (Config(iParameter).mode == "log") {
      lxy, 1, 1 ;
    } else {
      lxy, 0, 1 ;
    }
    xname = Config(iParameter).parameter ;
    for (i=1 ; i<=5 ; i++) xname = streplace(xname,strgrep("_", xname), " ") ;
    if (iParameter%nx == 1) {
      xytitles, xname, "!c^2", [0.02,0.02] ;
    } else {
      xytitles, xname, "", [0.02,0.02] ;
    }
    limits ; l = limits() ;
    limits, , , 0.9 * l(3), 2* l(3) ;
  }


  // Fichier pour eMCMC
  write, "Writing config file for eMCMC" ;
  f = open(running_dir+"/config.eMCMC","w") ;
  write, f, "# parameter  min max mode start_point"
  for (i=1 ; i<=N_params ; i++) {
    C = Config(i) ;
    write, f,  C.parameter, C.vmin, C.vmax, C.mode, BEST.Parameters(i) ;
  } ;
  close, f ;

  /*
  P = P/sum(P) ;

  gs_nm, win+1, 3, int(ceil((N_params*1.0)/3)), landscape=1, dx=0.07, dy = 0.07, square=1 ; fma ;
  N_bins = 100 ;
  for (iParameter=1 ; iParameter<=N_params ; iParameter++) {
    plsys, iParameter ;

    X = XX(iParameter,) ; order = sort(X) ;
    X = X(order) ;
    Proba = P(order) ;

    P_plot = X_plot = array(0.,N_bins) ;
    N_obs = numberof(X)/N_bins ;

    kmin = 1 ;
    for (i=1 ; i<=N_bins ; i++) {
      p = 0 ;
      x = 0 ;
      n = 0 ;

      for (k = kmin ; (n <= N_obs) && (k < numberof(X)) ; k++) {
        p += Proba(k) ;

        if (Config(iParameter).mode =="log") {
          x += log(X(k)) ;
        } else {
          x += X(k) ;
        }
        n++ ;
        kmin = k ;
      } // k
      P_plot(i) = p/n ; // todo, je pense qu'il faut diviser par le dx, pas sure car on fait des bins de meme taille
      if (Config(iParameter).mode == "log") {
        X_plot(i) = exp(x/n) ;
      } else {
        X_plot(i) = x/n ;
      }
    } // N_bins

    plh, P_plot, X_plot, color=color ;
    limits ; limits, min(X), max(X), 0, 1.1* max(P_plot);

    if (Config(iParameter).mode == "log") {
      lxy, 1, 0 ;
    } else {
      lxy, 0, 0 ;
    }
    xytitles, Config(iParameter).parameter, "Probability", [0.02,0.02] ;

  } // iParameter

  window, win ;
  */
  best_model = BEST ;

  write, "Downloading best model :" ;
  system, "rm -rf "+running_dir+"/best" ;
  system, "mkdir -p "+running_dir+"/best" ;
  system, "rsync -Pur "+cluster+":"+BEST.dir+"/* "+running_dir+"/best" ;
  write, "Done" ;
}

//************************************************************************

func sed_fitting_routine(model) {

  extern Av, OBS, Genetic_N_explored_param, Genetic_rt ;

  obs = OBS ;

  N_obs = numberof(obs) ;
  N_freedom = N_obs - (Genetic_N_explored_param + 2) ;

  ou = where(obs.flux > 0.) ; // donnees
  ou2 = where(obs.flux == 0.) ;  // upper limit

  // Fit en log sur les points de donnees
  // en lineaire sur les limites superieures
  obs(ou).erreur = max(obs(ou).erreur,obs(ou).flux * 0.01) ; // minimum error ;

  obs(ou).erreur = obs(ou).erreur / obs(ou).flux ;
  obs(ou).flux = log(obs(ou).flux) ;


  M = open_mcfost(model.dir,"th") ;
  if (!is_void(M)) {
    if (Genetic_rt) {
      n_incl = dimsof(M.sed_rt)(3) ;
    } else {
      n_incl = dimsof(M.sed)(3) ;
    }
    n_Av = numberof(Av) ;
    chi2=array(double, n_incl, n_Av);

    // Correction modeles par A_V
    extinction = OpenASCII(get_cwd()+"/SED/extinction_law.dat",prompt=0);

    ext_V = abs(extinction.lamb - 5.47e-01)(mnx);
    kext = extinction.kpa / (1.0-extinction.albedo) ;
    kext =  kext /  kext(ext_V);
    kext = interp( kext , extinction.lamb ,  M.lamb);

    tau_V=0.4*log(10)*Av;
    correct_Av = exp(-tau_V(,-:1:numberof(M.lamb)) * kext(-:1:n_Av,));

    // Boucle sur i et Av
    for (i=1 ; i <= n_incl ; i++) {
      for (j=1 ; j <= n_Av ; j++) {
         if (Genetic_rt) {
           m = interp(M.sed_rt(,i,1,1) * correct_Av(j,), M.lamb, obs.lamb) ;
         } else {
           m = interp(M.sed(,i,1,1) * correct_Av(j,), M.lamb, obs.lamb) ;
         }
        // Fit en log sur les points de donnees
        m(ou) = log(max(m(ou),1e-30)) ;
        C2 = ( (( m - obs.flux)/(obs.erreur))^2)(sum) ;
        chi2(i,j) = C2 ;
       }
     }


    best = chi2(*)(mnx);   best = indexof(chi2,best);
    f = open(model.dir+"/best","w") ;
    write, f, best(1), Av(best(2)) ;
    close, f ;

    // Proba marginale sur i et Av avec chi2 reduit
    if (N_freedom > 1) {
      P = sum(exp(-chi2/(2*N_freedom))) ;
    } else {
      P = -1 ;
    }

  } else {
    chi2 = 1e9 ; P = 0 ;
    f = open(model.dir+"/best","w") ;
    write, f,  0 , 0,  " #there was a pb with this model" ;
    close, f ;
   }

  //return chi2 ;
  return [min(chi2),P] ;
 }

//************************************************************************

func image_fitting_routine(model)  {
   extern img, mask, noise, psf ;

   S = sum(img)
   img = img / S ;
   noise = noise / S ;


   fft_img = fft(img) ;

   model_img = model.image_rt ;

   dims = dimsof(model_img) ;
   n_incl = dims(4) ;
   center = 0.5 * (dims(2:3)+1);

   chi2=array(double, n_incl) ;

   for (i=1 ; i <= n_incl ; i++) {
   //for (i=10 ; i <= 10 ; i++) {
     img_tmp = model_img(,,i,1,1) ;
     img_tmp = convol2df(img_tmp,psf) ;

     img_tmp = img_tmp/ sum(img_tmp) ;

     //window, 1, style="boxed.gs" ;
     //pli, img ;

     //window, 2, style="boxed.gs" ;
     //pli, img_tmp ;

     // compute correlation function
     CX = roll(abs( fft( (fft(img_tmp)) * conj(fft_img), -1) )) ;
     //cfitsWrite, "CX.fits", CX ;

     //window, 3, style="boxed.gs" ;
     //pli, CX ;

     // find peak position
     peak = indexof(CX,CX(*)(mxx)) ;
     offset = center - peak ;

     // apply offset accordingly
     img_tmp = fft_fine_shift(img_tmp, offset) ;

     //cfitsWrite, "mod.fits", img_tmp ;

     //window, 4, style="boxed.gs" ;
     //pli, img_tmp ;

     diff = ((img_tmp - img)/noise)^2 ;

     //window, 5 ; style="boxed.gs" ;
     //pli, diff ;
     //cfitsWrite, "diff.fits", diff ;

     chi2(i) = sum(diff) ;
   }

   return chi2 ;
}


//************************************************************************

func read_genetic_chi2(model) {

  chi2=float() ;
  f = open(model.dir+"/chi2.txt","r") ;
  read, f, chi2 ;
  close, f ;

  return [chi2,exp(-0.5*chi2)] ;
}


write, "MCFOST-genetic loaded" ;

//************************************************************************

func get_GA(grid_name,running_dir=,cluster=)  {

  extern default_cluster ;

  if (is_void(cluster)) cluster=default_cluster
  if (is_void(running_dir)) running_dir = "~/Simus/Genetic/"+grid_name ;

  system, "mkdir -p "+running_dir ;
  system, "rsync -Pur  "+cluster+":"+running_dir+"/config.genetic "+running_dir+"/" ;
  system, "rsync -Pur  \""+cluster+":"+running_dir+"/*.bin\" "+running_dir+"/" ;
}
