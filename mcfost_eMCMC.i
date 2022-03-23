require, "~/yorick/init/mcfost_genetic.i"
eMCMC_file = "~/yorick/init/mcfost_eMCMC.i" ;
eMCMC_N_explored_param  = (eMCMC_N_explored_param ? eMCMC_N_explored_param : 1) ; // will be updated

struct eMCMCModel {
  int iteration ;
  int id
  float chi2 ;
  float Proba ;
  float Parameters(eMCMC_N_explored_param) ;
  string dir ;
  double Z ;
  int out ;
  int n_accepted ;
} ;

/* // Same as genetic --> move to mcfost_struct ??
struct ClusterConfig {
  int n_nodes ;
  int omp_num_threads ;
  int walltime ;
  string running_dir ;
  string grid_name ;
  int n_options ;
  string options(10) ;
} ;
*/

struct eMCMCRun {
  ClusterConfig FostiConfig ;
  float A ;
  int iteration_max ;
  int iteration_reached ;
  int n_models ;
  string simu_config_file ;
  float seed ;
}

//************************************************************************

  func mcfost_eMCMC(grid_name,parameter_file,fitting_routine, simu_config_file=, n_nodes=,omp_num_threads=,walltime=,root_dir=, options=,n_models=,iteration_max=,no_mcfost_compute=,progressive_plot=,waiting_time=,rt=) {
/* DOCUMENT mcfost_genetic(grid_name,parameter_file,fitting_routine,
   simu_config_file=,n_nodes=,omp_num_threads=,walltime=,root_dir=,
   options=,n_models=,iteration_max=)

   Default values :
   - simu_config_file : root_dir+"/"+grid_name+"/config.eMCMC"
   - n_nodes : 20
   - omp_num_threads : 8
   - walltime (per iteration) : 24
   - waiting_time (to check if generation if finished) : 30s
   - root_dir : ~/Simus/eMCMC/
   - options : none
   - n_models : 100
   - iteration_max : 100
   - no_mcfost_compute = 0 : skip the mcfost calculations, for fitting routines that do not need mcfost
   - progressive_plot = 1

   SEE ALSO:
 */

  extern eMCMC_N_explored_param, eMCMC_file ;

  if (is_void(rt)) rt=0 ;
  Genetic_rt = rt ;

  // eMCMC config
  if (is_void(iteration_max)) iteration_max = 100 ;
  if (is_void(n_models)) n_models=100 ;
  if (is_void(no_mcfost_compute)) no_mcfost_compute=0 ;
  if (is_void(progressive_plot))  progressive_plot=1 ;


  // Cluster configuration
  if (is_void(omp_num_threads)) omp_num_threads = 8 ;
  if (is_void(walltime)) walltime = 24 ;
  if (is_void(waiting_time)) waiting_time = 30 ;
  if (is_void(n_nodes)) n_nodes = 20 ;
  if (is_void(root_dir)) {
    running_dir = "~/Simus/eMCMC/"+grid_name ;
  } else {
    running_dir = root_dir+"/"+grid_name ;
  }
  if (is_void(simu_config_file)) simu_config_file = running_dir+"/config.eMCMC" ;
  n_options = numberof(options) ;


  // Read configuration file
  SimuConfig = OpenASCII(simu_config_file,prompt=0) ;
  eMCMC_N_explored_param = numberof(SimuConfig) ;

  // Recharge la structure
  include, eMCMC_file, 1 ;

  if (n_models <= 2*eMCMC_N_explored_param) {
    write, "The number of walkers needs to be more than twice the" ;
    write, "dimension of your parameter space" ;
    write, "Exiting" ;
    return [] ;
  }

  // eMCMC structure
  eMCMC = array(eMCMCRun) ;
  eMCMC.iteration_max = iteration_max ;
  eMCMC.iteration_reached = 0 ;
  eMCMC.n_models = n_models ;
  eMCMC.simu_config_file = simu_config_file ;
  A = 2.0 ; eMCMC.A = A ;

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
  eMCMC.FostiConfig = FostiConfig ;


  // Read refence parameter file
  V = 2.17 ;
  ref = read_params(parameter_file) ;
  if (abs(V - ref.simu.version) > 1.e-4) {
    write, "Version ", V, "needed. Exiting !!!";
    return;
  }

  // Mise a jour des structures MCFOST
  include, struct_file ;
  pldefault, style="landscape.gs", marks=0, width=1, palette="heat.gp", legends=0, dpi=75;

  walkers = [] ; // TODO : le modele de reference

  // Routine de fit
  if (typeof(fitting_routine) == "string") {
    post_process = fitting_routine ;      // commande eexcute sur chaque noeud
    fitting_routine = read_genetic_chi2 ; // on ne lit que le chi2
  } else {
    post_process = [] ;
  }

  for (iteration=1 ; iteration <= iteration_max ; iteration++) {

    write, "------------------------------------------------------" ;
    write, "Iteration", iteration ;

    // Creation du dossier pour la iteration
    dir=running_dir+"/"+"it"+swrite(format="%4i",iteration) ;
    for (j=1 ; j<=3 ; j++) dir = streplace(dir,strgrep(" ",dir),"0");

    mkdirp, dir ;

    // Create the new population
    write, "Creating the new walkers population ..." ;
    if (iteration == 1) {
      write, "Initializing parameter file for first iteration" ;
      eMCMC_seed = randomize() ;  eMCMC.seed = eMCMC_seed ;
      write, "Random seed =", eMCMC_seed ;
    }
    proposed_walkers = propose_walkers(SimuConfig, walkers, A, iteration, n_models, dir) ;

    IN = where(proposed_walkers.out==0) ;

    if (!no_mcfost_compute) {
      // Create the corresponding parameter files
      write, "Creating the parameter files ..." ;
      make_parameter_files, SimuConfig, ref, proposed_walkers(IN), iteration ;

      // Run the models --> system call
      write, "Computing the MCFOST models ..." ;
      compute_mcfost_models, iteration, FostiConfig, options, post_process, dir ;
    }

    // Compute the chi2
    write, "computing fits ..."
    compute_fits, fitting_routine, proposed_walkers, iteration, dir ;

    // Select the models that will be used to generate the next iteration of models
    write, "Selecting walkers ..."
    walkers = select_walkers(walkers, proposed_walkers) ;

    // enregistrement de la pop en fichier binaire yorick
    f = createb(dir+".bin") ;
    save, f, proposed_walkers, walkers ;
    close, f ;

    // Saving the state of the eMCMC
    eMCMC.iteration_reached = iteration ;
    f = createb(running_dir+"/eMCMC.bin") ;
    save, f, eMCMC ;
    close, f ;
  }

  write, "eMCMC DONE" ;
  return [] ;
}

//************************************************************************

func propose_walkers(SimuConfig,walkers,A,iteration,n_models, dir) {
  /*
    Creates a table of models to compute for the iteration specified.

    If this is the iteration 0, then the parameters are sampled around
    the specified model.
    Otherwise, this method uses results from previous iterations to determine
    which models to run.
  */

  extern eMCMC_N_explored_param ;

  proposed_walkers = array(eMCMCModel,n_models) ;

  if (iteration == 1) {
    proposed_walkers.iteration() = 1 ;
    proposed_walkers.id() = indgen(n_models) ;

    // Creating the array of parameter values
    for (p=1 ; p<=eMCMC_N_explored_param ; p++) {
      CP = SimuConfig(p) ;

      // Gaussian ball around the initial parameter file
      r = Gauss_random(n_models) ;

      if (CP.mode == "linear") {
        sigma = 0.1 * (CP.vmax - CP.vmin)  ;
        proposed_walkers.Parameters(p,) = CP.start_point + r * sigma ;
      } else {
        sigma =  0.1 *(log10(CP.vmax) - log10(CP.vmin)) ;
        proposed_walkers.Parameters(p,) =  10^(log10(CP.start_point) + r * sigma)  ;
      }

      proposed_walkers(:).out = 0 ;
      for (i = 1 ; i<= n_models ; i++) {
        // Make sure we don't go above the limits
        if (proposed_walkers(i).Parameters(p) > CP.vmax) {
          //proposed_walkers(i).Parameters(p) =  CP.vmax ;
          proposed_walkers(i).out =  1 ;
        }
        if (proposed_walkers(i).Parameters(p) < CP.vmin) {
          //proposed_walkers(i).Parameters(p) =  CP.vmin ;
          proposed_walkers(i).out =  1 ;
        }

      } // i
    } // p

  } else { // iteration > 1

    // Creation new population
    for (i = 1 ; i<= n_models ; i++) {

      proposed_walkers(i).iteration = iteration ;
      proposed_walkers(i).id = i ;
      proposed_walkers(i).out = 0 ;


      // Select a random model /= i
      do {
        j = int_random(n_models) ;
      } while (j == i) ;

      // Stretch move
      Z = ( (A - 1.) * random() + 1)^2. / A ;
      proposed_walkers(i).Z = Z ;
      for (p=1 ; p<=eMCMC_N_explored_param ; p++) {
        CP = SimuConfig(p) ;
        if (CP.mode == "linear") {
          proposed_walkers(i).Parameters(p) = walkers(j).Parameters(p) + Z * (walkers(i).Parameters(p) -  walkers(j).Parameters(p)) ;
        } else {
          proposed_walkers(i).Parameters(p) = 10^( log10(walkers(j).Parameters(p)) + Z * (log10(walkers(i).Parameters(p)) -  log10(walkers(j).Parameters(p))) ) ;
        }

        // Make sure we don't go above the limits
        if (proposed_walkers(i).Parameters(p) > CP.vmax) {
          //proposed_walkers(i).Parameters(p) =  CP.vmax ;
          proposed_walkers(i).out = 1 ;
        }
        if (proposed_walkers(i).Parameters(p) < CP.vmin) {
          //proposed_walkers(i).Parameters(p) =  CP.vmin ;
          proposed_walkers(i).out = 1 ;
        }
      } // p
    } // i

  } // test 1er iteration

  // Dossier dans lequel va etre la simu
  n_chiffres = floor(log10(n_models) + 1);
  format = "%"+swrite(format="%i",int(n_chiffres))+"i";
  for (i=1 ; i<=n_models ; i++) {
    number=swrite(format=format,i);
    // Ajout des 0
    for (j=1 ; j<=n_chiffres -1 ; j++) {
      number = streplace(number,strgrep(" ",number),"0");
    }
    proposed_walkers(i).dir = dir+"/"+number ;
  }

  return proposed_walkers ;
}


//************************************************************************

func select_walkers(walkers,proposed_walkers) {

  if (is_void(walkers)) {
    return proposed_walkers ;
  } else {
    N_dim = eMCMC_N_explored_param ;
    R = random(numberof(proposed_walkers)) ;
    Z = proposed_walkers.Z ;

    // Decide whether or not the proposals should be accepted.
    accept = (R < Z^(N_dim-1) * (proposed_walkers.Proba / (walkers.Proba +1e-300))) ;
    ou = where(accept) ;

    accepted_walkers = walkers ;

    if (numberof(ou)) {
      accepted_walkers(ou) = proposed_walkers(ou) ;
      accepted_walkers(ou).n_accepted = walkers(ou).n_accepted + 1 ;
    }

    write, "Acceptance fraction =", (1.0*numberof(ou))/numberof(walkers) ;

    return accepted_walkers ;
  }
}

//************************************************************************

func acor(X,&mean,&sigma,&tau,warning=) {
/* DOCUMENT acor(X,&mean,&sigma,&tau)
     Acor estimates an error bar and the autocorrelation time of
     a time series that comes from a Markov chain.

     Adapted from C++ code : Jonathan Goodman, goodman@cims.nyu.edu, http://www.cims.nyu.edu/faculty/goodman

     X :  the time series to be analized.

     Return :
      - the mean of X, or possibly a slightly shorter sequence.
      - an estimate of the standard devation of the sample mean.
      - an estimate of the autocorrelation time.

   SEE ALSO:
 */

  if (is_void(warning)) warning=1 ;

  taumax = 2 ;
  winmult = 5 ;
  maxlag = taumax*winmult ;
  minfac = 5 ;

  // compute the mean and subtract it
  L = numberof(X) ;

  if ( (L < minfac * maxlag) ) {
    if (warning)
      write, "The autocorrelation time is too long relative to the variance" , L, minfac, maxlag ;
    return ;
  }

  mean = X(avg) ;
  X -= mean ;
  C = array(float,maxlag) ;

  // compute the autocovariance function
  iMax = L - maxlag ;
  C0 = sum(X^2)/iMax ;
  for (t=1 ; t<=maxlag ; t++) {
    C(t) = sum(X(1:iMax-1)*X(1+t:iMax-1+t))/iMax ;
  }

  // Diffusion coefficient : sum of autovariance
  D = C0 + 2*sum(C) ;
  if (D < 0) {
    write, "Warning : diffusion coeff is negative" ;
    D = - D ;
  }

  sigma = sqrt(D/L) ;
  tau = D/C0 ; // autocorrelation time

  if (tau < taumax) {
    return ;
  } else {
    Y = array(float, L/2) ;
    T2 = indgen(L/2)*2 ;
    X = X(T2-1) + X(T2) ;
    newMean = float() ;
    acor, X, newMean, sigma, tau, warning=0;
    D = 0.25 * sigma^2 * L ;
    tau = D/C0 ;
    sigma = sqrt(D/L) ;
  }

}


//************************************************************************

func eMCMC_stat(grid_name,running_dir=,color=,win=,iteration_max=,nx=,exclude_N_Tac=,nBars=,save_fits=) {
/* DOCUMENT eMCMC_stat(grid_name,running_dir=,color=,win=,iteration_max=,nx=,exclude_N_Tac=,nBars=,save_fits=) {

   SEE ALSO:
 */


  if (is_void(running_dir)) running_dir = "~/Simus/eMCMC/"+grid_name ;

  if (is_void(color)) color="black" ;
  if (is_void(nBars)) nBars=20 ;
  if (is_void(win)) win = 0;
  if (is_void(save_fits)) save_fits = 0;
  if (is_void(nx)) nx = 3 ;
  if (is_void(exclude_N_Tac)) exclude_N_Tac=3 ;
  pltitle_height = 12 ;

  Config = OpenASCII(running_dir+"/config.eMCMC",prompt=0) ;
  N_params = numberof(Config) ; // TODO : parametres virtuels : i et av par exemple
  window, win ; fma ; gs_nm, win, nx, int(ceil((N_params*1.0)/nx)), landscape=1, dx=0.02, dy = 0.07, square=0 ;

  f = openb(running_dir+"/eMCMC.bin") ;
  restore, f, eMCMC ;
  close, f ;

  if (is_void(iteration_max))   iteration_max = eMCMC.iteration_reached ;
  write, iteration_max, "iterations found  --> ", iteration_max * eMCMC.n_models, "models" ;


  P = [] ;
  Proba = [] ;
  chi2 = [] ;

  minchi2 = 1e300 ;
  for (iteration=1 ; iteration <= iteration_max ; iteration++) {
  //  i = 901
  //i = 903
  //for (iteration=i ; iteration <= i+15 ; iteration++) {
    pop = [] ;
    dir=running_dir+"/"+"it"+swrite(format="%4i",iteration) ;
    for (j=1 ; j<=3 ; j++) dir = streplace(dir,strgrep(" ",dir),"0");

    //write, "------------------------------------" ;
    write, "reading "+dir+".bin" ;
    f = openb(dir+".bin") ;
    restore, f, walkers, proposed_walkers ;
    close, f ;

    // TODO : avoid grow
    grow, chi2, [walkers.chi2] ;
    grow, Proba, [walkers.Proba] ;
    grow, P, [walkers.Parameters] ;  // TODO : pour test : walkers normalement
  } // iteration

  write, "Mean acceptance fraction:", avg((walkers.n_accepted*1.0)/walkers.iteration) ;
  write, " " ;

  best = Proba(*)(mxx) ;
  best = indexof(Proba,best);

  write, "Model with highest proba in eMCMC is model ", best(1), "in iteration ", best(2)
  write, " " ;


  tau_max = 0.
  for (iParameter=1 ; iParameter<=N_params ; iParameter++) {
    acor, P(iParameter,avg,) ,mean,sigma,tau ; // moyenne de la chaine a chaque pas de temps
    if (!is_void(tau)) {
      write, "Auto-correlation time for", Config(iParameter).parameter, tau ;
      if (tau > tau_max) tau_max = tau ;
    }
  }

  write, "Maximum auto-correlation time = ", tau_max ;
  iStart = int(ceil(tau_max * exclude_N_Tac)) ;
  write, "Excluding", iStart -1, "first chains for bayesian inference" ;
  write, " " ;


  if (save_fits) {
    write, "Exporting all the parameter values for the selected models :"
    for (iParameter=1 ; iParameter<=N_params ; iParameter++) {
      xname = Config(iParameter).parameter ;
      filename = running_dir+"/"+xname+".fits.gz" ;
      if (strpart(filename,1:1) == "~") filename=get_env("HOME")+strpart(filename,2:0);
      write, "writing "+filename ;
      cfitsWrite, filename, P(iParameter,,) ;
    }
  }

  P = P(,*,iStart:0) ;

  for (iParameter=1 ; iParameter<=N_params ; iParameter++) {
    plsys, iParameter ;
    if (Config(iParameter).mode == "log") {
      lxy, 1, 0 ;
      logHisto=1 ;
    } else {
      lxy, 0, 0 ;
      logHisto=0 ;
    }
    xname = Config(iParameter).parameter ;
    for (i=1 ; i<=5 ; i++) xname = streplace(xname,strgrep("_", xname), " ") ;
    if (iParameter%nx == 1) {
      xytitles, xname, "Probability density", [0.0,0.02] ;
    } else {
      xytitles, xname, "", [0.02,0.02] ;
    }

    yocoPlotHistogramme, P(iParameter,), nbBars = nBars, onlyBars=1, color=color, logHisto=logHisto ;
    if (logHisto) lxy, 1, 0 ;
    limits ; l = limits() ;
    limits, , , 0, 1.1*l(4) ;
  }

  window, win+1 ; fma ;

  //N_params = 2 ;
  gs_nm, win+1, N_params, N_params, dx = 0.01, dy = 0.01, square=1 ;

  for (iParameter=1 ; iParameter<=N_params ; iParameter++) {
    for (iParameter2=iParameter+1 ; iParameter2<=N_params ; iParameter2++) {
      if (iParameter2 == iParameter+1) {
        xname1 = Config(iParameter).parameter ;
        xname2 = Config(iParameter2).parameter ;
      } else {
        xname1 = "" ;
        xname2 = "" ;
      }

      plsys, iParameter2 + N_params * (iParameter-1) ;
      plp, P(iParameter,*), P(iParameter2,*), symbol = 6, size=0.05 ;


      if (Config(iParameter).mode == "log") {
        logy = 1 ;
      } else {
        logy = 0 ;
      }
      if (Config(iParameter2).mode == "log") {
        logx = 1 ;
      } else {
        logx = 0 ;
      }
      logxy, logx, logy ;
      xytitles, xname2, xname1, [0.04,0.04]

    }
  }
}

//************************************************************************

func get_eMCMC(grid_name,running_dir=)  {

  if (is_void(running_dir)) running_dir = "~/Simus/eMCMC/"+grid_name ;

  system, "mkdir -p "+running_dir ;
  system, "rsync -Pur  fosti:"+running_dir+"/config.eMCMC "+running_dir+"/" ;
  system, "rsync -Pur  \"fosti:"+running_dir+"/*.bin\" "+running_dir+"/" ;
}


write, "MCFOST-eMCMC loaded" ;

