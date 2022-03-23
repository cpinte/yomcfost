//Include other files;
//#include "plot_3d.i"

MCFOST_version = 3.0 ;

lambda_min=0.1;
lambda_max=3000.0;

// load dimension parameters used
nrad  = (nrad ? nrad : 60);
nz    = (nz   ? nz   : 20);
n_az  = (n_az   ? n_az   : 1);
nrad_in  = (nrad_in ? nrad_in : 10);
nthet  = (nthet ? nthet : 10);
RT_ntheta  = (RT_ntheta ? RT_ntheta : 1);
nphi  =  (nphi  ? nphi  : 1);
nlamb = (nlamb? nlamb: 100);
nlamb1 = (nlamb1? nlamb1: 51) ;
nx  =  (nx  ? nx  : 51);
ny  =  (ny  ? ny  : 51);
n_type  =  (n_type  ? n_type  : 7);
n_T = (n_T  ? n_T  : 1);
n_T_grain = (n_T_grain  ? n_T_grain  : 1);
data_dir = (data_dir ? data_dir :"undefined");
n_zones = (n_zones ? n_zones : 1) ;
n_especes = (n_especes ? n_especes : 7) ;
n_components = (n_components ? n_components : 3) ;
n_stars = (n_stars ? n_stars : 1) ;
nmol = (nmol ? nmol : 1) ;
n_lines_ray_tracing = (n_lines_ray_tracing ? n_lines_ray_tracing : 3);
n_Trans = (n_Trans ? n_Trans : 3) ;
n_speed = (n_speed ? n_speed : 50) ;

if (is_void(vect_n_phi)) vect_n_phi = indgen(nphi);
nphiV = numberof(vect_n_phi);

if (is_void(vect_n_sed))
  vect_n_sed = indgen(9); // Vector of num flux in sed_file
nsed = numberof(vect_n_sed);

if (is_void(vect_n_sed1))
  vect_n_sed1 = [1]; // Vector of num flux in sed_file
nsed1 = numberof(vect_n_sed1);


//struct_file = get_cwd()+"/init/mcfost_struct.i"; //put the complete path if this file is not in your yorick directory
struct_file = "~/yorick/init/mcfost_struct.i"; //put the complete path if this file is not in your yorick directory
sed1_file = ".sed_th.fits.gz";
sed2_file = "sed_mc.fits.gz";
sed_rt_file = "sed_rt.fits.gz";
image_file = "MC.fits.gz"
image_rt_file = "RT.fits.gz"
temperature_file = "Temperature.fits.gz";
spectre_file = "lines.fits.gz";

flux_adjust =  1.; // Correction eventuelle distance
//lamb_adjust =  1e-6;    // um to m
lamb_adjust =  1;    // um to m

// STRUCTURE DEFINITION
struct McfostPhoton {
  int n1;
  int n2;
  int nlambda;
  int nimage
}
struct McfostWavelength {
  string file;
  int nlamb;
  int nlamb1;
  float lambda_min
  float lambda_max
}
struct McfostDust {
  float sublim;
  string lstrat;
  string lmigration ;
  int strat_type ;
  float strat;
  float a_strat;
  float gas_dust ;
}
struct McfostDustPop {
  //float dens;
  float porosity;
  float amin;
  float amax;
  float aexp;
  string file(n_components) ;
  int type_grain;
  int n_grains;
  int n_components ;
  int mixing_rule ;
  float mass_fraction;
  float component_volume_fraction(n_components);
  string type;
  float dhs_maxf ;
}

struct McfostGrid {
  double R_grid(nrad);
  double H_grid(nrad);
}

struct McfostGridp {
  int nrad;
  int nz;
  int n_az;
  int nrad_in;
  int geometry;
}

struct McfostGridT {
  int n_T;
  int n_T_grain;
  float T_min;
  float T_max;
}

struct McfostMaps {
  int nthet;
  int nphi;
  int nx;
  int ny;
  int n_type;
  float zoom;
  int inc;
  int delta;
  float dist;
  float PA;
  float size_neb;
  float map_size ;
  float angle_interest;
  string lonly_bin_interest;
  float RT_imin ;
  float RT_imax ;
  float RT_az_min ;
  float RT_az_max ;
  int RT_ntheta ;
  int RT_naz ;
  string lRT_i_centered ;
  float RT_incl(RT_ntheta) ;
}
struct McfostZone {
  double mass;
  double gas_to_dust;
  double Ho;
  double Ro;
  double vert_exp ;
  double rout;
  double Rc;
  double rin;
  double edge;
  double beta;
  double dens;
  double gamma_exp;
  int geometry;
  double size_neb;
  int n_especes;
}

struct McfostCavity{
  string lcavity ;
  float H0 ;
  float R0 ;
  float beta ;
}


struct McfostMol{
  float vmax;
  float vturb;
  int n_speed;
  string lpop, lpop_precise, lLTE;
  float profile_width;
  string molecule_file;
  int level_max;
  string cst_abundance;
  float abundance;
  string abundance_file;
  string ray_tracing;
  int n_lines_ray_tracing;
  int transition_numbers(n_lines_ray_tracing);
}

struct McfostStars {
  float T;
  float R;
  float Mass;
  string is_bb;
  string spectre;
  float x, y, z;
  float fUV, slope_fUV ;
}

struct McfostSimu {
  float version;
  string lcheckpoint;
  int checkpoint_time;
  string ltemp, lsed, lsed_complete;
  string lem_disk_image, lsepar, lsepar_pola;
  float tau_seuil, wl_seuil;
  int scatt_method, aniso_method;
  string lsym_ima, lsym_centrale, lsym_axiale;
  string ldust_sublimation;
  float correct_Rsub ;
  string lchauff_int;
  string lhydro_eq ;
  float viscosity;
  int n_zones;
  int n_especes_max ;
  int n_components_max ;
  int n_stars;
  int nmol ;
}

struct McfostParams {
  string dir;
  McfostPhoton phot;
  McfostWavelength wave;
  McfostDust dust;
  McfostDustPop dust_pop(n_zones,n_especes);
  McfostGridp grid;
  McfostGridT gridT;
  McfostMaps map;
  McfostZone zones(n_zones);
  McfostCavity cavity ;
  McfostStars stars(n_stars);
  McfostSimu simu
  McfostMol mol(nmol)
}

struct McfostSED {
  McfostParams P;
  McfostGrid grid;
  float T(nrad, nz, n_az, n_T_grain);
  float lamb1(nlamb1);
  float lamb(nlamb);
  float sed1(nlamb1, nthet , nphi, nsed1);
  float sed(nlamb, nthet , nphi, nsed);
  float sed_rt(nlamb, RT_ntheta , nphi, n_type);
  string dir;
  float h_hydro(nrad)
}

struct McfostImage {
  McfostParams P;
  McfostGrid grid;
  float lamb;
  float image(nx, ny, nthet , nphi, n_type);
  float image_rt(nx, ny, RT_ntheta , nphi, n_type);
  string dir;
}

struct McfostSpectre {
  McfostParams P;
  McfostGrid grid;
  float Fline(nx, ny, 2*n_speed+1, n_lines_ray_tracing, RT_ntheta) ;
  float Fcont(nx, ny, n_lines_ray_tracing, RT_ntheta) ;
  int ifreq(n_lines_ray_tracing) ;
  float freq(n_lines_ray_tracing) ;
  float velocity(2*n_speed+1) ;
  string dir;
  string mol_name ;
  }


//***************************************************

n_params_tot = 53;
_n_valeurs_max = 40 ;

struct grid {
  string name ;

  // MCFOST
  int mcfost_version ;

  // Clusters
  int n_clusters ;
  string clusters_name(n_clusters) ;

  // Modeles
  int n_models ;
  McfostParams params(n_models) ;
  string model_name(n_models) ;
}


struct MCFOST_grid {
  int dims(n_params_tot);
  string names(n_params_tot);

  float values(n_params_tot,_n_valeurs_max) ;
  McfostStars stars(_n_valeurs_max);

  McfostParams ref ;
}

//***************************************************
