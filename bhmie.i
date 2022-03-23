func mie(amin,amax,aexp,n_grains,wl,optical_indices_file,nang,&kappa,&albedo,&S11,&S12,&gsca2,&q_ext,&q_sca,&g_sca,&tab_a,&nbre_grains,porosite=,P2=,color=)
{
/* DOCUMENT : mie(amin,amax,aexp,n_grains,wl,optical_indices_file,nang,&kappa,&albedo,&S11,&S12,&gsca2)
   Calcule les proppriétés optiques d'une distribution en tailles de grains
   arguments : amin,amax,aexp,n_grains,wl,mgrain,rho1g,nang,&kappa,&albedo,&S11,&gsca
   Unités :
      - amin, amax et [wl] en microns
      - aexp should be > 0
      - mgrains : indice complexe (partie imaginaire positive)
      - rho1g en g.cm^-3
      - nang  nbre angle de phases , 90 correspond a un pas de 1°
   Output :
      - kappa en cm²/g (g de poussiere : diviser par 100 pour avoir en g de gaz)
      - S11 normalise a 1 en integrant sur l'angle de diffusion (pas l'azimuth)
      12/2005
      SEE ALSO:
      bhmie, plot_phase_function, hg
 */

 if (is_void(porosite)) porosite=0 ;
 if (is_void(P2)) P2=0 ;

  dim_lambda = dimsof(wl);
  if (dim_lambda(1) == 0) {
    n_lambda=1;
  } else {
    n_lambda=dim_lambda(2);
  }

  S11 = S12 = array(double, 2*nang, n_lambda);
  kappa = albedo = gsca2 = array(double,n_lambda);
  tab_a = nbre_grains = array(double, n_grains);
  q_ext = q_sca = q_back =g_sca = array(double, n_grains, n_lambda);


  // Taille des grains
  exp_grains =  exp((1.0/double(n_grains)) * log(amax/amin));
  tab_a(1) = amin*sqrt(exp_grains);
  for (i=2 ; i<=n_grains ; i++) {
    tab_a(i) = tab_a(i-1) * exp_grains ;
  }

  // Distribution en taille
  nbre_grains() = tab_a()^(1.0-double(aexp));
  nbre_grains() = nbre_grains() / nbre_grains(sum);

  // Lecture des indices optiques et de la densite
  //ind = OpenASCII(optical_indices_file, prompt=0, valtype=double);
  f=open(optical_indices_file) ;
  tmp=ReadASCII(f,prompt=0,valtype=double,noheader=1);
  rho1g = tmp.X1;
  ind=ReadASCII(f,prompt=0,valtype=double);
  close, f;

  // on remet dans l'ordre
  order = sort(ind.X1);
  X1 = ind.X1(order); X2 = ind.X2(order); X3 = ind.X3(order);

  //  X2 = X2 + 1 ; write, "TMP for files directly from Draine's website"
  //min(X2)

  // porosite
  if (porosite) {
    // Cst matrice de matiere
    m = X2 + 1i * X3 ;
    e1 = m^2;
    f1 = 1.0 - porosite;

    // Cst inclusions de vide
    e2 = 1.0;
    f2 = porosite;

    // 80% de vide
    // arithmetique ==> Birchak
    //O.X2 = 0.2 * O.X2 + 0.8 * 1.0 ;
    //O.X3 = 0.2 * O.X3 + 0.8 * 0.0 ;

    // Maxwell-Garnet : BUG
    //coeff = (e0 - em)/(e0 + 2*em)  ;
    //coeff = -coeff // ????
    //eavg = em   * (1. * (3*f *coeff ) / (1-f*coeff));


    // Bruggeman
    //b = (1-3*f1)*e1 + (1-3*f2)*e2 ;
    b = 2*f1*e1 -f1*e2 + 2*f2*e2 -f2*e1 ;
    c = e1 * e2;
    a = -2.;
    delta = b^2 - 4*a*c ;
    eavg = (-b - sqrt(delta)) / (2.*a) ; // -sqrt(delta) sinon partie imaginaire de m < 0

    //    print, "verif", max(abs(f1 * (e1 - eavg)/(e1 + 2*eavg) + f2 * (e2 - eavg)/(e2 + 2*eavg))) ;


    // Mathis & Whieffen
    if (P2) eavg = f1 * e1 + f2 * e2;


    mavg = sqrt(eavg) ;

    X2 = mavg.re; X3 = mavg.im;
    rho1g = rho1g * (1.-porosite) ;
  }


  // on doit extrapoler ?
  if (max(wl) > X1(0)) { // A la main car interp n'extrapole pas
    frac = (log(max(wl))-log(X1(-1)))/(log(X1(0))-log(X1(-1)));
    grow, X2, exp(frac * log(X2(0)) + (1.-frac)*log(X2(-1)));
    grow, X3, exp(frac * log(X3(0)) + (1.-frac)*log(X3(-1)));
    grow, X1, wl(0);
    }
  // interpolation en log-log
  mgrain = exp(interp(log(X2), log(X1), log(wl))) + 1i * exp(interp(log(X3), log(X1), log(wl)));

  /*
  window, 20 ; lxy, 1, 1 ;
  plg, mgrain.re, wl, color=color ;

  window, 21 ; lxy, 1, 1 ;
  plg, mgrain.im, wl,color=color ;
  */

   // Pour normalisation fonction de phase
  dang=.5*pi/double(nang);
  tab_cos = cos((double(indgen(2*nang))-0.5)*dang);
  tab_sin = sin((double(indgen(2*nang))-0.5)*dang);


  for (l=1 ; l<= n_lambda ; l++) {
    // Calcul theorie de Mie et moyennage fonction de phase

    for (i=1 ; i<=n_grains ; i++) {
      x= 2*pi*tab_a(i)/wl(l) ;
      v = bhmie(x,mgrain(l),nang,s11,s12,qext,qsca,qback,gsca) ;
      q_ext(i,l) = qext ;
      q_sca(i,l) = qsca ;
      g_sca(i,l) = gsca ;
      S11(,l) = S11(,l) + s11 * nbre_grains(i);
      S12(,l) = S12(,l) + s12 * nbre_grains(i);
    }


    norme=sum(S11(,l)*tab_sin);
    //write, "PAS DE NORME" ; norme = 1.0 ;
    S11(:,l) = S11(:,l)/norme;
    S12(:,l) = S12(:,l)/norme;


    // G fonction de phase moyenne
    // Verifie egal au G moyen
    /*  norme = sum(S11() * tab_sin());
        gsca1 = sum(S11() * tab_sin() * tab_cos());
        gsca1 = gsca1/norme;
    */

    // G moyen
    sum0 = sum(tab_a()^2*q_sca(,l) * nbre_grains());
    norme = sum0;
    gsca2(l) = sum(tab_a()^2*q_sca(,l) * nbre_grains() * g_sca(,l));
    gsca2(l) = gsca2(l)/norme;

    // Opacite massique et albedo
    sum1= sum(tab_a()^2*q_ext(,l) * nbre_grains());
    sum2 =sum(tab_a()^3 * nbre_grains()) * 1.e-4 * (4.*rho1g / 3.) ;

    //sum2 ;

    // en cm².g^-1 si rho1g en g.cm^-3 et a en microns
    kappa(l) =  sum1/sum2 ;
    albedo(l) = sum0/sum1;


  }

  // Cas 1 lambda
  if (n_lambda==1) {
    S11 = S11(,1);
    S12 = S12(,1);
    kappa=kappa(1);
    albedo=albedo(1);
    gsca2=gsca2(1);
  }

  return ;
}



//***********************************************************************

func bhmie(xgrain,mgrain,nang,&s11,&s12,&qext,&qsca,&qback,&gsca)
/* DOCUMENT
   arguments : xgrain,mgrain,nang,&s11,&qext,&qsca,&qback,&gsca
   SEE ALSO:
   mie, plot_phse_function, hg
 */
{
  amu = PI = pi0 = pi1 = tau = array(double, nang) ;
/* DOCUMENT
     Routine de Mie Bohren & Hoffman
     traduction de la routine f90 de MCFOST
     12/2005
   SEE ALSO:
   mie
 */
  s1 = s2 = array(complex, 2*nang);

  y=xgrain*mgrain ;
  ymod=abs(y) ;
  /*   series expansion terminated after nstop terms
       logarithmic derivatives calculated from nmx on down  */
  nstop=xgrain+4.*xgrain^0.3333+2. ;
  nmx=max(nstop,ymod)+15 ;

  D = array(complex,int(nmx));

  // require nang.ge.1 in order to calculate scattering intensities
  dang=0. ;
  if (nang > 1) {
    dang=.5*pi/double(nang);
  }

  amu_0=1.0 ;
  amu=cos((double(indgen(nang))-0.5)*dang) ;

  pi0_0=0. ;
  pi1_0=1. ;
  pi0()=0. ;
  pi1()=1. ;

  s1_0= 0i ;
  s2_0= 0i ;
  nn=2*nang ;

  s1()=0.0+0.0*1i;
  s2()=0.0+0.0*1i;


  /* logarithmic derivative D(j) calculated by downward recurrence
     beginning with initial value (0.,0.) at j=nmx     */
  D(int(nmx))= 0 ;
  nn=nmx-1 ;
  for (n=1 ; n <= nn; n++) {
    i=int(nmx-n+1) ;
    C = double(i)/y
    D(i-1)= C - (1./(D(i)+C)) ;
  }


  /*  riccati-bessel functions with real argument x
      calculated by upward recurrence   */
  psi0=cos(xgrain) ;
  psi1=sin(xgrain) ;
  chi0=-sin(xgrain) ;
  chi1=cos(xgrain) ;
  xi1= psi1 -chi1 * 1i ;
  qsca=0. ;
  gsca=0. ;
  p= -1. ;
  for (n=1 ;  n <= nstop ; n++) {
    en=double(n) ;
    fn=(2.e0*en+1.)/(en*(en+1.)) ;
    /* for given n, psi  = psi_n        chi  = chi_n
                    psi1 = psi_{n-1}    chi1 = chi_{n-1}
                    psi0 = psi_{n-2}    chi0 = chi_{n-2}
       calculate psi_n and chi_n */
    psi=(2.*en-1.)*psi1/xgrain-psi0 ;
    chi=(2.*en-1.)*chi1/xgrain-chi0 ;
    xi= psi -chi *1i ;

    /* store previous values of an and bn for use
       in computation of g=<cos(theta)> */
    if(n > 1) {
      an1=an ;
      bn1=bn ;
    }

    //  compute an and bn:
    an=(D(n)/mgrain+en/xgrain)*psi-psi1 ;
    an=an/((D(n)/mgrain+en/xgrain)*xi-xi1) ;
    bn=(mgrain*D(n)+en/xgrain)*psi-psi1 ;
    bn=bn/((mgrain*D(n)+en/xgrain)*xi-xi1) ;

    // augment sums for qsca and g=<cos(theta)>
    qsca=qsca+(2.*en+1.)*(abs(an)^2+abs(bn)^2) ;
    gsca=gsca+((2.*en+1.)/(en*(en+1.)))*(double(an)*double(bn)+double(-1i*an)*double(-1i*bn)) ;
    if(n > 1) {
      gsca=gsca+((en-1.)*(en+1.)/en)*(double(an1)*double(an)+double(-1i*an1)*double(-1i*an) +double(bn1)*double(bn)+double(-1i*bn1)*double(-1i*bn)) ;
    }


    /* now calculate scattering intensity pattern
       0° calcule a part car les angles sont demi-entier
       et on a besoin de s1_0 pour qext */
    pi_0=pi1_0 ;
    tau_0=en*amu_0*pi_0-(en+1.)*pi0_0 ;
    s1_0=s1_0+fn*(an*pi_0+bn*tau_0) ;
    s2_0=s2_0+fn*(an*tau_0+bn*pi_0) ;
    //    first do angles from 0 to 90
    for (j=1 ;  j <= nang ; j++) {
      PI(j)=pi1(j) ;
      tau(j)=en*amu(j)*PI(j)-(en+1.)*pi0(j) ;
      s1(j)=s1(j)+fn*(an*PI(j)+bn*tau(j)) ;
      s2(j)=s2(j)+fn*(an*tau(j)+bn*PI(j)) ;
    }

    /* now do angles greater than 90 using pi and tau from
       angles less than 90.
       p=1 for n=1,3,...; p=-1 for n=2,4,... */
    p=-p ;
    for (j=1 ;  j <= nang ; j++) {
      jj=2*nang-j+1 ;
      s1(jj)=s1(jj)+fn*p*(an*PI(j)-bn*tau(j)) ;
      s2(jj)=s2(jj)+fn*p*(bn*PI(j)-an*tau(j)) ;
    }
    psi0=psi1 ;
    psi1=psi ;
    chi0=chi1 ;
    chi1=chi ;
    xi1= psi1 -chi1 * 1i;

    /*  compute pi_n for next value of n
        for each angle j, compute pi_n+1
        from pi = pi_n , pi0 = pi_n-1 */
    pi1_0=((2.*en+1.)*amu_0*pi_0-(en+1.)*pi0_0)/en ;
    pi0_0=pi_0 ;
    for (j=1 ;  j <= nang ; j++) {
      pi1(j)=((2.*en+1.)*amu(j)*PI(j)-(en+1.)*pi0(j))/en ;
      pi0(j)=PI(j) ;
    }
  } // boucle -> nstop

    /*  have summed sufficient terms.
        now compute qsca,qext,qback,and gsca*/
  gsca=2.*gsca/qsca ;
  qsca=(2./(xgrain*xgrain))*qsca ;
  qext=(4./(xgrain*xgrain))*double(s1_0) ;
  qback=(abs(s1(2*nang-1))/xgrain)^2/pi ;

  // fonction de phase normalisee de facon a ce que l'integrale face (0.5*xgrain^2*qsca)
  // (en integrant via sin(theta) dtheta)
  s11=0.5*(abs(s1)^2+abs(s2)^2);
  s12=0.5*(abs(s1)^2-abs(s2)^2);

}

//***********************************************************************

func hg(gsca,x) {
 /* DOCUMENT
    fonction de Heynyey-Greenstein
    arguments : gsca, x
    SEE ALSO:
    mie, bhmie
 */
  return ((1-gsca^2)/(2.0))*(1+gsca^2-2*gsca*cos(x/180.*pi))^(-1.5) ;
}

//***********************************************************************

func plot_phase_func(S11,HG=,width=,color1=,color2=) {
/* DOCUMENT
   Trace la fonction de phase et la Heynyey-Greenstein de même g
   argument : S11
   SEE ALSO:
   mie, bhmie, hg
 */
  if (is_void(width)) width=1 ;
  if (is_void(HG)) HG=0 ;
  if (is_void(color1)) color1="black" ;
  if (is_void(color2)) color2="black" ;


  // Tableau abscisses
  nang=dimsof(S11)(2)/2;
  dang=.5*pi/double(nang);
  tab_angle=(double(indgen(2*nang))-0.5)*dang;
  tab_cos = cos(tab_angle);
  tab_sin = sin(tab_angle);
  tab_angle=tab_angle*180./pi;

  // plot
  logxy, 0, 1;
  plg, S11 , tab_angle, color=color1,width=width;

  sum(S11 * tab_sin)

  if (HG) {
    // Calcul de g
    gsca = sum(S11() * tab_sin() * tab_cos());
    print, "g=", gsca;


    // Calcul de la HG
    fct_hg=hg(gsca,tab_angle);
    fct_hg = fct_hg/sum(fct_hg*tab_sin);

    plg, fct_hg, tab_angle, color=color2, type="dash",width=width;
  }
}


func plot_polarisabilite(S11,S12,width=,color=,type=) {
/* DOCUMENT
   Trace la fonction de phase et la Heynyey-Greenstein de même g
   argument : S11
   SEE ALSO:
   mie, bhmie, hg
 */

  if (is_void(width)) width=1;
  if (is_void(color)) color="black";


  // Tableau abscisses
  nang=dimsof(S11)(2)/2;
  dang=.5*pi/double(nang);
  tab_angle=(double(indgen(2*nang))-0.5)*dang;
  tab_angle=tab_angle*180./pi;

   // plot
  logxy, 0, 0;
  plg, S12/S11, tab_angle, width=width, color=color, type=type;
}
