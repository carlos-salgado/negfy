/*
 *                This source code is part of
 *
 *                     E  R  K  A  L  E
 *                             -
 *                       DFT from Hel
 *
 * Written by Susi Lehtola, 2010-2011
 * Copyright (c) 2010-2011, Susi Lehtola
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 */
#include "global.h"
#include <algorithm>
#include <armadillo>
#include <cstdio>
#include <cfloat>
#include "basis.h"
//#include "broyden.h"
#include "elements.h"
#include "dftfuncs.h"
#include "dftgrid.h"
//#include "diis.h"
#include "guess.h"
#include "linalg.h"
#include "mathf.h"
//#include "properties.h"
#include "scf.h"
#include "stringutil.h"
#include "timer.h"
//#include "trrh.h"
extern "C" {
#include <gsl/gsl_poly.h>
}
//# 1 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-force.cpp.in"
//# 1 "/usr/include/stdc-predef.h" 1

 
/*
 
 
//# 43 "/usr/include/stdc-predef.h"
//# 51 "/usr/include/stdc-predef.h"
 
 
//# 1 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-force.cpp.in" 2

 
 
arma::vec SCF::force_RHF(rscf_t & sol, const std::vector<double> & occs, double tol)
//# 33 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-force.cpp.in"
{
  arma::mat W;
  W=form_density(sol.E,sol.C,occs);
  
  arma::vec fpul_kin=basisp->kinetic_pulay(sol.P);
  
  arma::vec fpul_nuc=basisp->nuclear_pulay(sol.P);
  
  arma::vec fnuc=basisp->nuclear_der(sol.P);
  
  arma::vec forth=basisp->overlap_der(W);
  
  arma::vec frep=basisp->nuclear_force();
  
  double kfull=1.0;
  
  arma::vec fx_full;
  fx_full.zeros(fnuc.n_elem);
//# 90 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-force.cpp.in"
    if(!direct)
      scr.fill(basisp,intthr,verbose);
    fx_full=scr.forceJK(sol.P,tol,kfull);
//# 135 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-force.cpp.in"
  arma::vec ftot=fpul_kin+fpul_nuc+fnuc+forth+frep;
  arma::vec ffull=ftot+fx_full;
  
  
  
  return ffull;
}
//# 1 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-force.cpp.in"
//# 1 "/usr/include/stdc-predef.h" 1

 

 
 
//# 43 "/usr/include/stdc-predef.h"
//# 51 "/usr/include/stdc-predef.h"
 
 
//# 1 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-force.cpp.in" 2

 
 
//# 28 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-force.cpp.in"
arma::vec SCF::force_UHF(uscf_t & sol, const std::vector<double> & occa, const std::vector<double> & occb, double tol)
{
  arma::mat W;
  W=form_density(sol.Ea,sol.Ca,occa)+form_density(sol.Eb,sol.Cb,occb);
  
  arma::vec fpul_kin=basisp->kinetic_pulay(sol.P);
  
  arma::vec fpul_nuc=basisp->nuclear_pulay(sol.P);
  
  arma::vec fnuc=basisp->nuclear_der(sol.P);
  
  arma::vec forth=basisp->overlap_der(W);
  
  arma::vec frep=basisp->nuclear_force();
  
  double kfull=1.0;
  
  arma::vec fx_full;
  fx_full.zeros(fnuc.n_elem);
//# 90 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-force.cpp.in"
    if(!direct)
      scr.fill(basisp,intthr,verbose);
    fx_full=scr.forceJK(sol.Pa,sol.Pb,tol,kfull);
//# 135 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-force.cpp.in"
  arma::vec ftot=fpul_kin+fpul_nuc+fnuc+forth+frep;
  arma::vec ffull=ftot+fx_full;
  
  
  
  return ffull;
}
//# 1 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-force.cpp.in"
//# 1 "/usr/include/stdc-predef.h" 1

 

 
 
//# 43 "/usr/include/stdc-predef.h"
//# 51 "/usr/include/stdc-predef.h"
 
 
//# 1 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-force.cpp.in" 2

 
 
arma::vec SCF::force_RDFT(rscf_t & sol, const std::vector<double> & occs, const dft_t dft, DFTGrid & grid, DFTGrid & nlgrid, double tol)
//# 33 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-force.cpp.in"
{
  arma::mat W;
  W=form_density(sol.E,sol.C,occs);
  
  arma::vec fpul_kin=basisp->kinetic_pulay(sol.P);
  
  arma::vec fpul_nuc=basisp->nuclear_pulay(sol.P);
  
  arma::vec fnuc=basisp->nuclear_der(sol.P);
  
  arma::vec forth=basisp->overlap_der(W);
  
  arma::vec frep=basisp->nuclear_force();
  
  
  double omega, kfull, kshort;
  range_separation(dft.x_func,omega,kfull,kshort);
  
  arma::vec fx_full;
  fx_full.zeros(fnuc.n_elem);
  
  arma::vec fx_short;
  fx_short.zeros(fnuc.n_elem);
  if(kfull==0.0) {
    if(densityfit) {
      if(kshort!=0.0 || kfull!=0.0)
	throw std::runtime_error("Forces not implemented for density fitting of exact exchange.\n");
      fx_full=dfit.forceJ(sol.P);
    } else {
      if(!direct)
	scr.fill(basisp,intthr,verbose);
      fx_full=scr.forceJ(sol.P,tol);
    }
  } else {
    if(!direct)
      scr.fill(basisp,intthr,verbose);
    fx_full=scr.forceJK(sol.P,tol,kfull);
  }
  if(omega != 0.0) {
    
    scr_rs.set_range_separation(omega,0.0,1.0);
    scr_rs.fill(basisp,intthr,verbose);
    
    fx_short=scr_rs.forceK(sol.P,tol,kshort);
  }
  
  arma::vec fxc;
  fxc.zeros(fnuc.n_elem);
  if(dft.x_func>0 || dft.c_func>0) {
    fxc=grid.eval_force(dft.x_func,dft.c_func,sol.P);
  }
  
  if(dft.nl) {
    arma::vec vv10f(grid.eval_VV10_force(nlgrid,dft.vv10_b,dft.vv10_C,sol.P));
    
    fxc+=vv10f;
  }
  arma::vec ftot=fpul_kin+fpul_nuc+fnuc+forth+frep;
  arma::vec ffull=ftot+fx_full;
  
  ffull+=fxc+fx_short;
  
  
  
  return ffull;
}
//# 1 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-force.cpp.in"
//# 1 "/usr/include/stdc-predef.h" 1


 

 
 
//# 43 "/usr/include/stdc-predef.h"
//# 51 "/usr/include/stdc-predef.h"
 
 
//# 1 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-force.cpp.in" 2

 */




 
//# 25 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-force.cpp.in"
arma::vec SCF::force_UDFT(uscf_t & sol, const std::vector<double> & occa, const std::vector<double> & occb, const dft_t dft, DFTGrid & grid, DFTGrid & nlgrid, double tol)
//# 33 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-force.cpp.in"
{
  arma::mat W;
  W=form_density(sol.Ea,sol.Ca,occa)+form_density(sol.Eb,sol.Cb,occb);
  W.print(std::cout);
  
  arma::vec fpul_kin=basisp->kinetic_pulay(sol.P);
  
  arma::vec fpul_nuc=basisp->nuclear_pulay(sol.P);
  
  arma::vec fnuc=basisp->nuclear_der(sol.P);
  
  arma::vec forth=basisp->overlap_der(W);
  
  arma::vec frep=basisp->nuclear_force();
  
  
  double omega, kfull, kshort;
  range_separation(dft.x_func,omega,kfull,kshort);
  
  arma::vec fx_full;
  fx_full.zeros(fnuc.n_elem);
  
  arma::vec fx_short;
  fx_short.zeros(fnuc.n_elem);
  /*
  if(kfull==0.0) {
    if(densityfit) {
      if(kshort!=0.0 || kfull!=0.0)
	throw std::runtime_error("Forces not implemented for density fitting of exact exchange.\n");
      fx_full=dfit.forceJ(sol.P);
    } else {
      if(!direct)
	scr.fill(basisp,intthr,verbose);
      fx_full=scr.forceJ(sol.P,tol);
    }
  } else {
    if(!direct)
      scr.fill(basisp,intthr,verbose);
    fx_full=scr.forceJK(sol.Pa,sol.Pb,tol,kfull);
  }
  if(omega != 0.0) {
    
    scr_rs.set_range_separation(omega,0.0,1.0);
    scr_rs.fill(basisp,intthr,verbose);
    
    fx_short=scr_rs.forceK(sol.Pa,sol.Pb,tol,kshort);
  }
  */
  
  arma::vec fxc;
  fxc.zeros(fnuc.n_elem);
  if(dft.x_func>0 || dft.c_func>0) {
    fxc=grid.eval_force(dft.x_func,dft.c_func,sol.Pa,sol.Pb);
  }
  
  if(dft.nl) {
    arma::vec vv10f(grid.eval_VV10_force(nlgrid,dft.vv10_b,dft.vv10_C,sol.P));
    
    fxc+=vv10f;
  }
  arma::vec ftot=fpul_kin+fpul_nuc+fnuc+forth+frep;
  arma::vec ffull=ftot+fx_full;
  
  ffull+=fxc+fx_short;
  
  
  
  return ffull;
}

 //# 25 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-force.cpp.in"
 arma::vec SCF::force_UDFT_ANT(uscf_t & sol, const std::vector<double> & occa, const std::vector<double> & occb, const dft_t dft, DFTGrid & grid, DFTGrid & nlgrid, double tol)
 //# 33 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-force.cpp.in"
 {
   std::cout << "\nENTER arma::vec SCF::force_UDFT_ANT 1st version\n";
   arma::mat W;
   W=form_density(sol.Ea,sol.Ca,occa)+form_density(sol.Eb,sol.Cb,occb);
   W.print(std::cout);

   arma::vec fpul_kin=basisp->kinetic_pulay(sol.P);

   arma::vec fpul_nuc=basisp->nuclear_pulay(sol.P);

   arma::vec fnuc=basisp->nuclear_der(sol.P);

   arma::vec forth=basisp->overlap_der(W);

   arma::vec frep=basisp->nuclear_force();


   double omega, kfull, kshort;
   range_separation(dft.x_func,omega,kfull,kshort);

   arma::vec fx_full;
   fx_full.zeros(fnuc.n_elem);

   arma::vec fx_short;
   fx_short.zeros(fnuc.n_elem);
   /*
   if(kfull==0.0) {
     if(densityfit) {
       if(kshort!=0.0 || kfull!=0.0)
 	throw std::runtime_error("Forces not implemented for density fitting of exact exchange.\n");
       fx_full=dfit.forceJ(sol.P);
     } else {
       if(!direct)
 	scr.fill(basisp,intthr,verbose);
       fx_full=scr.forceJ(sol.P,tol);
     }
   } else {
     if(!direct)
       scr.fill(basisp,intthr,verbose);
     fx_full=scr.forceJK(sol.Pa,sol.Pb,tol,kfull);
   }
   if(omega != 0.0) {

     scr_rs.set_range_separation(omega,0.0,1.0);
     scr_rs.fill(basisp,intthr,verbose);

     fx_short=scr_rs.forceK(sol.Pa,sol.Pb,tol,kshort);
   }
   */

   arma::vec fxc;
   fxc.zeros(fnuc.n_elem);
   if(dft.x_func>0 || dft.c_func>0) {
     fxc=grid.eval_force(dft.x_func,dft.c_func,sol.Pa,sol.Pb);
   }

   if(dft.nl) {
     arma::vec vv10f(grid.eval_VV10_force(nlgrid,dft.vv10_b,dft.vv10_C,sol.P));

     fxc+=vv10f;
   }
   arma::vec ftot=fpul_kin+fpul_nuc+fnuc+forth+frep;
   //arma::vec ffull=ftot+fx_full; // ORIGINAL CODE COMMENTED BY C.SALGADO ON 2019-02-26 TO RETURN ONLY EXCHANGE CORELATION EFFECTS.
   arma::vec ffull=fx_full;

   ffull+=fxc+fx_short;

   std::cout << "\nEXIT arma::vec SCF::force_UDFT_ANT 1st version\n";
   return ffull;
 }

 //# 25 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-force.cpp.in"
 //arma::vec SCF::force_UDFT_ANT(uscf_t & sol, const std::vector<double> & occa, const std::vector<double> & occb, const dft_t dft, DFTGrid & grid, DFTGrid & nlgrid, double tol, const std::vector< std::vector<arma::mat> > & fmat)
 arma::vec SCF::force_UDFT_ANT(uscf_t & sol, const std::vector<double> & occa, const std::vector<double> & occb, const dft_t dft, DFTGrid & grid, DFTGrid & nlgrid, double tol, std::vector< std::vector<arma::mat> > & fmat)
 //# 33 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-force.cpp.in"
 {
   std::cout << "\nENTER arma::vec SCF::force_UDFT_ANT 2nd version\n";
   arma::mat W;
   W=form_density(sol.Ea,sol.Ca,occa)+form_density(sol.Eb,sol.Cb,occb);
   W.print(std::cout);

   /*
   arma::vec fpul_kin=basisp->kinetic_pulay(sol.P);

   arma::vec fpul_nuc=basisp->nuclear_pulay(sol.P);
   */
   arma::vec fnuc=basisp->nuclear_der(sol.P);
   /*
   arma::vec forth=basisp->overlap_der(W);

   arma::vec frep=basisp->nuclear_force();

   std::cout << "\nAFTER basisp QUANTITIES\n";
   */

   basisp->print(true);


   double omega, kfull, kshort;
   range_separation(dft.x_func,omega,kfull,kshort);

   std::cout << "\nAFTER range_separation(dft.x_func,omega,kfull,kshort);\n";

   arma::vec fx_full;
   fx_full.zeros(fnuc.n_elem);

   arma::vec fx_short;
   fx_short.zeros(fnuc.n_elem);
   /*
   if(kfull==0.0) {
     if(densityfit) {
       if(kshort!=0.0 || kfull!=0.0)
 	throw std::runtime_error("Forces not implemented for density fitting of exact exchange.\n");
       fx_full=dfit.forceJ(sol.P);
     } else {
       if(!direct)
 	scr.fill(basisp,intthr,verbose);
       fx_full=scr.forceJ(sol.P,tol);
     }
   } else {
     if(!direct)
       scr.fill(basisp,intthr,verbose);
     fx_full=scr.forceJK(sol.Pa,sol.Pb,tol,kfull);
   }
   if(omega != 0.0) {

     scr_rs.set_range_separation(omega,0.0,1.0);
     scr_rs.fill(basisp,intthr,verbose);

     fx_short=scr_rs.forceK(sol.Pa,sol.Pb,tol,kshort);
   }
   */

   std::cout << "\ndft.x_func = "<<dft.x_func<<"; dft.c_func = "<<dft.c_func<<"\n";
   arma::vec fxc;
   fxc.zeros(fnuc.n_elem);
   if(dft.x_func>0 || dft.c_func>0) {
     //fxc=grid.eval_force(dft.x_func,dft.c_func,sol.Pa,sol.Pb);
	 fxc=grid.eval_force(dft.x_func,dft.c_func,sol.Pa,sol.Pb,fmat); //3*basp->get_Nnuc()
   }

   if(dft.nl) {
     arma::vec vv10f(grid.eval_VV10_force(nlgrid,dft.vv10_b,dft.vv10_C,sol.P));

     fxc+=vv10f;
   }
   //arma::vec ftot=fpul_kin+fpul_nuc+fnuc+forth+frep; // ORIGINAL CODE COMMENTED BY C.SALGADO ON 2019-02-26 TO RETURN ONLY EXCHANGE CORELATION EFFECTS.
   //arma::vec ffull=ftot+fx_full; // ORIGINAL CODE COMMENTED BY C.SALGADO ON 2019-02-26 TO RETURN ONLY EXCHANGE CORELATION EFFECTS.
   arma::vec ffull=fx_full;

   ffull+=fxc+fx_short;

   std::cout << "\nPRINT OUTPUT ffull\n";
   std::cout << "\nffull EN scf-force DEBERÍA COINCIDIR CON F2 en HartreeFockClass 2-body Erkale forces\n";
   ffull.print("ffull");
   std::cout << "\nEXIT arma::vec SCF::force_UDFT_ANT 2nd version\n";
   return ffull;
 }
