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
#include "broyden.h"
#include "elements.h"
#include "dftfuncs.h"
#include "dftgrid.h"
#include "diis.h"
#include "guess.h"
#include "linalg.h"
#include "mathf.h"
#include "properties.h"
#include "scf.h"
#include "stringutil.h"
#include "timer.h"
#include "trrh.h"

extern "C" {
#include <gsl/gsl_poly.h>
}
//# 1 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
//# 1 "/usr/include/stdc-predef.h" 1















 










 






 

//# 43 "/usr/include/stdc-predef.h"

//# 51 "/usr/include/stdc-predef.h"



 


 


//# 1 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in" 2














 





void SCF::RHF(rscf_t & sol, const std::vector<double> & occs, double convthr)

//# 38 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
{

  
//# 47 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

//# 65 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

  




  
//# 79 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

  
//# 93 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"








  



  size_t nocc;
  for(nocc=occs.size()-1;nocc<occs.size();nocc--)
    if(occs[nocc]>0)
      break;
  nocc++;
//# 124 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

  int nfock=0;

  Timer t;
  Timer ttot;

  
  double diiserr=DBL_MAX;
  
  bool diissucc=false;

  
  arma::mat orbs;
  arma::mat Horth;


  
  rDIIS diis(S,Sinvh,usediis,diiseps,diisthr,useadiis,verbose,diisorder);
  
  Broyden broyd(verbose);

  
  arma::mat J(Nbf,Nbf), K(Nbf,Nbf);
  J.zeros();
  K.zeros();
//# 158 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

  
  arma::vec dipmom;

  
  double deltaE=0;

  
  double rmsdiff=0.0, maxdiff=0.0;






  if(sol.P.n_rows!=Nbf)



      {
	throw std::runtime_error("No starting guess provided for SCF solver!\n");
      }





  
  if(verbose) {

    if(sol.E.n_elem)
      print_E(sol.E,occs,false);
//# 198 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
    fflush(stdout);
  }

//# 226 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"


  
  if(verbose && maxiter>0) {
    fprintf(stderr,"Running ");

    fprintf(stderr,"restricted ");
//# 242 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"


    fprintf(stderr,"HF ");
//# 252 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

    fprintf(stderr,"calculation");
    if(densityfit) {
      fprintf(stderr," with density fitting");
    }

//# 274 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
    fprintf(stderr,".\n%4s %16s %10s %9s %9s %9s %10s\n","iter","E","dE","RMS dens","MAX dens","DIIS err","titer (s)");
  }
  fflush(stdout);

//# 302 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

  
  bool converged=false;

  

  rscf_t oldsol;



  oldsol.en.E=0.0;

  
//# 322 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
  std::vector<double> occr(occs);
  while(occr.size()<Sinvh.n_cols)
    occr.push_back(0.0);


  
  
  if(maxiter>0) {
    chkptp->open();
    chkptp->write("tol",intthr);
    chkptp->write("P",sol.P);
    chkptp->write(sol.en);
//# 350 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
    chkptp->write("C",sol.C);
    chkptp->write("E",sol.E);
    chkptp->write("occs",occr);
    
    chkptp->write("Restricted",1);

    chkptp->write("Converged",converged);
    chkptp->close();

  } else {
    Timer tg;
    
    Fock_RHF(sol,occs); nfock++;;
    
    if(verbose) {
      if(shift==0.0)
	printf("Solving SCF equations ... ");
      else
	printf("Solving SCF equations with level shift %.3f ... ",shift);
      fflush(stdout);
      t.set();
    }

    
    dipmom=dipole_moment(sol.P,*basisp);
    
    diagonalize(sol,shift);

    if(verbose)
      printf("done (%s)\n",t.elapsed().c_str());
  }

  
  int iiter=1;

//# 391 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

  while(iiter<=maxiter) {
    Timer titer;

    if(verbose)
      printf("\n ******* Iteration %4i ********\n",iiter);

    
    Fock_RHF(sol,occs); nfock++;;

    
    deltaE=sol.en.E-oldsol.en.E;


    
    diis.update(sol.H,sol.P,sol.en.E,diiserr);





    if(iiter>1 && usebroyden) {

      
      broyd.push_x(MatToVec(oldsol.H));
      broyd.push_f(MatToVec(oldsol.H-sol.H));
//# 432 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
    }

    
    try {

      
      diis.solve_F(sol.H);




      diissucc=true;
    } catch(std::runtime_error &) {
      diissucc=false;
    }

    
    if(usebroyden && !diissucc && iiter>1) {

      if(verbose) {
	printf("Performing Broyden interpolation of Fock operator ... ");
	fflush(stdout);
	t.set();
      }


      
      sol.H=VecToMat(broyd.update_x(),Nbf,Nbf);
//# 468 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

      if(verbose)
	printf("done (%s)\n",t.elapsed().c_str());
    }

    
    oldsol=sol;

    if(usetrrh) {
      
      if(verbose) {
	printf("\nSolving TRRH equations.\n");
	fflush(stdout);
	t.set();
      }


      arma::mat Cnew;
      arma::vec Enew;
      TRRH_update(sol.H,sol.C,S,Cnew,Enew,nocc,verbose,trrhmins);

      
      sol.C=Cnew;
      sol.E=Enew;

      
      check_orth(sol.C,S,false);
//# 515 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

      if(verbose)
	printf("TRRH solved in %s.\n\n",t.elapsed().c_str());

    } else {
      
      if(verbose) {
	if(shift==0.0)
	  printf("\nSolving SCF equations ... ");
	else
	  printf("\nSolving SCF equations with level shift %.3f ... ",shift);
	fflush(stdout);
	t.set();
      }

      
      diagonalize(sol,shift);

      if(verbose)
	printf("done (%s)\n",t.elapsed().c_str());
    }

    
    form_density(sol,occs);;

    
    arma::mat deltaP=sol.P-oldsol.P;





    
    dipmom=dipole_moment(sol.P,*basisp);

    
    maxdiff=max_abs(deltaP/2.0);
    rmsdiff=rms_norm(deltaP/2.0);

    
    double maxdiff_cvd(maxdiff);
    double rmsdiff_cvd(maxdiff);

//# 568 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

    
    if(verbose) {
      printf("\n");
      printf("%-30s: % .16e\n","Total energy",sol.en.E);
      printf("%-30s: % e\n","DIIS error",diiserr);
      printf("%-30s: % e\n","Energy change",deltaE);
      printf("%-30s: % e\n","Max total density change",maxdiff);
      printf("%-30s: % e\n","Max rms   density change",rmsdiff);
//# 583 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
      printf("Dipole mu = (% 08.8f, % 08.8f, % 08.8f) D\n",dipmom(0)/AUINDEBYE,dipmom(1)/AUINDEBYE,dipmom(2)/AUINDEBYE);

      printf("\nIteration took %s.\n",titer.elapsed().c_str());
      fflush(stdout);
    }

    
    if(verbose)
      fprintf(stderr,"%4i % 16.8f % 10.3e %9.3e %9.3e %9.3e %10.3f\n",iiter,sol.en.E,deltaE,rmsdiff_cvd,maxdiff_cvd,diiserr,titer.get());

//# 605 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
    if(diiserr < convthr) {
      converged=true;
      if(verbose)
        printf("\n ******* Convergence achieved ********\n");
    }

    
    
    chkptp->open();
    chkptp->write("P",sol.P);
    chkptp->write(sol.en);
//# 635 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
    chkptp->write("C",sol.C);
    chkptp->write("E",sol.E);
    chkptp->write("H",sol.H);
    chkptp->write("occs",occr);
    
    chkptp->write("Restricted",1);

    chkptp->write("Converged",converged);
    chkptp->close();

    
    if(converged)
      break;

    iiter++;
  } 

  if(verbose) {

    if(converged) {
      std::string method=

	"R"
//# 666 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
	"HF"

	;

      printf("Solution of %s took %s.\n",method.c_str(),ttot.elapsed().c_str());
      fprintf(stderr,"%s converged in %s.\n",method.c_str(),ttot.elapsed().c_str());
    }

    
    printf("\n");
    printf("%-21s energy: % .16e\n","Kinetic",sol.en.Ekin);
    printf("%-21s energy: % .16e\n","Nuclear attraction",sol.en.Enuca);
    printf("%-21s energy: % .16e\n","Total one-electron",sol.en.Eone);
    printf("%-21s energy: % .16e\n","Nuclear repulsion",sol.en.Enucr);
    printf("%-21s energy: % .16e\n","Coulomb",sol.en.Ecoul);

    printf("%-21s energy: % .16e\n","Exchange",sol.en.Exc);




    printf("-----------------------------------------------------\n");
    printf("%28s: % .16e\n","Total energy",sol.en.E);
    printf("%28s: % .16e\n","Virial factor",-sol.en.E/sol.en.Ekin);

    printf("\nDipole mu = (% 08.8f, % 08.8f, % 08.8f) D\n",dipmom(0)/AUINDEBYE,dipmom(1)/AUINDEBYE,dipmom(2)/AUINDEBYE);

    printf("\n");
    

    print_E(sol.E,occs,true);
//# 703 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"


    if(dimcalc) {
      
      arma::ivec mv(basisp->get_m_values());
      arma::umat occ;


      occ.zeros(arma::max(arma::abs(mv))+1,1);

      
      arma::ivec mc(m_classify(sol.C,basisp->get_m_values()));
      for(size_t i=0;i<nocc;i++)
        occ(std::abs(mc(i)),0)++;
//# 728 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

      printf("\nOrbital occupations by symmetry\n");
      for(size_t m=0;m<occ.n_rows;m++) {

        if(occ(m,0)>0)
          printf("m = %i: %2i\n",(int) m,(int) occ(m,0));




      }
    }


  }

  if(converged) {
    if(doforce) {
      arma::vec f;

      f=force_RHF(sol,occs,intthr);
//# 763 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

      chkptp->write("Force",f);
    }

  } else if(maxiter>0) {
    std::ostringstream oss;
    oss << "Error in function " << __FUNCTION__ << " (file " << "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in" << ", near line " << 769 << "): SCF did not converge in "<<maxiter<<" iterations!\n";
    throw std::runtime_error(oss.str());
  }




}
//# 1 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
//# 1 "/usr/include/stdc-predef.h" 1















 










 






 

//# 43 "/usr/include/stdc-predef.h"

//# 51 "/usr/include/stdc-predef.h"



 


 


//# 1 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in" 2














 

//# 27 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
void SCF::UHF(uscf_t & sol, const std::vector<double> & occa, const std::vector<double> & occb, double convthr)

//# 38 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
{

  
//# 49 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
  

//# 65 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

  




  
//# 79 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

  




















  


//# 111 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
  size_t nocca;
  for(nocca=occa.size()-1;nocca<occa.size();nocca--)
    if(occa[nocca]>0)
      break;
  nocca++;

  size_t noccb;
  for(noccb=occb.size()-1;noccb<occb.size();noccb--)
    if(occb[noccb]>0)
      break;
  noccb++;



  int nfock=0;

  Timer t;
  Timer ttot;

  
  double diiserr=DBL_MAX;
  
  bool diissucc=false;

  
  arma::mat orbs;
  arma::mat Horth;

//# 150 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
  
  uDIIS diis(S,Sinvh,diiscomb,usediis,diiseps,diisthr,useadiis,verbose,diisorder);

  
  Broyden broyd_sum(verbose);
  Broyden broyd_diff(verbose);



  
  arma::vec dipmom;

  
  double deltaE=0;

  
  double rmsdiff=0.0, maxdiff=0.0;

  double rmsdiffa=0.0, maxdiffa=0.0;
  double rmsdiffb=0.0, maxdiffb=0.0;





    if(sol.Pa.n_rows!=Nbf || sol.Pb.n_rows!=Nbf)

      {
	throw std::runtime_error("No starting guess provided for SCF solver!\n");
      }





  
  if(verbose) {




    if(sol.Ea.n_elem) {
      printf("alpha: ");
      print_E(sol.Ea,occa,false);
      printf("beta:  ");
      print_E(sol.Eb,occb,false);
    }

    fflush(stdout);
  }

//# 226 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"


  
  if(verbose && maxiter>0) {
    fprintf(stderr,"Running ");
//# 238 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
    fprintf(stderr,"unrestricted ");





    fprintf(stderr,"HF ");
//# 252 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

    fprintf(stderr,"calculation");
    if(densityfit) {
      fprintf(stderr," with density fitting");
    }

//# 274 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
    fprintf(stderr,".\n%4s %16s %10s %9s %9s %9s %10s\n","iter","E","dE","RMS dens","MAX dens","DIIS err","titer (s)");
  }
  fflush(stdout);

//# 302 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

  
  bool converged=false;

  



  uscf_t oldsol;

  oldsol.en.E=0.0;

  

  std::vector<double> occar(occa), occbr(occb);
  while(occar.size()<Sinvh.n_cols)
    occar.push_back(0.0);
  while(occbr.size()<Sinvh.n_cols)
    occbr.push_back(0.0);






  
  
  if(maxiter>0) {
    chkptp->open();
    chkptp->write("tol",intthr);
    chkptp->write("P",sol.P);
    chkptp->write(sol.en);

    chkptp->write("Ca",sol.Ca);
    chkptp->write("Cb",sol.Cb);

    chkptp->write("Ea",sol.Ea);
    chkptp->write("Eb",sol.Eb);

    chkptp->write("Pa",sol.Pa);
    chkptp->write("Pb",sol.Pb);

    chkptp->write("occa",occar);
    chkptp->write("occb",occbr);

    
    chkptp->write("Restricted",0);
//# 356 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
    chkptp->write("Converged",converged);
    chkptp->close();

  } else {
    Timer tg;
    
    Fock_UHF(sol,occa,occb); nfock++;;
    
    if(verbose) {
      if(shift==0.0)
	printf("Solving SCF equations ... ");
      else
	printf("Solving SCF equations with level shift %.3f ... ",shift);
      fflush(stdout);
      t.set();
    }

    
    dipmom=dipole_moment(sol.P,*basisp);
    
    diagonalize(sol,shift);

    if(verbose)
      printf("done (%s)\n",t.elapsed().c_str());
  }

  
  int iiter=1;

//# 391 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

  while(iiter<=maxiter) {
    Timer titer;

    if(verbose)
      printf("\n ******* Iteration %4i ********\n",iiter);

    
    Fock_UHF(sol,occa,occb); nfock++;;

    
    deltaE=sol.en.E-oldsol.en.E;





    
    diis.update(sol.Ha,sol.Hb,sol.Pa,sol.Pb,sol.en.E,diiserr);


    if(iiter>1 && usebroyden) {





      
      arma::mat Hs=sol.Ha+sol.Hb;
      arma::mat Hd=sol.Ha-sol.Hb;

      arma::mat Hsold=oldsol.Ha+oldsol.Hb;
      arma::mat Hdold=oldsol.Ha-oldsol.Hb;

      
      broyd_sum.push_x(MatToVec(Hsold));
      broyd_sum.push_f(MatToVec(Hsold-Hs));

      broyd_diff.push_x(MatToVec(Hdold));
      broyd_diff.push_f(MatToVec(Hdold-Hd));

    }

    
    try {




      
      diis.solve_F(sol.Ha,sol.Hb);

      diissucc=true;
    } catch(std::runtime_error &) {
      diissucc=false;
    }

    
    if(usebroyden && !diissucc && iiter>1) {

      if(verbose) {
	printf("Performing Broyden interpolation of Fock operator ... ");
	fflush(stdout);
	t.set();
      }





      arma::mat Hs=VecToMat(broyd_sum.update_x(),Nbf,Nbf);
      arma::mat Hd=VecToMat(broyd_diff.update_x(),Nbf,Nbf);

      
      sol.Ha=0.5*(Hs+Hd);
      sol.Hb=0.5*(Hs-Hd);


      if(verbose)
	printf("done (%s)\n",t.elapsed().c_str());
    }

    
    oldsol=sol;

    if(usetrrh) {
      
      if(verbose) {
	printf("\nSolving TRRH equations.\n");
	fflush(stdout);
	t.set();
      }

//# 496 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
      arma::mat Canew;
      arma::vec Eanew;
      TRRH_update(sol.Ha,sol.Ca,S,Canew,Eanew,nocca,verbose,trrhmins);

      arma::mat Cbnew;
      arma::vec Ebnew;
      TRRH_update(sol.Hb,sol.Cb,S,Cbnew,Ebnew,noccb,verbose,trrhmins);

      
      sol.Ca=Canew;
      sol.Cb=Cbnew;

      sol.Ea=Eanew;
      sol.Eb=Ebnew;

      
      check_orth(sol.Ca,S,false);
      check_orth(sol.Cb,S,false);


      if(verbose)
	printf("TRRH solved in %s.\n\n",t.elapsed().c_str());

    } else {
      
      if(verbose) {
	if(shift==0.0)
	  printf("\nSolving SCF equations ... ");
	else
	  printf("\nSolving SCF equations with level shift %.3f ... ",shift);
	fflush(stdout);
	t.set();
      }

      
      diagonalize(sol,shift);

      if(verbose)
	printf("done (%s)\n",t.elapsed().c_str());
    }

    
    form_density(sol,occa,occb);;

    
    arma::mat deltaP=sol.P-oldsol.P;

    arma::mat deltaPa=sol.Pa-oldsol.Pa;
    arma::mat deltaPb=sol.Pb-oldsol.Pb;


    
    dipmom=dipole_moment(sol.P,*basisp);

    
    maxdiff=max_abs(deltaP/2.0);
    rmsdiff=rms_norm(deltaP/2.0);

    
    double maxdiff_cvd(maxdiff);
    double rmsdiff_cvd(maxdiff);


    maxdiffa=max_abs(deltaPa);
    maxdiffb=max_abs(deltaPb);

    rmsdiffa=rms_norm(deltaPa);
    rmsdiffb=rms_norm(deltaPb);

    maxdiff_cvd=std::max(maxdiffa,maxdiffb);
    rmsdiff_cvd=std::max(rmsdiffa,rmsdiffb);


    
    if(verbose) {
      printf("\n");
      printf("%-30s: % .16e\n","Total energy",sol.en.E);
      printf("%-30s: % e\n","DIIS error",diiserr);
      printf("%-30s: % e\n","Energy change",deltaE);
      printf("%-30s: % e\n","Max total density change",maxdiff);
      printf("%-30s: % e\n","Max rms   density change",rmsdiff);

      printf("%-30s: % e\n","Max total alpha density change",maxdiffa);
      printf("%-30s: % e\n","Max rms   alpha density change",rmsdiffa);
      printf("%-30s: % e\n","Max total beta  density change",maxdiffb);
      printf("%-30s: % e\n","Max rms   beta  density change",rmsdiffb);

      printf("Dipole mu = (% 08.8f, % 08.8f, % 08.8f) D\n",dipmom(0)/AUINDEBYE,dipmom(1)/AUINDEBYE,dipmom(2)/AUINDEBYE);

      printf("\nIteration took %s.\n",titer.elapsed().c_str());
      fflush(stdout);
    }

    
    if(verbose)
      fprintf(stderr,"%4i % 16.8f % 10.3e %9.3e %9.3e %9.3e %10.3f\n",iiter,sol.en.E,deltaE,rmsdiff_cvd,maxdiff_cvd,diiserr,titer.get());

//# 605 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
    if(diiserr < convthr) {
      converged=true;
      if(verbose)
        printf("\n ******* Convergence achieved ********\n");
    }

    
    
    chkptp->open();
    chkptp->write("P",sol.P);
    chkptp->write(sol.en);

    chkptp->write("Ca",sol.Ca);
    chkptp->write("Cb",sol.Cb);

    chkptp->write("Ea",sol.Ea);
    chkptp->write("Eb",sol.Eb);

    chkptp->write("Ha",sol.Ha);
    chkptp->write("Hb",sol.Hb);

    chkptp->write("Pa",sol.Pa);
    chkptp->write("Pb",sol.Pb);

    chkptp->write("occa",occar);
    chkptp->write("occb",occbr);

    
    chkptp->write("Restricted",0);
//# 642 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
    chkptp->write("Converged",converged);
    chkptp->close();

    
    if(converged)
      break;

    iiter++;
  } 

  if(verbose) {

    if(converged) {
      std::string method=





	"U"




	"HF"

	;

      printf("Solution of %s took %s.\n",method.c_str(),ttot.elapsed().c_str());
      fprintf(stderr,"%s converged in %s.\n",method.c_str(),ttot.elapsed().c_str());
    }

    
    printf("\n");
    printf("%-21s energy: % .16e\n","Kinetic",sol.en.Ekin);
    printf("%-21s energy: % .16e\n","Nuclear attraction",sol.en.Enuca);
    printf("%-21s energy: % .16e\n","Total one-electron",sol.en.Eone);
    printf("%-21s energy: % .16e\n","Nuclear repulsion",sol.en.Enucr);
    printf("%-21s energy: % .16e\n","Coulomb",sol.en.Ecoul);

    printf("%-21s energy: % .16e\n","Exchange",sol.en.Exc);




    printf("-----------------------------------------------------\n");
    printf("%28s: % .16e\n","Total energy",sol.en.E);
    printf("%28s: % .16e\n","Virial factor",-sol.en.E/sol.en.Ekin);

    printf("\nDipole mu = (% 08.8f, % 08.8f, % 08.8f) D\n",dipmom(0)/AUINDEBYE,dipmom(1)/AUINDEBYE,dipmom(2)/AUINDEBYE);

    printf("\n");
    



    printf("alpha: ");
    print_E(sol.Ea,occa,true);
    printf("beta:  ");
    print_E(sol.Eb,occb,true);



    if(dimcalc) {
      
      arma::ivec mv(basisp->get_m_values());
      arma::umat occ;

//# 718 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
      occ.zeros(arma::max(arma::abs(mv))+1,2);

      
      arma::ivec mca(m_classify(sol.Ca,mv));
      for(size_t i=0;i<nocca;i++)
        occ(std::abs(mca(i)),0)++;
      arma::ivec mcb(m_classify(sol.Cb,mv));
      for(size_t i=0;i<noccb;i++)
        occ(std::abs(mcb(i)),1)++;


      printf("\nOrbital occupations by symmetry\n");
      for(size_t m=0;m<occ.n_rows;m++) {




        if(occ(m,0)+occ(m,1)>0)
          printf("m = %i: %2i %2i\n",(int) m,(int) occ(m,0),(int) occ(m,1));

      }
    }


  }

  if(converged) {
    if(doforce) {
      arma::vec f;




      f=force_UHF(sol,occa,occb,intthr);
//# 763 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

      chkptp->write("Force",f);
    }

  } else if(maxiter>0) {
    std::ostringstream oss;
    oss << "Error in function " << __FUNCTION__ << " (file " << "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in" << ", near line " << 769 << "): SCF did not converge in "<<maxiter<<" iterations!\n";
    throw std::runtime_error(oss.str());
  }




}
//# 1 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
//# 1 "/usr/include/stdc-predef.h" 1















 










 






 

//# 43 "/usr/include/stdc-predef.h"

//# 51 "/usr/include/stdc-predef.h"



 


 


//# 1 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in" 2














 

//# 30 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
void SCF::ROHF(uscf_t & sol, const std::vector<double> & occa, const std::vector<double> & occb, double convthr)

//# 38 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
{

  
//# 49 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
  

//# 65 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

  




  
//# 79 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

  




















  


//# 111 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
  size_t nocca;
  for(nocca=occa.size()-1;nocca<occa.size();nocca--)
    if(occa[nocca]>0)
      break;
  nocca++;

  size_t noccb;
  for(noccb=occb.size()-1;noccb<occb.size();noccb--)
    if(occb[noccb]>0)
      break;
  noccb++;



  int nfock=0;

  Timer t;
  Timer ttot;

  
  double diiserr=DBL_MAX;
  
  bool diissucc=false;

  
  arma::mat orbs;
  arma::mat Horth;

//# 150 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
  
  uDIIS diis(S,Sinvh,diiscomb,usediis,diiseps,diisthr,useadiis,verbose,diisorder);

  
  Broyden broyd_sum(verbose);
  Broyden broyd_diff(verbose);



  
  arma::vec dipmom;

  
  double deltaE=0;

  
  double rmsdiff=0.0, maxdiff=0.0;

  double rmsdiffa=0.0, maxdiffa=0.0;
  double rmsdiffb=0.0, maxdiffb=0.0;





    if(sol.Pa.n_rows!=Nbf || sol.Pb.n_rows!=Nbf)

      {
	throw std::runtime_error("No starting guess provided for SCF solver!\n");
      }





  
  if(verbose) {




    if(sol.Ea.n_elem) {
      printf("alpha: ");
      print_E(sol.Ea,occa,false);
      printf("beta:  ");
      print_E(sol.Eb,occb,false);
    }

    fflush(stdout);
  }

//# 226 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"


  
  if(verbose && maxiter>0) {
    fprintf(stderr,"Running ");




    fprintf(stderr,"restricted open-shell ");
//# 242 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"


    fprintf(stderr,"HF ");
//# 252 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

    fprintf(stderr,"calculation");
    if(densityfit) {
      fprintf(stderr," with density fitting");
    }

//# 274 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
    fprintf(stderr,".\n%4s %16s %10s %9s %9s %9s %10s\n","iter","E","dE","RMS dens","MAX dens","DIIS err","titer (s)");
  }
  fflush(stdout);

//# 302 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

  
  bool converged=false;

  



  uscf_t oldsol;

  oldsol.en.E=0.0;

  

  std::vector<double> occar(occa), occbr(occb);
  while(occar.size()<Sinvh.n_cols)
    occar.push_back(0.0);
  while(occbr.size()<Sinvh.n_cols)
    occbr.push_back(0.0);






  
  
  if(maxiter>0) {
    chkptp->open();
    chkptp->write("tol",intthr);
    chkptp->write("P",sol.P);
    chkptp->write(sol.en);

    chkptp->write("Ca",sol.Ca);
    chkptp->write("Cb",sol.Cb);

    chkptp->write("Ea",sol.Ea);
    chkptp->write("Eb",sol.Eb);

    chkptp->write("Pa",sol.Pa);
    chkptp->write("Pb",sol.Pb);

    chkptp->write("occa",occar);
    chkptp->write("occb",occbr);

    
    chkptp->write("Restricted",0);
//# 356 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
    chkptp->write("Converged",converged);
    chkptp->close();

  } else {
    Timer tg;
    
    Fock_ROHF(sol,occa,occb); nfock++;;
    
    if(verbose) {
      if(shift==0.0)
	printf("Solving SCF equations ... ");
      else
	printf("Solving SCF equations with level shift %.3f ... ",shift);
      fflush(stdout);
      t.set();
    }

    
    dipmom=dipole_moment(sol.P,*basisp);
    
    diagonalize(sol,shift);

    if(verbose)
      printf("done (%s)\n",t.elapsed().c_str());
  }

  
  int iiter=1;

//# 391 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

  while(iiter<=maxiter) {
    Timer titer;

    if(verbose)
      printf("\n ******* Iteration %4i ********\n",iiter);

    
    Fock_ROHF(sol,occa,occb); nfock++;;

    
    deltaE=sol.en.E-oldsol.en.E;





    
    diis.update(sol.Ha,sol.Hb,sol.Pa,sol.Pb,sol.en.E,diiserr);


    if(iiter>1 && usebroyden) {





      
      arma::mat Hs=sol.Ha+sol.Hb;
      arma::mat Hd=sol.Ha-sol.Hb;

      arma::mat Hsold=oldsol.Ha+oldsol.Hb;
      arma::mat Hdold=oldsol.Ha-oldsol.Hb;

      
      broyd_sum.push_x(MatToVec(Hsold));
      broyd_sum.push_f(MatToVec(Hsold-Hs));

      broyd_diff.push_x(MatToVec(Hdold));
      broyd_diff.push_f(MatToVec(Hdold-Hd));

    }

    
    try {




      
      diis.solve_F(sol.Ha,sol.Hb);

      diissucc=true;
    } catch(std::runtime_error &) {
      diissucc=false;
    }

    
    if(usebroyden && !diissucc && iiter>1) {

      if(verbose) {
	printf("Performing Broyden interpolation of Fock operator ... ");
	fflush(stdout);
	t.set();
      }





      arma::mat Hs=VecToMat(broyd_sum.update_x(),Nbf,Nbf);
      arma::mat Hd=VecToMat(broyd_diff.update_x(),Nbf,Nbf);

      
      sol.Ha=0.5*(Hs+Hd);
      sol.Hb=0.5*(Hs-Hd);


      if(verbose)
	printf("done (%s)\n",t.elapsed().c_str());
    }

    
    oldsol=sol;

    if(usetrrh) {
      
      if(verbose) {
	printf("\nSolving TRRH equations.\n");
	fflush(stdout);
	t.set();
      }

//# 496 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
      arma::mat Canew;
      arma::vec Eanew;
      TRRH_update(sol.Ha,sol.Ca,S,Canew,Eanew,nocca,verbose,trrhmins);

      arma::mat Cbnew;
      arma::vec Ebnew;
      TRRH_update(sol.Hb,sol.Cb,S,Cbnew,Ebnew,noccb,verbose,trrhmins);

      
      sol.Ca=Canew;
      sol.Cb=Cbnew;

      sol.Ea=Eanew;
      sol.Eb=Ebnew;

      
      check_orth(sol.Ca,S,false);
      check_orth(sol.Cb,S,false);


      if(verbose)
	printf("TRRH solved in %s.\n\n",t.elapsed().c_str());

    } else {
      
      if(verbose) {
	if(shift==0.0)
	  printf("\nSolving SCF equations ... ");
	else
	  printf("\nSolving SCF equations with level shift %.3f ... ",shift);
	fflush(stdout);
	t.set();
      }

      
      diagonalize(sol,shift);

      if(verbose)
	printf("done (%s)\n",t.elapsed().c_str());
    }

    
    form_density(sol,occa,occb);;

    
    arma::mat deltaP=sol.P-oldsol.P;

    arma::mat deltaPa=sol.Pa-oldsol.Pa;
    arma::mat deltaPb=sol.Pb-oldsol.Pb;


    
    dipmom=dipole_moment(sol.P,*basisp);

    
    maxdiff=max_abs(deltaP/2.0);
    rmsdiff=rms_norm(deltaP/2.0);

    
    double maxdiff_cvd(maxdiff);
    double rmsdiff_cvd(maxdiff);


    maxdiffa=max_abs(deltaPa);
    maxdiffb=max_abs(deltaPb);

    rmsdiffa=rms_norm(deltaPa);
    rmsdiffb=rms_norm(deltaPb);

    maxdiff_cvd=std::max(maxdiffa,maxdiffb);
    rmsdiff_cvd=std::max(rmsdiffa,rmsdiffb);


    
    if(verbose) {
      printf("\n");
      printf("%-30s: % .16e\n","Total energy",sol.en.E);
      printf("%-30s: % e\n","DIIS error",diiserr);
      printf("%-30s: % e\n","Energy change",deltaE);
      printf("%-30s: % e\n","Max total density change",maxdiff);
      printf("%-30s: % e\n","Max rms   density change",rmsdiff);

      printf("%-30s: % e\n","Max total alpha density change",maxdiffa);
      printf("%-30s: % e\n","Max rms   alpha density change",rmsdiffa);
      printf("%-30s: % e\n","Max total beta  density change",maxdiffb);
      printf("%-30s: % e\n","Max rms   beta  density change",rmsdiffb);

      printf("Dipole mu = (% 08.8f, % 08.8f, % 08.8f) D\n",dipmom(0)/AUINDEBYE,dipmom(1)/AUINDEBYE,dipmom(2)/AUINDEBYE);

      printf("\nIteration took %s.\n",titer.elapsed().c_str());
      fflush(stdout);
    }

    
    if(verbose)
      fprintf(stderr,"%4i % 16.8f % 10.3e %9.3e %9.3e %9.3e %10.3f\n",iiter,sol.en.E,deltaE,rmsdiff_cvd,maxdiff_cvd,diiserr,titer.get());

//# 605 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
    if(diiserr < convthr) {
      converged=true;
      if(verbose)
        printf("\n ******* Convergence achieved ********\n");
    }

    
    
    chkptp->open();
    chkptp->write("P",sol.P);
    chkptp->write(sol.en);

    chkptp->write("Ca",sol.Ca);
    chkptp->write("Cb",sol.Cb);

    chkptp->write("Ea",sol.Ea);
    chkptp->write("Eb",sol.Eb);

    chkptp->write("Ha",sol.Ha);
    chkptp->write("Hb",sol.Hb);

    chkptp->write("Pa",sol.Pa);
    chkptp->write("Pb",sol.Pb);

    chkptp->write("occa",occar);
    chkptp->write("occb",occbr);

    
    chkptp->write("Restricted",0);
//# 642 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
    chkptp->write("Converged",converged);
    chkptp->close();

    
    if(converged)
      break;

    iiter++;
  } 

  if(verbose) {

    if(converged) {
      std::string method=



	"RO"
//# 666 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
	"HF"

	;

      printf("Solution of %s took %s.\n",method.c_str(),ttot.elapsed().c_str());
      fprintf(stderr,"%s converged in %s.\n",method.c_str(),ttot.elapsed().c_str());
    }

    
    printf("\n");
    printf("%-21s energy: % .16e\n","Kinetic",sol.en.Ekin);
    printf("%-21s energy: % .16e\n","Nuclear attraction",sol.en.Enuca);
    printf("%-21s energy: % .16e\n","Total one-electron",sol.en.Eone);
    printf("%-21s energy: % .16e\n","Nuclear repulsion",sol.en.Enucr);
    printf("%-21s energy: % .16e\n","Coulomb",sol.en.Ecoul);

    printf("%-21s energy: % .16e\n","Exchange",sol.en.Exc);




    printf("-----------------------------------------------------\n");
    printf("%28s: % .16e\n","Total energy",sol.en.E);
    printf("%28s: % .16e\n","Virial factor",-sol.en.E/sol.en.Ekin);

    printf("\nDipole mu = (% 08.8f, % 08.8f, % 08.8f) D\n",dipmom(0)/AUINDEBYE,dipmom(1)/AUINDEBYE,dipmom(2)/AUINDEBYE);

    printf("\n");
    



    printf("alpha: ");
    print_E(sol.Ea,occa,true);
    printf("beta:  ");
    print_E(sol.Eb,occb,true);



    if(dimcalc) {
      
      arma::ivec mv(basisp->get_m_values());
      arma::umat occ;

//# 718 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
      occ.zeros(arma::max(arma::abs(mv))+1,2);

      
      arma::ivec mca(m_classify(sol.Ca,mv));
      for(size_t i=0;i<nocca;i++)
        occ(std::abs(mca(i)),0)++;
      arma::ivec mcb(m_classify(sol.Cb,mv));
      for(size_t i=0;i<noccb;i++)
        occ(std::abs(mcb(i)),1)++;


      printf("\nOrbital occupations by symmetry\n");
      for(size_t m=0;m<occ.n_rows;m++) {




        if(occ(m,0)+occ(m,1)>0)
          printf("m = %i: %2i %2i\n",(int) m,(int) occ(m,0),(int) occ(m,1));

      }
    }


  }

  if(converged) {
    if(doforce) {
      arma::vec f;
//# 760 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
      ERROR_INFO();
      throw std::runtime_error("Forces not supported for this method.\n");


      chkptp->write("Force",f);
    }

  } else if(maxiter>0) {
    std::ostringstream oss;
    oss << "Error in function " << __FUNCTION__ << " (file " << "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in" << ", near line " << 769 << "): SCF did not converge in "<<maxiter<<" iterations!\n";
    throw std::runtime_error(oss.str());
  }




}
//# 1 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
//# 1 "/usr/include/stdc-predef.h" 1















 










 






 

//# 43 "/usr/include/stdc-predef.h"

//# 51 "/usr/include/stdc-predef.h"



 


 


//# 1 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in" 2














 


void SCF::RDFT(rscf_t & sol, const std::vector<double> & occs, double convthr, const dft_t dft0)

//# 38 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
{

  
//# 47 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

//# 65 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

  




  
//# 79 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

  
//# 93 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"








  



  size_t nocc;
  for(nocc=occs.size()-1;nocc<occs.size();nocc--)
    if(occs[nocc]>0)
      break;
  nocc++;
//# 124 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

  int nfock=0;

  Timer t;
  Timer ttot;

  
  double diiserr=DBL_MAX;
  
  bool diissucc=false;

  
  arma::mat orbs;
  arma::mat Horth;


  
  rDIIS diis(S,Sinvh,usediis,diiseps,diisthr,useadiis,verbose,diisorder);
  
  Broyden broyd(verbose);

  
  arma::mat J(Nbf,Nbf), K(Nbf,Nbf);
  J.zeros();
  K.zeros();
//# 158 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

  
  arma::vec dipmom;

  
  double deltaE=0;

  
  double rmsdiff=0.0, maxdiff=0.0;






  if(sol.P.n_rows!=Nbf)



      {
	throw std::runtime_error("No starting guess provided for SCF solver!\n");
      }





  
  if(verbose) {

    if(sol.E.n_elem)
      print_E(sol.E,occs,false);
//# 198 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
    fflush(stdout);
  }


  
  dft_t dft(dft0);

  
  double omega, kfull, kshort;
  range_separation(dft.x_func,omega,kfull,kshort);

  if(verbose) {
    if(omega!=0.0) {
      printf("\nUsing range separated exchange with range separation constant omega = % .3f.\n",omega);
      printf("Using % .3f %% short range and % .3f %% long range exchange.\n",(kfull+kshort)*100,kfull*100);
    } else if(kfull!=0.0)
      printf("\nUsing hybrid exchange with % .3f %% of exact exchange.\n",kfull*100);
    else
      printf("\nA pure exchange functional used, no exact exchange.\n");
  }

  
  if(omega!=0.0)
    fill_rs(omega);

  DFTGrid grid(basisp,verbose,dft.lobatto);
  DFTGrid nlgrid(basisp,verbose,dft.lobatto);



  
  if(verbose && maxiter>0) {
    fprintf(stderr,"Running ");

    fprintf(stderr,"restricted ");
//# 242 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"




    if(dft.c_func>0) {
      
      fprintf(stderr,"%s-%s ",get_keyword(dft.x_func).c_str(),get_keyword(dft.c_func).c_str());
    } else
      fprintf(stderr,"%s ",get_keyword(dft.x_func).c_str());


    fprintf(stderr,"calculation");
    if(densityfit) {
      fprintf(stderr," with density fitting");
    }

//# 272 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"


    fprintf(stderr,".\n%4s %16s %10s %9s %9s %9s %10s\n","iter","E","dE","RMS dens","MAX dens","DIIS err","titer (s)");
  }
  fflush(stdout);


  if(dft.x_func>0 || dft.c_func>0) {
    if(dft.adaptive) {

      
      grid.construct(sol.P,dft.gridtol,dft.x_func,dft.c_func);




    } else {
      
      grid.construct(dft.nrad,dft.lmax,dft.x_func,dft.c_func,strictint);
      
      if(dft.nl)
	nlgrid.construct(dft.nlnrad,dft.nllmax,true,false,false,strictint,true);
    }

    if(verbose) {
      fflush(stdout);
      fprintf(stderr,"%-65s %10.3f\n","    DFT grid formation",t.get());
    }
  }


  
  bool converged=false;

  

  rscf_t oldsol;



  oldsol.en.E=0.0;

  
//# 322 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
  std::vector<double> occr(occs);
  while(occr.size()<Sinvh.n_cols)
    occr.push_back(0.0);


  
  
  if(maxiter>0) {
    chkptp->open();
    chkptp->write("tol",intthr);
    chkptp->write("P",sol.P);
    chkptp->write(sol.en);
//# 350 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
    chkptp->write("C",sol.C);
    chkptp->write("E",sol.E);
    chkptp->write("occs",occr);
    
    chkptp->write("Restricted",1);

    chkptp->write("Converged",converged);
    chkptp->close();

  } else {
    Timer tg;
    
    Fock_RDFT(sol,occs,dft,grid,nlgrid); nfock++;;
    
    if(verbose) {
      if(shift==0.0)
	printf("Solving SCF equations ... ");
      else
	printf("Solving SCF equations with level shift %.3f ... ",shift);
      fflush(stdout);
      t.set();
    }

    
    dipmom=dipole_moment(sol.P,*basisp);
    
    diagonalize(sol,shift);

    if(verbose)
      printf("done (%s)\n",t.elapsed().c_str());
  }

  
  int iiter=1;


  
  
  if(dft.nl)
    dft.nl=false;


  while(iiter<=maxiter) {
    Timer titer;

    if(verbose)
      printf("\n ******* Iteration %4i ********\n",iiter);

    
    Fock_RDFT(sol,occs,dft,grid,nlgrid); nfock++;;

    
    deltaE=sol.en.E-oldsol.en.E;


    
    diis.update(sol.H,sol.P,sol.en.E,diiserr);





    if(iiter>1 && usebroyden) {

      
      broyd.push_x(MatToVec(oldsol.H));
      broyd.push_f(MatToVec(oldsol.H-sol.H));
//# 432 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
    }

    
    try {

      
      diis.solve_F(sol.H);




      diissucc=true;
    } catch(std::runtime_error &) {
      diissucc=false;
    }

    
    if(usebroyden && !diissucc && iiter>1) {

      if(verbose) {
	printf("Performing Broyden interpolation of Fock operator ... ");
	fflush(stdout);
	t.set();
      }


      
      sol.H=VecToMat(broyd.update_x(),Nbf,Nbf);
//# 468 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

      if(verbose)
	printf("done (%s)\n",t.elapsed().c_str());
    }

    
    oldsol=sol;

    if(usetrrh) {
      
      if(verbose) {
	printf("\nSolving TRRH equations.\n");
	fflush(stdout);
	t.set();
      }


      arma::mat Cnew;
      arma::vec Enew;
      TRRH_update(sol.H,sol.C,S,Cnew,Enew,nocc,verbose,trrhmins);

      
      sol.C=Cnew;
      sol.E=Enew;

      
      check_orth(sol.C,S,false);
//# 515 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

      if(verbose)
	printf("TRRH solved in %s.\n\n",t.elapsed().c_str());

    } else {
      
      if(verbose) {
	if(shift==0.0)
	  printf("\nSolving SCF equations ... ");
	else
	  printf("\nSolving SCF equations with level shift %.3f ... ",shift);
	fflush(stdout);
	t.set();
      }

      
      diagonalize(sol,shift);

      if(verbose)
	printf("done (%s)\n",t.elapsed().c_str());
    }

    
    form_density(sol,occs);;

    
    arma::mat deltaP=sol.P-oldsol.P;





    
    dipmom=dipole_moment(sol.P,*basisp);

    
    maxdiff=max_abs(deltaP/2.0);
    rmsdiff=rms_norm(deltaP/2.0);

    
    double maxdiff_cvd(maxdiff);
    double rmsdiff_cvd(maxdiff);

//# 568 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

    
    if(verbose) {
      printf("\n");
      printf("%-30s: % .16e\n","Total energy",sol.en.E);
      printf("%-30s: % e\n","DIIS error",diiserr);
      printf("%-30s: % e\n","Energy change",deltaE);
      printf("%-30s: % e\n","Max total density change",maxdiff);
      printf("%-30s: % e\n","Max rms   density change",rmsdiff);
//# 583 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
      printf("Dipole mu = (% 08.8f, % 08.8f, % 08.8f) D\n",dipmom(0)/AUINDEBYE,dipmom(1)/AUINDEBYE,dipmom(2)/AUINDEBYE);

      printf("\nIteration took %s.\n",titer.elapsed().c_str());
      fflush(stdout);
    }

    
    if(verbose)
      fprintf(stderr,"%4i % 16.8f % 10.3e %9.3e %9.3e %9.3e %10.3f\n",iiter,sol.en.E,deltaE,rmsdiff_cvd,maxdiff_cvd,diiserr,titer.get());


    if(dft0.nl && !dft.nl && diiserr<=1e-3) {
      if(verbose) {
	printf("Turning on non-local correlation contributions\n");
	fprintf(stderr,"Turning on non-local correlation contributions\n");
      }
      dft.nl=true;
      diis.clear();
      continue;
    }
    else

    if(diiserr < convthr) {
      converged=true;
      if(verbose)
        printf("\n ******* Convergence achieved ********\n");
    }

    
    
    chkptp->open();
    chkptp->write("P",sol.P);
    chkptp->write(sol.en);
//# 635 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
    chkptp->write("C",sol.C);
    chkptp->write("E",sol.E);
    chkptp->write("H",sol.H);
    chkptp->write("occs",occr);
    
    chkptp->write("Restricted",1);

    chkptp->write("Converged",converged);
    chkptp->close();

    
    if(converged)
      break;

    iiter++;
  } 

  if(verbose) {

    if(converged) {
      std::string method=

	"R"
//# 664 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
	"DFT"



	;

      printf("Solution of %s took %s.\n",method.c_str(),ttot.elapsed().c_str());
      fprintf(stderr,"%s converged in %s.\n",method.c_str(),ttot.elapsed().c_str());
    }

    
    printf("\n");
    printf("%-21s energy: % .16e\n","Kinetic",sol.en.Ekin);
    printf("%-21s energy: % .16e\n","Nuclear attraction",sol.en.Enuca);
    printf("%-21s energy: % .16e\n","Total one-electron",sol.en.Eone);
    printf("%-21s energy: % .16e\n","Nuclear repulsion",sol.en.Enucr);
    printf("%-21s energy: % .16e\n","Coulomb",sol.en.Ecoul);



    printf("%-21s energy: % .16e\n","Exchange-correlation",sol.en.Exc);
    printf("%-21s energy: % .16e\n","Non-local correlation",sol.en.Enl);

    printf("-----------------------------------------------------\n");
    printf("%28s: % .16e\n","Total energy",sol.en.E);
    printf("%28s: % .16e\n","Virial factor",-sol.en.E/sol.en.Ekin);

    printf("\nDipole mu = (% 08.8f, % 08.8f, % 08.8f) D\n",dipmom(0)/AUINDEBYE,dipmom(1)/AUINDEBYE,dipmom(2)/AUINDEBYE);

    printf("\n");
    

    print_E(sol.E,occs,true);
//# 703 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"


    if(dimcalc) {
      
      arma::ivec mv(basisp->get_m_values());
      arma::umat occ;


      occ.zeros(arma::max(arma::abs(mv))+1,1);

      
      arma::ivec mc(m_classify(sol.C,basisp->get_m_values()));
      for(size_t i=0;i<nocc;i++)
        occ(std::abs(mc(i)),0)++;
//# 728 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

      printf("\nOrbital occupations by symmetry\n");
      for(size_t m=0;m<occ.n_rows;m++) {

        if(occ(m,0)>0)
          printf("m = %i: %2i\n",(int) m,(int) occ(m,0));




      }
    }


  }

  if(converged) {
    if(doforce) {
      arma::vec f;
//# 754 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
      f=force_RDFT(sol,occs,dft,grid,nlgrid,intthr);
//# 763 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

      chkptp->write("Force",f);
    }

  } else if(maxiter>0) {
    std::ostringstream oss;
    oss << "Error in function " << __FUNCTION__ << " (file " << "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in" << ", near line " << 769 << "): SCF did not converge in "<<maxiter<<" iterations!\n";
    throw std::runtime_error(oss.str());
  }




}
//# 1 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
//# 1 "/usr/include/stdc-predef.h" 1















 










 






 

//# 43 "/usr/include/stdc-predef.h"

//# 51 "/usr/include/stdc-predef.h"



 


 


//# 1 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in" 2














 

//# 24 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
void SCF::UDFT(uscf_t & sol, const std::vector<double> & occa, const std::vector<double> & occb, double convthr, const dft_t dft0)

//# 38 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
{

  
//# 49 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
  

//# 65 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

  




  
//# 79 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"

  




















  


//# 111 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
  size_t nocca;
  for(nocca=occa.size()-1;nocca<occa.size();nocca--)
    if(occa[nocca]>0)
      break;
  nocca++;

  size_t noccb;
  for(noccb=occb.size()-1;noccb<occb.size();noccb--)
    if(occb[noccb]>0)
      break;
  noccb++;



  int nfock=0;

  Timer t;
  Timer ttot;

  
  double diiserr=DBL_MAX;
  
  bool diissucc=false;

  
  arma::mat orbs;
  arma::mat Horth;

//# 150 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
  
  uDIIS diis(S,Sinvh,diiscomb,usediis,diiseps,diisthr,useadiis,verbose,diisorder);

  
  Broyden broyd_sum(verbose);
  Broyden broyd_diff(verbose);



  
  arma::vec dipmom;

  
  double deltaE=0;

  
  double rmsdiff=0.0, maxdiff=0.0;

  double rmsdiffa=0.0, maxdiffa=0.0;
  double rmsdiffb=0.0, maxdiffb=0.0;





    if(sol.Pa.n_rows!=Nbf || sol.Pb.n_rows!=Nbf)

      {
	throw std::runtime_error("No starting guess provided for SCF solver!\n");
      }





  
  if(verbose) {




    if(sol.Ea.n_elem) {
      printf("alpha: ");
      print_E(sol.Ea,occa,false);
      printf("beta:  ");
      print_E(sol.Eb,occb,false);
    }

    fflush(stdout);
  }


  
  dft_t dft(dft0);

  
  double omega, kfull, kshort;
  range_separation(dft.x_func,omega,kfull,kshort);

  if(verbose) {
    if(omega!=0.0) {
      printf("\nUsing range separated exchange with range separation constant omega = % .3f.\n",omega);
      printf("Using % .3f %% short range and % .3f %% long range exchange.\n",(kfull+kshort)*100,kfull*100);
    } else if(kfull!=0.0)
      printf("\nUsing hybrid exchange with % .3f %% of exact exchange.\n",kfull*100);
    else
      printf("\nA pure exchange functional used, no exact exchange.\n");
  }

  
  if(omega!=0.0)
    fill_rs(omega);

  DFTGrid grid(basisp,verbose,dft.lobatto);
  DFTGrid nlgrid(basisp,verbose,dft.lobatto);



  
  if(verbose && maxiter>0) {
    fprintf(stderr,"Running ");
//# 238 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
    fprintf(stderr,"unrestricted ");







    if(dft.c_func>0) {
      
      fprintf(stderr,"%s-%s ",get_keyword(dft.x_func).c_str(),get_keyword(dft.c_func).c_str());
    } else
      fprintf(stderr,"%s ",get_keyword(dft.x_func).c_str());


    fprintf(stderr,"calculation");
    if(densityfit) {
      fprintf(stderr," with density fitting");
    }

//# 272 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"


    fprintf(stderr,".\n%4s %16s %10s %9s %9s %9s %10s\n","iter","E","dE","RMS dens","MAX dens","DIIS err","titer (s)");
  }
  fflush(stdout);


  if(dft.x_func>0 || dft.c_func>0) {
    if(dft.adaptive) {




      
      grid.construct(sol.Pa,sol.Pb,dft.gridtol,dft.x_func,dft.c_func);

    } else {
      
      grid.construct(dft.nrad,dft.lmax,dft.x_func,dft.c_func,strictint);
      
      if(dft.nl)
	nlgrid.construct(dft.nlnrad,dft.nllmax,true,false,false,strictint,true);
    }

    if(verbose) {
      fflush(stdout);
      fprintf(stderr,"%-65s %10.3f\n","    DFT grid formation",t.get());
    }
  }


  
  bool converged=false;

  



  uscf_t oldsol;

  oldsol.en.E=0.0;

  

  std::vector<double> occar(occa), occbr(occb);
  while(occar.size()<Sinvh.n_cols)
    occar.push_back(0.0);
  while(occbr.size()<Sinvh.n_cols)
    occbr.push_back(0.0);






  
  
  if(maxiter>0) {
    chkptp->open();
    chkptp->write("tol",intthr);
    chkptp->write("P",sol.P);
    chkptp->write(sol.en);

    chkptp->write("Ca",sol.Ca);
    chkptp->write("Cb",sol.Cb);

    chkptp->write("Ea",sol.Ea);
    chkptp->write("Eb",sol.Eb);

    chkptp->write("Pa",sol.Pa);
    chkptp->write("Pb",sol.Pb);

    chkptp->write("occa",occar);
    chkptp->write("occb",occbr);

    
    chkptp->write("Restricted",0);
//# 356 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
    chkptp->write("Converged",converged);
    chkptp->close();

  } else {
    Timer tg;
    
    Fock_UDFT(sol,occa,occb,dft,grid,nlgrid); nfock++;;
    
    if(verbose) {
      if(shift==0.0)
	printf("Solving SCF equations ... ");
      else
	printf("Solving SCF equations with level shift %.3f ... ",shift);
      fflush(stdout);
      t.set();
    }

    
    dipmom=dipole_moment(sol.P,*basisp);
    
    diagonalize(sol,shift);

    if(verbose)
      printf("done (%s)\n",t.elapsed().c_str());
  }

  
  int iiter=1;


  
  
  if(dft.nl)
    dft.nl=false;


  while(iiter<=maxiter) {
    Timer titer;

    if(verbose)
      printf("\n ******* Iteration %4i ********\n",iiter);

    
    Fock_UDFT(sol,occa,occb,dft,grid,nlgrid); nfock++;;

    
    deltaE=sol.en.E-oldsol.en.E;





    
    diis.update(sol.Ha,sol.Hb,sol.Pa,sol.Pb,sol.en.E,diiserr);


    if(iiter>1 && usebroyden) {





      
      arma::mat Hs=sol.Ha+sol.Hb;
      arma::mat Hd=sol.Ha-sol.Hb;

      arma::mat Hsold=oldsol.Ha+oldsol.Hb;
      arma::mat Hdold=oldsol.Ha-oldsol.Hb;

      
      broyd_sum.push_x(MatToVec(Hsold));
      broyd_sum.push_f(MatToVec(Hsold-Hs));

      broyd_diff.push_x(MatToVec(Hdold));
      broyd_diff.push_f(MatToVec(Hdold-Hd));

    }

    
    try {




      
      diis.solve_F(sol.Ha,sol.Hb);

      diissucc=true;
    } catch(std::runtime_error &) {
      diissucc=false;
    }

    
    if(usebroyden && !diissucc && iiter>1) {

      if(verbose) {
	printf("Performing Broyden interpolation of Fock operator ... ");
	fflush(stdout);
	t.set();
      }





      arma::mat Hs=VecToMat(broyd_sum.update_x(),Nbf,Nbf);
      arma::mat Hd=VecToMat(broyd_diff.update_x(),Nbf,Nbf);

      
      sol.Ha=0.5*(Hs+Hd);
      sol.Hb=0.5*(Hs-Hd);


      if(verbose)
	printf("done (%s)\n",t.elapsed().c_str());
    }

    
    oldsol=sol;

    if(usetrrh) {
      
      if(verbose) {
	printf("\nSolving TRRH equations.\n");
	fflush(stdout);
	t.set();
      }

//# 496 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
      arma::mat Canew;
      arma::vec Eanew;
      TRRH_update(sol.Ha,sol.Ca,S,Canew,Eanew,nocca,verbose,trrhmins);

      arma::mat Cbnew;
      arma::vec Ebnew;
      TRRH_update(sol.Hb,sol.Cb,S,Cbnew,Ebnew,noccb,verbose,trrhmins);

      
      sol.Ca=Canew;
      sol.Cb=Cbnew;

      sol.Ea=Eanew;
      sol.Eb=Ebnew;

      
      check_orth(sol.Ca,S,false);
      check_orth(sol.Cb,S,false);


      if(verbose)
	printf("TRRH solved in %s.\n\n",t.elapsed().c_str());

    } else {
      
      if(verbose) {
	if(shift==0.0)
	  printf("\nSolving SCF equations ... ");
	else
	  printf("\nSolving SCF equations with level shift %.3f ... ",shift);
	fflush(stdout);
	t.set();
      }

      
      diagonalize(sol,shift);

      if(verbose)
	printf("done (%s)\n",t.elapsed().c_str());
    }

    
    form_density(sol,occa,occb);;

    
    arma::mat deltaP=sol.P-oldsol.P;

    arma::mat deltaPa=sol.Pa-oldsol.Pa;
    arma::mat deltaPb=sol.Pb-oldsol.Pb;


    
    dipmom=dipole_moment(sol.P,*basisp);

    
    maxdiff=max_abs(deltaP/2.0);
    rmsdiff=rms_norm(deltaP/2.0);

    
    double maxdiff_cvd(maxdiff);
    double rmsdiff_cvd(maxdiff);


    maxdiffa=max_abs(deltaPa);
    maxdiffb=max_abs(deltaPb);

    rmsdiffa=rms_norm(deltaPa);
    rmsdiffb=rms_norm(deltaPb);

    maxdiff_cvd=std::max(maxdiffa,maxdiffb);
    rmsdiff_cvd=std::max(rmsdiffa,rmsdiffb);


    
    if(verbose) {
      printf("\n");
      printf("%-30s: % .16e\n","Total energy",sol.en.E);
      printf("%-30s: % e\n","DIIS error",diiserr);
      printf("%-30s: % e\n","Energy change",deltaE);
      printf("%-30s: % e\n","Max total density change",maxdiff);
      printf("%-30s: % e\n","Max rms   density change",rmsdiff);

      printf("%-30s: % e\n","Max total alpha density change",maxdiffa);
      printf("%-30s: % e\n","Max rms   alpha density change",rmsdiffa);
      printf("%-30s: % e\n","Max total beta  density change",maxdiffb);
      printf("%-30s: % e\n","Max rms   beta  density change",rmsdiffb);

      printf("Dipole mu = (% 08.8f, % 08.8f, % 08.8f) D\n",dipmom(0)/AUINDEBYE,dipmom(1)/AUINDEBYE,dipmom(2)/AUINDEBYE);

      printf("\nIteration took %s.\n",titer.elapsed().c_str());
      fflush(stdout);
    }

    
    if(verbose)
      fprintf(stderr,"%4i % 16.8f % 10.3e %9.3e %9.3e %9.3e %10.3f\n",iiter,sol.en.E,deltaE,rmsdiff_cvd,maxdiff_cvd,diiserr,titer.get());


    if(dft0.nl && !dft.nl && diiserr<=1e-3) {
      if(verbose) {
	printf("Turning on non-local correlation contributions\n");
	fprintf(stderr,"Turning on non-local correlation contributions\n");
      }
      dft.nl=true;
      diis.clear();
      continue;
    }
    else

    if(diiserr < convthr) {
      converged=true;
      if(verbose)
        printf("\n ******* Convergence achieved ********\n");
    }

    
    
    chkptp->open();
    chkptp->write("P",sol.P);
    chkptp->write(sol.en);

    chkptp->write("Ca",sol.Ca);
    chkptp->write("Cb",sol.Cb);

    chkptp->write("Ea",sol.Ea);
    chkptp->write("Eb",sol.Eb);

    chkptp->write("Ha",sol.Ha);
    chkptp->write("Hb",sol.Hb);

    chkptp->write("Pa",sol.Pa);
    chkptp->write("Pb",sol.Pb);

    chkptp->write("occa",occar);
    chkptp->write("occb",occbr);

    
    chkptp->write("Restricted",0);
//# 642 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
    chkptp->write("Converged",converged);
    chkptp->close();

    
    if(converged)
      break;

    iiter++;
  } 

  if(verbose) {

    if(converged) {
      std::string method=





	"U"


	"DFT"



	;

      printf("Solution of %s took %s.\n",method.c_str(),ttot.elapsed().c_str());
      fprintf(stderr,"%s converged in %s.\n",method.c_str(),ttot.elapsed().c_str());
    }

    
    printf("\n");
    printf("%-21s energy: % .16e\n","Kinetic",sol.en.Ekin);
    printf("%-21s energy: % .16e\n","Nuclear attraction",sol.en.Enuca);
    printf("%-21s energy: % .16e\n","Total one-electron",sol.en.Eone);
    printf("%-21s energy: % .16e\n","Nuclear repulsion",sol.en.Enucr);
    printf("%-21s energy: % .16e\n","Coulomb",sol.en.Ecoul);



    printf("%-21s energy: % .16e\n","Exchange-correlation",sol.en.Exc);
    printf("%-21s energy: % .16e\n","Non-local correlation",sol.en.Enl);

    printf("-----------------------------------------------------\n");
    printf("%28s: % .16e\n","Total energy",sol.en.E);
    printf("%28s: % .16e\n","Virial factor",-sol.en.E/sol.en.Ekin);

    printf("\nDipole mu = (% 08.8f, % 08.8f, % 08.8f) D\n",dipmom(0)/AUINDEBYE,dipmom(1)/AUINDEBYE,dipmom(2)/AUINDEBYE);

    printf("\n");
    



    printf("alpha: ");
    print_E(sol.Ea,occa,true);
    printf("beta:  ");
    print_E(sol.Eb,occb,true);



    if(dimcalc) {
      
      arma::ivec mv(basisp->get_m_values());
      arma::umat occ;

//# 718 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
      occ.zeros(arma::max(arma::abs(mv))+1,2);

      
      arma::ivec mca(m_classify(sol.Ca,mv));
      for(size_t i=0;i<nocca;i++)
        occ(std::abs(mca(i)),0)++;
      arma::ivec mcb(m_classify(sol.Cb,mv));
      for(size_t i=0;i<noccb;i++)
        occ(std::abs(mcb(i)),1)++;


      printf("\nOrbital occupations by symmetry\n");
      for(size_t m=0;m<occ.n_rows;m++) {




        if(occ(m,0)+occ(m,1)>0)
          printf("m = %i: %2i %2i\n",(int) m,(int) occ(m,0),(int) occ(m,1));

      }
    }


  }

  if(converged) {
    if(doforce) {
      arma::vec f;
//# 757 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in"
      f=force_UDFT(sol,occa,occb,dft,grid,nlgrid,intthr);






      chkptp->write("Force",f);
    }

  } else if(maxiter>0) {
    std::ostringstream oss;
    oss << "Error in function " << __FUNCTION__ << " (file " << "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-solvers.cpp.in" << ", near line " << 769 << "): SCF did not converge in "<<maxiter<<" iterations!\n";
    throw std::runtime_error(oss.str());
  }




}
