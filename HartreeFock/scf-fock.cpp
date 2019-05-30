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
//# 1 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
//# 1 "/usr/include/stdc-predef.h" 1







 





 



 
//# 43 "/usr/include/stdc-predef.h"
//# 51 "/usr/include/stdc-predef.h"

 

 

//# 1 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in" 2







 
//# 1 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/global.h" 1







 






















 





































/*



//# 140 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/global.h"
//# 147 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/global.h"
//# 18 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in" 2


 


void SCF::Fock_RHF(rscf_t & sol, const std::vector<double> & occs) const
//# 44 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
{
  Timer t;
  
  if(sol . P . n_rows != Nbf || sol . P . n_cols != Nbf) { std::ostringstream oss; oss << "sol.P" << " should be " << Nbf << " x " << Nbf << " but is " << sol . P . n_rows << " x " << sol . P . n_cols << "!\n"; throw std::runtime_error(oss . str());};


  
  sol.J.zeros(Nbf,Nbf);
  sol.K.zeros(Nbf,Nbf);


//# 72 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"



  
  arma::mat Kfull, Kfull_im, Kshort, Kshort_im;
  Kfull.zeros(Nbf,Nbf);
  Kfull_im.zeros(Nbf,Nbf);
  Kshort.zeros(Nbf,Nbf);
  Kshort_im.zeros(Nbf,Nbf);
//# 97 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
  
  memset(&sol.en, 0, sizeof(energy_t));
  
  if(densityfit) {
    if(verbose) {
      printf("Forming density fitted Coulomb matrix ... ");
      fflush(stdout);
      t.set();
    }
    sol.J=dfit.calcJ(sol.P);
    if(verbose) {
      printf("done (%s)\n",t.elapsed().c_str());
      fflush(stdout);
    }


      if(verbose) {
	printf("Forming density fitted exchange matrix ... ");
	fflush(stdout);
	t.set();
      }

      if(sol.P_im.n_rows == sol.P.n_rows && sol.P_im.n_cols == sol.P.n_cols) {
	
	arma::cx_mat cK(dfit.calcK(sol.cC,occs,fitmem));
	Kfull=arma::real(cK);
	Kfull_im=arma::imag(cK);
      } else
	Kfull=dfit.calcK(sol.C,occs,fitmem);
//# 150 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
      if(verbose) {
	printf("done (%s)\n",t.elapsed().c_str());
	fflush(stdout);
      }
//# 194 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"



  } else { 
    if(cholesky) {
      if(verbose) {
	printf("Forming Cholesky Coulomb matrix ... ");
	fflush(stdout);
	t.set();
      }
      sol.J=chol.calcJ(sol.P);
      if(verbose) {
	printf("done (%s)\n",t.elapsed().c_str());
	fflush(stdout);
      }


	if(verbose) {
	  printf("Forming Cholesky exchange matrix ... ");
	  fflush(stdout);
	  t.set();
	}

	if(sol.P_im.n_rows == sol.P.n_rows && sol.P_im.n_cols == sol.P.n_cols) {
	  arma::cx_mat cK(chol.calcK(sol.cC,occs));
	  Kfull=arma::real(cK);
	  Kfull_im=arma::imag(cK);
	} else
	  Kfull=chol.calcK(sol.C,occs);
//# 248 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
	if(verbose) {
	  printf("done (%s)\n",t.elapsed().c_str());
	  fflush(stdout);
	}


//# 295 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
    } else {
      if(direct) {


	  if(verbose) {
	    printf("Computing HF Coulomb and exchange matrices.\nScreening integrals with tolerance %.3e ... ",intthr);
	    fflush(stdout);
	    t.set();
	  }
	  if(sol.P_im.n_rows == sol.P.n_rows && sol.P_im.n_cols == sol.P.n_cols) {
	    
	    arma::cx_mat cP(sol.P*std::complex<double>(1.0,0.0) + sol.P_im*std::complex<double>(0.0,1.0));
	    arma::cx_mat cK;
	    if(!decfock) {
	      scr.calcJK(cP,sol.J,cK,intthr);
	    } else {
	      
	      arma::cx_mat Phlp=decconv*cP*arma::trans(decconv);
	      
	      scr.calcJK(Phlp,sol.J,cK,intthr);
	      
	      sol.J=arma::trans(decconv)*sol.J*decconv;
	      cK=arma::trans(decconv)*cK*decconv;
	    }
	    Kfull=arma::real(cK);
	    Kfull_im=arma::imag(cK);
	  } else {
	    if(!decfock) {
	      scr.calcJK(sol.P,sol.J,Kfull,intthr);
	    } else {
	      
	      arma::mat Phlp=decconv*sol.P*arma::trans(decconv);
	      
	      scr.calcJK(Phlp,sol.J,Kfull,intthr);
	      
	      sol.J=arma::trans(decconv)*sol.J*decconv;
	      Kfull=arma::trans(decconv)*Kfull*decconv;
	    }
	  }
	  if(verbose) {
	    printf("done (%s)\n",t.elapsed().c_str());
	    fflush(stdout);
	  }
//# 416 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"

//# 555 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
      } else {
	
	if(verbose) {
	  printf("Forming HF Coulomb matrix ... ");
	  fflush(stdout);
	  t.set();
	}
	sol.J=tab.calcJ(sol.P);
	if(verbose) {
	  printf("done (%s)\n",t.elapsed().c_str());
	  fflush(stdout);
	}


	  if(verbose) {
	    printf("Forming HF exchange matrix ... ");
	    fflush(stdout);
	    t.set();
	  }

	  if(sol.P_im.n_rows == sol.P.n_rows && sol.P_im.n_cols == sol.P.n_cols) {
	    
	    arma::cx_mat cP(sol.P*std::complex<double>(1.0,0.0) + sol.P_im*std::complex<double>(0.0,1.0));
	    arma::cx_mat cK(tab.calcK(cP));
	    Kfull=arma::real(cK);
	    Kfull_im=arma::imag(cK);
	  } else
	    Kfull=tab.calcK(sol.P);
//# 610 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
	  if(verbose) {
	    printf("done (%s)\n",t.elapsed().c_str());
	    fflush(stdout);
	  }


//# 662 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
      }
    }
  }

//# 689 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
  
  sol.K=Kfull;



//# 763 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
  
//# 777 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
  sol.H=Hcore+sol.J-0.5*sol.K;
  sol.en.Exc=-0.25*arma::trace(sol.P*sol.K);
  if(sol.P.n_rows == sol.P_im.n_rows && sol.P.n_cols == sol.P_im.n_cols)
    sol.en.Exc+=0.25*arma::trace(sol.P_im*sol.K_im);

//# 811 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"

//# 822 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
  
  if(freeze.size()>0) {
    freeze_orbs(freeze,sol.C,S,sol.H,verbose);


  }
  fflush(stdout);
  if(dimcalc) {
    
    std::vector<nucleus_t> nuclei(basisp->get_nuclei());
    for(size_t i=0;i<nuclei.size();i++)
      if(nuclei[i].r.x!=0.0 || nuclei[i].r.y!=0.0)
        throw std::logic_error("Nuclei must be on z axis for dimer calculation!\n");
    
    arma::ivec mvals(basisp->get_m_values());
    
    sol.H=block_m(sol.H,mvals);


  }
  
  sol.en.Ekin=arma::trace(sol.P*T);
  sol.en.Enuca=arma::trace(sol.P*Vnuc);
  sol.en.Enucr=Enuc;
  sol.en.Eone=arma::trace(sol.P*Hcore);
  sol.en.Ecoul=0.5*arma::trace(sol.P*sol.J);
  
  sol.en.Eel=sol.en.Ecoul+sol.en.Exc+sol.en.Eone+sol.en.Enl;
  sol.en.E=sol.en.Eel+sol.en.Enucr;
  
  if(!arma::is_finite(sol.H)) {
//# 873 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
    throw std::runtime_error("Fock operator is not finite.\n");
  }
//# 903 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
  if(!std::isfinite(sol.en.E)) {
//# 935 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
    std::ostringstream oss;
    oss << "\nSomething wrong with total energy " << sol.en.E <<"?\nEnding program.\n";
    throw std::runtime_error(oss.str());
  }
}
//# 1 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
//# 1 "/usr/include/stdc-predef.h" 1







 





 



 
//# 43 "/usr/include/stdc-predef.h"
//# 51 "/usr/include/stdc-predef.h"

 

 

//# 1 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in" 2







 
//# 1 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/global.h" 1







 






















 














































//# 140 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/global.h"
//# 147 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/global.h"
//# 18 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in" 2


 
//# 33 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
void SCF::Fock_UHF(uscf_t & sol, const std::vector<double> & occa, const std::vector<double> & occb) const
//# 44 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
{
  Timer t;
  
  if(sol . P . n_rows != Nbf || sol . P . n_cols != Nbf) { std::ostringstream oss; oss << "sol.P" << " should be " << Nbf << " x " << Nbf << " but is " << sol . P . n_rows << " x " << sol . P . n_cols << "!\n"; throw std::runtime_error(oss . str());};
  if(sol . Pa . n_rows != Nbf || sol . Pa . n_cols != Nbf) { std::ostringstream oss; oss << "sol.Pa" << " should be " << Nbf << " x " << Nbf << " but is " << sol . Pa . n_rows << " x " << sol . Pa . n_cols << "!\n"; throw std::runtime_error(oss . str());};
  if(sol . Pb . n_rows != Nbf || sol . Pb . n_cols != Nbf) { std::ostringstream oss; oss << "sol.Pb" << " should be " << Nbf << " x " << Nbf << " but is " << sol . Pb . n_rows << " x " << sol . Pb . n_cols << "!\n"; throw std::runtime_error(oss . str());};

  
  sol.J.zeros(Nbf,Nbf);

  sol.Ka.zeros(Nbf,Nbf);
  sol.Kb.zeros(Nbf,Nbf);

//# 72 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"



  
//# 87 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
  arma::mat Kafull, Kafull_im, Kbfull, Kbfull_im, Kashort, Kashort_im, Kbshort, Kbshort_im;
  Kafull.zeros(Nbf,Nbf);
  Kafull_im.zeros(Nbf,Nbf);
  Kbfull.zeros(Nbf,Nbf);
  Kbfull_im.zeros(Nbf,Nbf);
  Kashort.zeros(Nbf,Nbf);
  Kashort_im.zeros(Nbf,Nbf);
  Kbshort.zeros(Nbf,Nbf);
  Kbshort_im.zeros(Nbf,Nbf);

  
  memset(&sol.en, 0, sizeof(energy_t));
  
  if(densityfit) {
    if(verbose) {
      printf("Forming density fitted Coulomb matrix ... ");
      fflush(stdout);
      t.set();
    }
    sol.J=dfit.calcJ(sol.P);
    if(verbose) {
      printf("done (%s)\n",t.elapsed().c_str());
      fflush(stdout);
    }


      if(verbose) {
	printf("Forming density fitted exchange matrix ... ");
	fflush(stdout);
	t.set();
      }
//# 134 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
      if(sol.Pa_im.n_rows == sol.Pa.n_rows && sol.Pa_im.n_cols == sol.Pa.n_cols) {
	
	arma::cx_mat cKa(dfit.calcK(sol.cCa,occa,fitmem));
	Kafull=arma::real(cKa);
	Kafull_im=arma::imag(cKa);
      } else
	Kafull=dfit.calcK(sol.Ca,occa,fitmem);
      if(sol.Pb_im.n_rows == sol.Pb.n_rows && sol.Pb_im.n_cols == sol.Pb.n_cols) {
	
	arma::cx_mat cKb(dfit.calcK(sol.cCb,occb,fitmem));
	Kbfull=arma::real(cKb);
	Kbfull_im=arma::imag(cKb);
      } else
	Kbfull=dfit.calcK(sol.Cb,occb,fitmem);

      if(verbose) {
	printf("done (%s)\n",t.elapsed().c_str());
	fflush(stdout);
      }
//# 194 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"



  } else { 
    if(cholesky) {
      if(verbose) {
	printf("Forming Cholesky Coulomb matrix ... ");
	fflush(stdout);
	t.set();
      }
      sol.J=chol.calcJ(sol.P);
      if(verbose) {
	printf("done (%s)\n",t.elapsed().c_str());
	fflush(stdout);
      }


	if(verbose) {
	  printf("Forming Cholesky exchange matrix ... ");
	  fflush(stdout);
	  t.set();
	}
//# 234 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
	if(sol.Pa_im.n_rows == sol.Pa.n_rows && sol.Pa_im.n_cols == sol.Pa.n_cols) {
	  arma::cx_mat cKa(chol.calcK(sol.cCa,occa));
	  Kafull=arma::real(cKa);
	  Kafull_im=arma::imag(cKa);
	} else
	  Kafull=chol.calcK(sol.Ca,occa);
	if(sol.Pb_im.n_rows == sol.Pb.n_rows && sol.Pb_im.n_cols == sol.Pb.n_cols) {
	  arma::cx_mat cKb(chol.calcK(sol.cCb,occb));
	  Kbfull=arma::real(cKb);
	  Kbfull_im=arma::imag(cKb);
	} else
	  Kbfull=chol.calcK(sol.Cb,occb);

	if(verbose) {
	  printf("done (%s)\n",t.elapsed().c_str());
	  fflush(stdout);
	}


//# 295 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
    } else {
      if(direct) {
//# 419 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
	

	  if(verbose) {
	    printf("Computing HF Coulomb and exchange matrices.\nScreening integrals with tolerance %.3e ... ",intthr);
	    fflush(stdout);
	    t.set();
	  }
	  if(sol.Pa_im.n_rows == sol.Pa.n_rows && sol.Pa_im.n_cols == sol.Pa.n_cols) {
	    
	    arma::cx_mat cPa(sol.Pa*std::complex<double>(1.0,0.0) + sol.Pa_im*std::complex<double>(0.0,1.0));
	    arma::cx_mat cPb(sol.Pb*std::complex<double>(1.0,0.0) + sol.Pb_im*std::complex<double>(0.0,1.0));
	    arma::cx_mat cKa, cKb;
	    if(!decfock) {
	      scr.calcJK(cPa,cPb,sol.J,cKa,cKb,intthr);
	    } else {
	      
	      arma::cx_mat Pahlp=decconv*cPa*arma::trans(decconv);
	      arma::cx_mat Pbhlp=decconv*cPb*arma::trans(decconv);
	      
	      scr.calcJK(Pahlp,Pbhlp,sol.J,cKa,cKb,intthr);
	      
	      sol.J=arma::trans(decconv)*sol.J*decconv;
	      cKa=arma::trans(decconv)*cKa*decconv;
	      cKb=arma::trans(decconv)*cKb*decconv;
	    }
	    Kafull=arma::real(cKa);
	    Kafull_im=arma::imag(cKa);
	    Kbfull=arma::real(cKb);
	    Kbfull_im=arma::imag(cKb);
	  } else {
	    if(!decfock) {
	      scr.calcJK(sol.Pa,sol.Pb,sol.J,Kafull,Kbfull,intthr);
	    } else {
	      
	      arma::mat Pahlp=decconv*sol.Pa*arma::trans(decconv);
	      arma::mat Pbhlp=decconv*sol.Pb*arma::trans(decconv);
	      
	      scr.calcJK(Pahlp,Pbhlp,sol.J,Kafull,Kbfull,intthr);
	      
	      sol.J=arma::trans(decconv)*sol.J*decconv;
	      Kafull=arma::trans(decconv)*Kafull*decconv;
	      Kbfull=arma::trans(decconv)*Kbfull*decconv;
	    }
	  }
	  if(verbose) {
	    printf("done (%s)\n",t.elapsed().c_str());
	    fflush(stdout);
	  }

//# 553 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"

      } else {
	
	if(verbose) {
	  printf("Forming HF Coulomb matrix ... ");
	  fflush(stdout);
	  t.set();
	}
	sol.J=tab.calcJ(sol.P);
	if(verbose) {
	  printf("done (%s)\n",t.elapsed().c_str());
	  fflush(stdout);
	}


	  if(verbose) {
	    printf("Forming HF exchange matrix ... ");
	    fflush(stdout);
	    t.set();
	  }
//# 592 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
	  if(sol.Pa_im.n_rows == sol.Pa.n_rows && sol.Pa_im.n_cols == sol.Pa.n_cols) {
	    
	    arma::cx_mat cPa(sol.Pa*std::complex<double>(1.0,0.0) + sol.Pa_im*std::complex<double>(0.0,1.0));
	    arma::cx_mat cKa(tab.calcK(cPa));
	    Kafull=arma::real(cKa);
	    Kafull_im=arma::imag(cKa);
	  } else
	    Kafull=tab.calcK(sol.Pa);
	  if(sol.Pb_im.n_rows == sol.Pb.n_rows && sol.Pb_im.n_cols == sol.Pb.n_cols) {
	    
	    arma::cx_mat cPb(sol.Pb*std::complex<double>(1.0,0.0) + sol.Pb_im*std::complex<double>(0.0,1.0));
	    arma::cx_mat cKb(tab.calcK(cPb));
	    Kbfull=arma::real(cKb);
	    Kbfull_im=arma::imag(cKb);
	  } else
	    Kbfull=tab.calcK(sol.Pb);

	  if(verbose) {
	    printf("done (%s)\n",t.elapsed().c_str());
	    fflush(stdout);
	  }


//# 662 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
      }
    }
  }

//# 689 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
  

  sol.Ka=Kafull;
  sol.Kb=Kbfull;

//# 763 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
  
//# 785 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
//# 797 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
  sol.Ha=Hcore+sol.J-sol.Ka;
  sol.Hb=Hcore+sol.J-sol.Kb;
  sol.en.Exc=-0.5*(arma::trace(sol.Pa*sol.Ka)+arma::trace(sol.Pb*sol.Kb));
  if(sol.Pa.n_rows == sol.Pa_im.n_rows && sol.Pa.n_cols == sol.Pa_im.n_cols)
    sol.en.Exc+=0.5*(arma::trace(sol.Pa_im*sol.Ka_im)+arma::trace(sol.Pb_im*sol.Kb_im));





//# 822 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
  
  if(freeze.size()>0) {

    freeze_orbs(freeze,sol.Ca,S,sol.Ha,verbose);
    freeze_orbs(freeze,sol.Cb,S,sol.Hb,verbose);
  }
  fflush(stdout);
  if(dimcalc) {
    
    std::vector<nucleus_t> nuclei(basisp->get_nuclei());
    for(size_t i=0;i<nuclei.size();i++)
      if(nuclei[i].r.x!=0.0 || nuclei[i].r.y!=0.0)
        throw std::logic_error("Nuclei must be on z axis for dimer calculation!\n");
    
    arma::ivec mvals(basisp->get_m_values());
    

    sol.Ha=block_m(sol.Ha,mvals);
    sol.Hb=block_m(sol.Hb,mvals);
  }
  
  sol.en.Ekin=arma::trace(sol.P*T);
  sol.en.Enuca=arma::trace(sol.P*Vnuc);
  sol.en.Enucr=Enuc;
  sol.en.Eone=arma::trace(sol.P*Hcore);
  sol.en.Ecoul=0.5*arma::trace(sol.P*sol.J);
  
  sol.en.Eel=sol.en.Ecoul+sol.en.Exc+sol.en.Eone+sol.en.Enl;
  sol.en.E=sol.en.Eel+sol.en.Enucr;
  
//# 876 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
  if(!arma::is_finite(sol.Ha)) {
//# 887 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
    throw std::runtime_error("Alpha Fock operator is not finite.\n");
  }
  if(!arma::is_finite(sol.Hb)) {
//# 900 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
    throw std::runtime_error("Beta Fock operator is not finite.\n");
  }

  if(!std::isfinite(sol.en.E)) {
//# 935 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
    std::ostringstream oss;
    oss << "\nSomething wrong with total energy " << sol.en.E <<"?\nEnding program.\n";
    throw std::runtime_error(oss.str());
  }
}
//# 1 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
//# 1 "/usr/include/stdc-predef.h" 1







 





 



 
//# 43 "/usr/include/stdc-predef.h"
//# 51 "/usr/include/stdc-predef.h"

 

 

//# 1 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in" 2







 
//# 1 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/global.h" 1









 






















 









































//# 140 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/global.h"
//# 147 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/global.h"
//# 18 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in" 2


 
//# 36 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
void SCF::Fock_ROHF(uscf_t & sol, const std::vector<double> & occa, const std::vector<double> & occb) const
//# 44 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
{
  Timer t;
  
  if(sol . P . n_rows != Nbf || sol . P . n_cols != Nbf) { std::ostringstream oss; oss << "sol.P" << " should be " << Nbf << " x " << Nbf << " but is " << sol . P . n_rows << " x " << sol . P . n_cols << "!\n"; throw std::runtime_error(oss . str());};
  if(sol . Pa . n_rows != Nbf || sol . Pa . n_cols != Nbf) { std::ostringstream oss; oss << "sol.Pa" << " should be " << Nbf << " x " << Nbf << " but is " << sol . Pa . n_rows << " x " << sol . Pa . n_cols << "!\n"; throw std::runtime_error(oss . str());};
  if(sol . Pb . n_rows != Nbf || sol . Pb . n_cols != Nbf) { std::ostringstream oss; oss << "sol.Pb" << " should be " << Nbf << " x " << Nbf << " but is " << sol . Pb . n_rows << " x " << sol . Pb . n_cols << "!\n"; throw std::runtime_error(oss . str());};

  
  sol.J.zeros(Nbf,Nbf);

  sol.Ka.zeros(Nbf,Nbf);
  sol.Kb.zeros(Nbf,Nbf);

//# 72 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"



  
//# 87 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
  arma::mat Kafull, Kafull_im, Kbfull, Kbfull_im, Kashort, Kashort_im, Kbshort, Kbshort_im;
  Kafull.zeros(Nbf,Nbf);
  Kafull_im.zeros(Nbf,Nbf);
  Kbfull.zeros(Nbf,Nbf);
  Kbfull_im.zeros(Nbf,Nbf);
  Kashort.zeros(Nbf,Nbf);
  Kashort_im.zeros(Nbf,Nbf);
  Kbshort.zeros(Nbf,Nbf);
  Kbshort_im.zeros(Nbf,Nbf);

  
  memset(&sol.en, 0, sizeof(energy_t));
  
  if(densityfit) {
    if(verbose) {
      printf("Forming density fitted Coulomb matrix ... ");
      fflush(stdout);
      t.set();
    }
    sol.J=dfit.calcJ(sol.P);
    if(verbose) {
      printf("done (%s)\n",t.elapsed().c_str());
      fflush(stdout);
    }


      if(verbose) {
	printf("Forming density fitted exchange matrix ... ");
	fflush(stdout);
	t.set();
      }
//# 134 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
      if(sol.Pa_im.n_rows == sol.Pa.n_rows && sol.Pa_im.n_cols == sol.Pa.n_cols) {
	
	arma::cx_mat cKa(dfit.calcK(sol.cCa,occa,fitmem));
	Kafull=arma::real(cKa);
	Kafull_im=arma::imag(cKa);
      } else
	Kafull=dfit.calcK(sol.Ca,occa,fitmem);
      if(sol.Pb_im.n_rows == sol.Pb.n_rows && sol.Pb_im.n_cols == sol.Pb.n_cols) {
	
	arma::cx_mat cKb(dfit.calcK(sol.cCb,occb,fitmem));
	Kbfull=arma::real(cKb);
	Kbfull_im=arma::imag(cKb);
      } else
	Kbfull=dfit.calcK(sol.Cb,occb,fitmem);

      if(verbose) {
	printf("done (%s)\n",t.elapsed().c_str());
	fflush(stdout);
      }
//# 194 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"



  } else { 
    if(cholesky) {
      if(verbose) {
	printf("Forming Cholesky Coulomb matrix ... ");
	fflush(stdout);
	t.set();
      }
      sol.J=chol.calcJ(sol.P);
      if(verbose) {
	printf("done (%s)\n",t.elapsed().c_str());
	fflush(stdout);
      }


	if(verbose) {
	  printf("Forming Cholesky exchange matrix ... ");
	  fflush(stdout);
	  t.set();
	}
//# 234 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
	if(sol.Pa_im.n_rows == sol.Pa.n_rows && sol.Pa_im.n_cols == sol.Pa.n_cols) {
	  arma::cx_mat cKa(chol.calcK(sol.cCa,occa));
	  Kafull=arma::real(cKa);
	  Kafull_im=arma::imag(cKa);
	} else
	  Kafull=chol.calcK(sol.Ca,occa);
	if(sol.Pb_im.n_rows == sol.Pb.n_rows && sol.Pb_im.n_cols == sol.Pb.n_cols) {
	  arma::cx_mat cKb(chol.calcK(sol.cCb,occb));
	  Kbfull=arma::real(cKb);
	  Kbfull_im=arma::imag(cKb);
	} else
	  Kbfull=chol.calcK(sol.Cb,occb);

	if(verbose) {
	  printf("done (%s)\n",t.elapsed().c_str());
	  fflush(stdout);
	}


//# 295 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
    } else {
      if(direct) {
//# 419 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
	

	  if(verbose) {
	    printf("Computing HF Coulomb and exchange matrices.\nScreening integrals with tolerance %.3e ... ",intthr);
	    fflush(stdout);
	    t.set();
	  }
	  if(sol.Pa_im.n_rows == sol.Pa.n_rows && sol.Pa_im.n_cols == sol.Pa.n_cols) {
	    
	    arma::cx_mat cPa(sol.Pa*std::complex<double>(1.0,0.0) + sol.Pa_im*std::complex<double>(0.0,1.0));
	    arma::cx_mat cPb(sol.Pb*std::complex<double>(1.0,0.0) + sol.Pb_im*std::complex<double>(0.0,1.0));
	    arma::cx_mat cKa, cKb;
	    if(!decfock) {
	      scr.calcJK(cPa,cPb,sol.J,cKa,cKb,intthr);
	    } else {
	      
	      arma::cx_mat Pahlp=decconv*cPa*arma::trans(decconv);
	      arma::cx_mat Pbhlp=decconv*cPb*arma::trans(decconv);
	      
	      scr.calcJK(Pahlp,Pbhlp,sol.J,cKa,cKb,intthr);
	      
	      sol.J=arma::trans(decconv)*sol.J*decconv;
	      cKa=arma::trans(decconv)*cKa*decconv;
	      cKb=arma::trans(decconv)*cKb*decconv;
	    }
	    Kafull=arma::real(cKa);
	    Kafull_im=arma::imag(cKa);
	    Kbfull=arma::real(cKb);
	    Kbfull_im=arma::imag(cKb);
	  } else {
	    if(!decfock) {
	      scr.calcJK(sol.Pa,sol.Pb,sol.J,Kafull,Kbfull,intthr);
	    } else {
	      
	      arma::mat Pahlp=decconv*sol.Pa*arma::trans(decconv);
	      arma::mat Pbhlp=decconv*sol.Pb*arma::trans(decconv);
	      
	      scr.calcJK(Pahlp,Pbhlp,sol.J,Kafull,Kbfull,intthr);
	      
	      sol.J=arma::trans(decconv)*sol.J*decconv;
	      Kafull=arma::trans(decconv)*Kafull*decconv;
	      Kbfull=arma::trans(decconv)*Kbfull*decconv;
	    }
	  }
	  if(verbose) {
	    printf("done (%s)\n",t.elapsed().c_str());
	    fflush(stdout);
	  }

//# 553 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"

      } else {
	
	if(verbose) {
	  printf("Forming HF Coulomb matrix ... ");
	  fflush(stdout);
	  t.set();
	}
	sol.J=tab.calcJ(sol.P);
	if(verbose) {
	  printf("done (%s)\n",t.elapsed().c_str());
	  fflush(stdout);
	}


	  if(verbose) {
	    printf("Forming HF exchange matrix ... ");
	    fflush(stdout);
	    t.set();
	  }
//# 592 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
	  if(sol.Pa_im.n_rows == sol.Pa.n_rows && sol.Pa_im.n_cols == sol.Pa.n_cols) {
	    
	    arma::cx_mat cPa(sol.Pa*std::complex<double>(1.0,0.0) + sol.Pa_im*std::complex<double>(0.0,1.0));
	    arma::cx_mat cKa(tab.calcK(cPa));
	    Kafull=arma::real(cKa);
	    Kafull_im=arma::imag(cKa);
	  } else
	    Kafull=tab.calcK(sol.Pa);
	  if(sol.Pb_im.n_rows == sol.Pb.n_rows && sol.Pb_im.n_cols == sol.Pb.n_cols) {
	    
	    arma::cx_mat cPb(sol.Pb*std::complex<double>(1.0,0.0) + sol.Pb_im*std::complex<double>(0.0,1.0));
	    arma::cx_mat cKb(tab.calcK(cPb));
	    Kbfull=arma::real(cKb);
	    Kbfull_im=arma::imag(cKb);
	  } else
	    Kbfull=tab.calcK(sol.Pb);

	  if(verbose) {
	    printf("done (%s)\n",t.elapsed().c_str());
	    fflush(stdout);
	  }


//# 662 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
      }
    }
  }

//# 689 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
  

  sol.Ka=Kafull;
  sol.Kb=Kbfull;

//# 763 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
  
//# 785 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
//# 797 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
  sol.Ha=Hcore+sol.J-sol.Ka;
  sol.Hb=Hcore+sol.J-sol.Kb;
  sol.en.Exc=-0.5*(arma::trace(sol.Pa*sol.Ka)+arma::trace(sol.Pb*sol.Kb));
  if(sol.Pa.n_rows == sol.Pa_im.n_rows && sol.Pa.n_cols == sol.Pa_im.n_cols)
    sol.en.Exc+=0.5*(arma::trace(sol.Pa_im*sol.Ka_im)+arma::trace(sol.Pb_im*sol.Kb_im));

  
  ROHF_update(sol.Ha,sol.Hb,sol.P,S,occa,occb,verbose);



//# 822 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
  
  if(freeze.size()>0) {

    freeze_orbs(freeze,sol.Ca,S,sol.Ha,verbose);
    freeze_orbs(freeze,sol.Cb,S,sol.Hb,verbose);
  }
  fflush(stdout);
  if(dimcalc) {
    
    std::vector<nucleus_t> nuclei(basisp->get_nuclei());
    for(size_t i=0;i<nuclei.size();i++)
      if(nuclei[i].r.x!=0.0 || nuclei[i].r.y!=0.0)
        throw std::logic_error("Nuclei must be on z axis for dimer calculation!\n");
    
    arma::ivec mvals(basisp->get_m_values());
    

    sol.Ha=block_m(sol.Ha,mvals);
    sol.Hb=block_m(sol.Hb,mvals);
  }
  
  sol.en.Ekin=arma::trace(sol.P*T);
  sol.en.Enuca=arma::trace(sol.P*Vnuc);
  sol.en.Enucr=Enuc;
  sol.en.Eone=arma::trace(sol.P*Hcore);
  sol.en.Ecoul=0.5*arma::trace(sol.P*sol.J);
  
  sol.en.Eel=sol.en.Ecoul+sol.en.Exc+sol.en.Eone+sol.en.Enl;
  sol.en.E=sol.en.Eel+sol.en.Enucr;
  
//# 876 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
  if(!arma::is_finite(sol.Ha)) {
//# 887 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
    throw std::runtime_error("Alpha Fock operator is not finite.\n");
  }
  if(!arma::is_finite(sol.Hb)) {
//# 900 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
    throw std::runtime_error("Beta Fock operator is not finite.\n");
  }

  if(!std::isfinite(sol.en.E)) {
//# 935 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
    std::ostringstream oss;
    oss << "\nSomething wrong with total energy " << sol.en.E <<"?\nEnding program.\n";
    throw std::runtime_error(oss.str());
  }
}
//# 1 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
//# 1 "/usr/include/stdc-predef.h" 1







 





 



 
//# 43 "/usr/include/stdc-predef.h"
//# 51 "/usr/include/stdc-predef.h"

 

 

//# 1 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in" 2







 
//# 1 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/global.h" 1







 



























 



































//# 140 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/global.h"
//# 147 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/global.h"
//# 18 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in" 2


 
void SCF::Fock_RDFT(rscf_t & sol, const std::vector<double> & occs, const dft_t dft, DFTGrid & grid, DFTGrid & nlgrid) const
//# 44 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
{
  Timer t;
  
  if(sol . P . n_rows != Nbf || sol . P . n_cols != Nbf) { std::ostringstream oss; oss << "sol.P" << " should be " << Nbf << " x " << Nbf << " but is " << sol . P . n_rows << " x " << sol . P . n_cols << "!\n"; throw std::runtime_error(oss . str());};


  
  sol.J.zeros(Nbf,Nbf);
  sol.K.zeros(Nbf,Nbf);



  
  sol.XC.zeros(Nbf,Nbf);
//# 72 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"

  
  double omega, kfull, kshort;
  range_separation(dft.x_func,omega,kfull,kshort);

  
  arma::mat Kfull, Kfull_im, Kshort, Kshort_im;
  Kfull.zeros(Nbf,Nbf);
  Kfull_im.zeros(Nbf,Nbf);
  Kshort.zeros(Nbf,Nbf);
  Kshort_im.zeros(Nbf,Nbf);
//# 97 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
  
  memset(&sol.en, 0, sizeof(energy_t));
  
  if(densityfit) {
    if(verbose) {
      printf("Forming density fitted Coulomb matrix ... ");
      fflush(stdout);
      t.set();
    }
    sol.J=dfit.calcJ(sol.P);
    if(verbose) {
      printf("done (%s)\n",t.elapsed().c_str());
      fflush(stdout);
    }

    if(kfull!=0.0) {
      if(verbose) {
	printf("Forming density fitted exchange matrix ... ");
	fflush(stdout);
	t.set();
      }

      if(sol.P_im.n_rows == sol.P.n_rows && sol.P_im.n_cols == sol.P.n_cols) {
	
	arma::cx_mat cK(dfit.calcK(sol.cC,occs,fitmem));
	Kfull=arma::real(cK);
	Kfull_im=arma::imag(cK);
      } else
	Kfull=dfit.calcK(sol.C,occs,fitmem);
//# 150 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
      if(verbose) {
	printf("done (%s)\n",t.elapsed().c_str());
	fflush(stdout);
      }

      if(kshort!=0.0) {
	if(verbose) {
	  printf("Forming density fitted short-range exchange matrix ... ");
	  fflush(stdout);
	  t.set();
	}

	if(sol.P_im.n_rows == sol.P.n_rows && sol.P_im.n_cols == sol.P.n_cols) {
	  arma::cx_mat cK(dfit_rs.calcK(sol.cC,occs,fitmem));
	  Kshort=arma::real(cK);
	  Kshort_im=arma::imag(cK);
	} else
	  Kshort=dfit_rs.calcK(sol.C,occs,fitmem);
//# 187 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
	if(verbose) {
	  printf("done (%s)\n",t.elapsed().c_str());
	  fflush(stdout);
	}
      }


    }

  } else { 
    if(cholesky) {
      if(verbose) {
	printf("Forming Cholesky Coulomb matrix ... ");
	fflush(stdout);
	t.set();
      }
      sol.J=chol.calcJ(sol.P);
      if(verbose) {
	printf("done (%s)\n",t.elapsed().c_str());
	fflush(stdout);
      }

      if(kfull!=0.0) {

	if(verbose) {
	  printf("Forming Cholesky exchange matrix ... ");
	  fflush(stdout);
	  t.set();
	}

	if(sol.P_im.n_rows == sol.P.n_rows && sol.P_im.n_cols == sol.P.n_cols) {
	  arma::cx_mat cK(chol.calcK(sol.cC,occs));
	  Kfull=arma::real(cK);
	  Kfull_im=arma::imag(cK);
	} else
	  Kfull=chol.calcK(sol.C,occs);
//# 248 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
	if(verbose) {
	  printf("done (%s)\n",t.elapsed().c_str());
	  fflush(stdout);
	}
      }

      if(kshort!=0.0) {
	if(verbose) {
	  printf("Forming short-range Cholesky exchange matrix ... ");
	  fflush(stdout);
	  t.set();
	}

	if(sol.P_im.n_rows == sol.P.n_rows && sol.P_im.n_cols == sol.P.n_cols) {
	  arma::cx_mat cK(chol_rs.calcK(sol.cC,occs));
	  Kshort=arma::real(cK);
	  Kshort_im=arma::imag(cK);
	} else
	  Kshort=chol_rs.calcK(sol.C,occs);
//# 288 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
	if(verbose) {
	  printf("done (%s)\n",t.elapsed().c_str());
	  fflush(stdout);
	}
      }

    } else {
      if(direct) {

	if(kfull!=0.0) {
	  if(verbose) {
	    printf("Computing HF Coulomb and exchange matrices.\nScreening integrals with tolerance %.3e ... ",intthr);
	    fflush(stdout);
	    t.set();
	  }
	  if(sol.P_im.n_rows == sol.P.n_rows && sol.P_im.n_cols == sol.P.n_cols) {
	    
	    arma::cx_mat cP(sol.P*std::complex<double>(1.0,0.0) + sol.P_im*std::complex<double>(0.0,1.0));
	    arma::cx_mat cK;
	    if(!decfock) {
	      scr.calcJK(cP,sol.J,cK,intthr);
	    } else {
	      
	      arma::cx_mat Phlp=decconv*cP*arma::trans(decconv);
	      
	      scr.calcJK(Phlp,sol.J,cK,intthr);
	      
	      sol.J=arma::trans(decconv)*sol.J*decconv;
	      cK=arma::trans(decconv)*cK*decconv;
	    }
	    Kfull=arma::real(cK);
	    Kfull_im=arma::imag(cK);
	  } else {
	    if(!decfock) {
	      scr.calcJK(sol.P,sol.J,Kfull,intthr);
	    } else {
	      
	      arma::mat Phlp=decconv*sol.P*arma::trans(decconv);
	      
	      scr.calcJK(Phlp,sol.J,Kfull,intthr);
	      
	      sol.J=arma::trans(decconv)*sol.J*decconv;
	      Kfull=arma::trans(decconv)*Kfull*decconv;
	    }
	  }
	  if(verbose) {
	    printf("done (%s)\n",t.elapsed().c_str());
	    fflush(stdout);
	  }

	} else {
	  if(verbose) {
	    printf("Computing HF Coulomb matrix.\nScreening integrals with tolerance %.3e ... ",intthr);
	    fflush(stdout);
	    t.set();
	  }
	  if(!decfock)
	    sol.J=scr.calcJ(sol.P,intthr);
	  else {
	    arma::mat Phlp=decconv*sol.P*arma::trans(decconv);
	    arma::mat Jhlp=scr.calcJ(Phlp,intthr);
	    sol.J=arma::trans(decconv)*Jhlp*decconv;
	  }
	  if(verbose) {
	    printf("done (%s)\n",t.elapsed().c_str());
	    fflush(stdout);
	  }
	}
	if(kshort!=0.0) {
	  if(verbose) {
	    printf("Computing screened exchange matrix.\nScreening integrals with tolerance %.3e ... ",intthr);
	    fflush(stdout);
	    t.set();
	  }
	  if(sol.P_im.n_rows == sol.P.n_rows && sol.P_im.n_cols == sol.P.n_cols) {
	    
	    arma::cx_mat cP(sol.P*std::complex<double>(1.0,0.0) + sol.P_im*std::complex<double>(0.0,1.0));
	    arma::cx_mat cK;
	    if(!decfock) {
	      cK=scr_rs.calcK(cP,intthr);
	    } else {
	      
	      arma::cx_mat Phlp=decconv*cP*arma::trans(decconv);
	      
	      cK=scr_rs.calcK(Phlp,intthr);
	      
	      cK=arma::trans(decconv)*cK*decconv;
	    }
	    Kshort=arma::real(cK);
	    Kshort_im=arma::imag(cK);
	  } else {
	    if(!decfock) {
	      Kshort=scr_rs.calcK(sol.P,intthr);
	    } else {
	      
	      arma::mat Phlp=decconv*sol.P*arma::trans(decconv);
	      
	      Kshort=scr_rs.calcK(Phlp,intthr);
	      
	      Kshort=arma::trans(decconv)*Kshort*decconv;
	    }
	  }
	  if(verbose) {
	    printf("done (%s)\n",t.elapsed().c_str());
	    fflush(stdout);
	  }
	}

//# 555 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
      } else {
	
	if(verbose) {
	  printf("Forming HF Coulomb matrix ... ");
	  fflush(stdout);
	  t.set();
	}
	sol.J=tab.calcJ(sol.P);
	if(verbose) {
	  printf("done (%s)\n",t.elapsed().c_str());
	  fflush(stdout);
	}

	if(kfull!=0.0) {

	  if(verbose) {
	    printf("Forming HF exchange matrix ... ");
	    fflush(stdout);
	    t.set();
	  }

	  if(sol.P_im.n_rows == sol.P.n_rows && sol.P_im.n_cols == sol.P.n_cols) {
	    
	    arma::cx_mat cP(sol.P*std::complex<double>(1.0,0.0) + sol.P_im*std::complex<double>(0.0,1.0));
	    arma::cx_mat cK(tab.calcK(cP));
	    Kfull=arma::real(cK);
	    Kfull_im=arma::imag(cK);
	  } else
	    Kfull=tab.calcK(sol.P);
//# 610 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
	  if(verbose) {
	    printf("done (%s)\n",t.elapsed().c_str());
	    fflush(stdout);
	  }
	}

	if(kshort!=0.0) {
	  if(verbose) {
	    printf("Forming short-range HF exchange matrix ... ");
	    fflush(stdout);
	    t.set();
	  }

	  if(sol.P_im.n_rows == sol.P.n_rows && sol.P_im.n_cols == sol.P.n_cols) {
	    
	    arma::cx_mat cP(sol.P*std::complex<double>(1.0,0.0) + sol.P_im*std::complex<double>(0.0,1.0));
	    arma::cx_mat cK(tab_rs.calcK(cP));
	    Kshort=arma::real(cK);
	    Kshort_im=arma::imag(cK);
	  } else
	    Kshort=tab_rs.calcK(sol.P);
//# 656 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
	  if(verbose) {
	    printf("done (%s)\n",t.elapsed().c_str());
	    fflush(stdout);
	  }
	}
      }
    }
  }

  
  sol.K=kfull*Kfull + kshort*Kshort;
  if(sol.P.n_rows == sol.P_im.n_rows && sol.P.n_cols == sol.P_im.n_cols)
    sol.K_im=kfull*Kfull_im + kshort*Kshort_im;
  else
    sol.K_im.clear();
//# 697 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"

  sol.en.Exc=0.0;
  if(dft.x_func>0 || dft.c_func>0) {
    if(verbose) {
      printf("Computing DFT exchange and correlation ... ");
      fflush(stdout);
      t.set();
    }
    double Nelnum; 
    grid.eval_Fxc(dft.x_func,dft.c_func,sol.P,sol.XC,sol.en.Exc,Nelnum);




    double rel_diff=(Nelnum-Nel)*100.0/Nel;

    if(verbose) {
      printf("done (%s)\n",t.elapsed().c_str());
      printf("Numerically integrated density is %.5f (%+.4f %%).\n",Nelnum,rel_diff);
    }
    if(fabs(rel_diff)>1e-2) {
      std::ostringstream oss;
      
      oss << "Warning - numerically integrated density seems inaccurate.\n";
      if(verbose)
	std::cout << oss.str();
      
    }
  }
  
  sol.en.Esic=0.0;
  
  if(dft.nl) {
    if(verbose) {
      printf("Computing non-local correlation ... ");
      fflush(stdout);
      t.set();
    }

    grid.eval_VV10(nlgrid,dft.vv10_b,dft.vv10_C,sol.P,sol.XC,sol.en.Enl);
//# 757 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
    if(verbose) {
      printf("done (%s)\n",t.elapsed().c_str());
    }
  }

  

  sol.H=Hcore+sol.J+sol.XC;
  
  if(kfull!=0.0 || kshort!=0.0) {
    sol.H-=0.5*sol.K;
    sol.en.Exc-=0.25*arma::trace(sol.P*sol.K);
    if(sol.P.n_rows == sol.P_im.n_rows && sol.P.n_cols == sol.P_im.n_cols)
      sol.en.Exc+=0.25*arma::trace(sol.P_im*sol.K_im);
  }
//# 782 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"

//# 811 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"

//# 822 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
  
  if(freeze.size()>0) {
    freeze_orbs(freeze,sol.C,S,sol.H,verbose);


  }
  fflush(stdout);
  if(dimcalc) {
    
    std::vector<nucleus_t> nuclei(basisp->get_nuclei());
    for(size_t i=0;i<nuclei.size();i++)
      if(nuclei[i].r.x!=0.0 || nuclei[i].r.y!=0.0)
        throw std::logic_error("Nuclei must be on z axis for dimer calculation!\n");
    
    arma::ivec mvals(basisp->get_m_values());
    
    sol.H=block_m(sol.H,mvals);


  }
  
  sol.en.Ekin=arma::trace(sol.P*T);
  sol.en.Enuca=arma::trace(sol.P*Vnuc);
  sol.en.Enucr=Enuc;
  sol.en.Eone=arma::trace(sol.P*Hcore);
  sol.en.Ecoul=0.5*arma::trace(sol.P*sol.J);
  
  sol.en.Eel=sol.en.Ecoul+sol.en.Exc+sol.en.Eone+sol.en.Enl;
  sol.en.E=sol.en.Eel+sol.en.Enucr;
  
  if(!arma::is_finite(sol.H)) {
//# 873 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
    throw std::runtime_error("Fock operator is not finite.\n");
  }
//# 903 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
  if(!std::isfinite(sol.en.E)) {
//# 935 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
    std::ostringstream oss;
    oss << "\nSomething wrong with total energy " << sol.en.E <<"?\nEnding program.\n";
    throw std::runtime_error(oss.str());
  }
}
//# 1 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
//# 1 "/usr/include/stdc-predef.h" 1


*/




 





 



 
//# 43 "/usr/include/stdc-predef.h"
//# 51 "/usr/include/stdc-predef.h"

 

 

//# 1 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in" 2







 
//# 1 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/global.h" 1








 






















 






































//# 140 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/global.h"
//# 147 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/global.h"
//# 18 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in" 2



//# 30 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
void SCF::Fock_UDFT(uscf_t & sol, const std::vector<double> & occa, const std::vector<double> & occb, const dft_t dft, DFTGrid & grid, DFTGrid & nlgrid) const
//# 44 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
{

	std::cout << "\nENTERING SCF::Fock_UDFT\n";
  Timer t;
  
  if(sol . P . n_rows != Nbf || sol . P . n_cols != Nbf) { std::ostringstream oss; oss << "sol.P" << " should be " << Nbf << " x " << Nbf << " but is " << sol . P . n_rows << " x " << sol . P . n_cols << "!\n"; throw std::runtime_error(oss . str());};
  if(sol . Pa . n_rows != Nbf || sol . Pa . n_cols != Nbf) { std::ostringstream oss; oss << "sol.Pa" << " should be " << Nbf << " x " << Nbf << " but is " << sol . Pa . n_rows << " x " << sol . Pa . n_cols << "!\n"; throw std::runtime_error(oss . str());};
  if(sol . Pb . n_rows != Nbf || sol . Pb . n_cols != Nbf) { std::ostringstream oss; oss << "sol.Pb" << " should be " << Nbf << " x " << Nbf << " but is " << sol . Pb . n_rows << " x " << sol . Pb . n_cols << "!\n"; throw std::runtime_error(oss . str());};

  
  sol.J.zeros(Nbf,Nbf);

  sol.Ka.zeros(Nbf,Nbf);
  sol.Kb.zeros(Nbf,Nbf);



  
  sol.XCa.zeros(Nbf,Nbf);
  sol.XCb.zeros(Nbf,Nbf);

  
  double omega, kfull, kshort;
  std::cout << "\nBEFORE range_separation(dft.x_func,omega,kfull,kshort);\n";
  range_separation(dft.x_func,omega,kfull,kshort);
  std::cout << "\nAFTER range_separation(dft.x_func,omega,kfull,kshort);\n";
  
//# 87 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
  //arma::mat Kafull, Kafull_im, Kbfull, Kbfull_im, Kashort, Kashort_im, Kbshort, Kbshort_im;
  std::cout << "\nAFTER arma::mat Kafull, Kafull_im, Kbfull, Kbfull_im, Kashort, Kashort_im, Kbshort, Kbshort_im;\n";

  printf("Nbf = %i\n",Nbf);

  /*
  Kafull.zeros(Nbf,Nbf);
  Kafull_im.zeros(Nbf,Nbf);
  Kbfull.zeros(Nbf,Nbf);
  Kbfull_im.zeros(Nbf,Nbf);
  Kashort.zeros(Nbf,Nbf);
  Kashort_im.zeros(Nbf,Nbf);
  Kbshort.zeros(Nbf,Nbf);
  Kbshort_im.zeros(Nbf,Nbf);
  */

  std::cout << "\nAFTER Kbshort_im.zeros(Nbf,Nbf);\n";


  printf("sizeof(energy_t) = %i\n",sizeof(energy_t));
  
  memset(&sol.en, 0, sizeof(energy_t));
  

  std::cout << "\nAFTER memset(&sol.en, 0, sizeof(energy_t));\n";

  /* ESTE TROZO DE update_density PARA ACTUALIZAR rho ME LO HE INVENTADO Y ESTÁ MAL PORQUE ES DE AngularGrid y no de DFGrid.
  grid.upda update_density(sol.Pa,sol.Pb);
  rho.print();
  */

  /*
  if(densityfit) {
    if(verbose) {
      printf("Forming density fitted Coulomb matrix ... ");
      fflush(stdout);
      t.set();
    }
    sol.J=dfit.calcJ(sol.P);
    if(verbose) {
      printf("done (%s)\n",t.elapsed().c_str());
      fflush(stdout);
    }

    if(kfull!=0.0) {
      if(verbose) {
	printf("Forming density fitted exchange matrix ... ");
	fflush(stdout);
	t.set();
      }
//# 134 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
      if(sol.Pa_im.n_rows == sol.Pa.n_rows && sol.Pa_im.n_cols == sol.Pa.n_cols) {
	
	arma::cx_mat cKa(dfit.calcK(sol.cCa,occa,fitmem));
	Kafull=arma::real(cKa);
	Kafull_im=arma::imag(cKa);
      } else
	Kafull=dfit.calcK(sol.Ca,occa,fitmem);
      if(sol.Pb_im.n_rows == sol.Pb.n_rows && sol.Pb_im.n_cols == sol.Pb.n_cols) {
	
	arma::cx_mat cKb(dfit.calcK(sol.cCb,occb,fitmem));
	Kbfull=arma::real(cKb);
	Kbfull_im=arma::imag(cKb);
      } else
	Kbfull=dfit.calcK(sol.Cb,occb,fitmem);

      if(verbose) {
	printf("done (%s)\n",t.elapsed().c_str());
	fflush(stdout);
      }

      if(kshort!=0.0) {
	if(verbose) {
	  printf("Forming density fitted short-range exchange matrix ... ");
	  fflush(stdout);
	  t.set();
	}
//# 173 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
	if(sol.Pa_im.n_rows == sol.Pa.n_rows && sol.Pa_im.n_cols == sol.Pa.n_cols) {
	  arma::cx_mat cKa(dfit_rs.calcK(sol.cCa,occa,fitmem));
	  Kashort=arma::real(cKa);
	  Kashort_im=arma::imag(cKa);
	} else
	  Kashort=dfit_rs.calcK(sol.Ca,occa,fitmem);
	if(sol.Pb_im.n_rows == sol.Pb.n_rows && sol.Pb_im.n_cols == sol.Pb.n_cols) {
	  arma::cx_mat cKb(dfit_rs.calcK(sol.cCb,occb,fitmem));
	  Kbshort=arma::real(cKb);
	  Kbshort_im=arma::imag(cKb);
	} else
	  Kbshort=dfit_rs.calcK(sol.Cb,occb,fitmem);

	if(verbose) {
	  printf("done (%s)\n",t.elapsed().c_str());
	  fflush(stdout);
	}
      }


    }

  } else { 
    if(cholesky) {
      if(verbose) {
	printf("Forming Cholesky Coulomb matrix ... ");
	fflush(stdout);
	t.set();
      }
      sol.J=chol.calcJ(sol.P);
      if(verbose) {
	printf("done (%s)\n",t.elapsed().c_str());
	fflush(stdout);
      }

      if(kfull!=0.0) {

	if(verbose) {
	  printf("Forming Cholesky exchange matrix ... ");
	  fflush(stdout);
	  t.set();
	}
//# 234 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
	if(sol.Pa_im.n_rows == sol.Pa.n_rows && sol.Pa_im.n_cols == sol.Pa.n_cols) {
	  arma::cx_mat cKa(chol.calcK(sol.cCa,occa));
	  Kafull=arma::real(cKa);
	  Kafull_im=arma::imag(cKa);
	} else
	  Kafull=chol.calcK(sol.Ca,occa);
	if(sol.Pb_im.n_rows == sol.Pb.n_rows && sol.Pb_im.n_cols == sol.Pb.n_cols) {
	  arma::cx_mat cKb(chol.calcK(sol.cCb,occb));
	  Kbfull=arma::real(cKb);
	  Kbfull_im=arma::imag(cKb);
	} else
	  Kbfull=chol.calcK(sol.Cb,occb);

	if(verbose) {
	  printf("done (%s)\n",t.elapsed().c_str());
	  fflush(stdout);
	}
      }

      if(kshort!=0.0) {
	if(verbose) {
	  printf("Forming short-range Cholesky exchange matrix ... ");
	  fflush(stdout);
	  t.set();
	}
//# 274 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
	if(sol.Pa_im.n_rows == sol.Pa.n_rows && sol.Pa_im.n_cols == sol.Pa.n_cols) {
	  arma::cx_mat cKa(chol_rs.calcK(sol.cCa,occa));
	  Kashort=arma::real(cKa);
	  Kashort_im=arma::imag(cKa);
	} else
	  Kashort=chol_rs.calcK(sol.Ca,occa);
	if(sol.Pb_im.n_rows == sol.Pb.n_rows && sol.Pb_im.n_cols == sol.Pb.n_cols) {
	  arma::cx_mat cKb(chol_rs.calcK(sol.cCb,occb));
	  Kbshort=arma::real(cKb);
	  Kbshort_im=arma::imag(cKb);
	} else
	  Kbshort=chol_rs.calcK(sol.Cb,occb);

	if(verbose) {
	  printf("done (%s)\n",t.elapsed().c_str());
	  fflush(stdout);
	}
      }

    } else {
      if(direct) {
//# 419 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
	
	if(kfull!=0.0) {
	  if(verbose) {
	    printf("Computing HF Coulomb and exchange matrices.\nScreening integrals with tolerance %.3e ... ",intthr);
	    fflush(stdout);
	    t.set();
	  }
	  if(sol.Pa_im.n_rows == sol.Pa.n_rows && sol.Pa_im.n_cols == sol.Pa.n_cols) {
	    
	    arma::cx_mat cPa(sol.Pa*std::complex<double>(1.0,0.0) + sol.Pa_im*std::complex<double>(0.0,1.0));
	    arma::cx_mat cPb(sol.Pb*std::complex<double>(1.0,0.0) + sol.Pb_im*std::complex<double>(0.0,1.0));
	    arma::cx_mat cKa, cKb;
	    if(!decfock) {
	      scr.calcJK(cPa,cPb,sol.J,cKa,cKb,intthr);
	    } else {
	      
	      arma::cx_mat Pahlp=decconv*cPa*arma::trans(decconv);
	      arma::cx_mat Pbhlp=decconv*cPb*arma::trans(decconv);
	      
	      scr.calcJK(Pahlp,Pbhlp,sol.J,cKa,cKb,intthr);
	      
	      sol.J=arma::trans(decconv)*sol.J*decconv;
	      cKa=arma::trans(decconv)*cKa*decconv;
	      cKb=arma::trans(decconv)*cKb*decconv;
	    }
	    Kafull=arma::real(cKa);
	    Kafull_im=arma::imag(cKa);
	    Kbfull=arma::real(cKb);
	    Kbfull_im=arma::imag(cKb);
	  } else {
	    if(!decfock) {
	      scr.calcJK(sol.Pa,sol.Pb,sol.J,Kafull,Kbfull,intthr);
	    } else {
	      
	      arma::mat Pahlp=decconv*sol.Pa*arma::trans(decconv);
	      arma::mat Pbhlp=decconv*sol.Pb*arma::trans(decconv);
	      
	      scr.calcJK(Pahlp,Pbhlp,sol.J,Kafull,Kbfull,intthr);
	      
	      sol.J=arma::trans(decconv)*sol.J*decconv;
	      Kafull=arma::trans(decconv)*Kafull*decconv;
	      Kbfull=arma::trans(decconv)*Kbfull*decconv;
	    }
	  }
	  if(verbose) {
	    printf("done (%s)\n",t.elapsed().c_str());
	    fflush(stdout);
	  }

	} else {
	  if(verbose) {
	    printf("Computing HF Coulomb matrix.\nScreening integrals with tolerance %.3e ... ",intthr);
	    fflush(stdout);
	    t.set();
	  }
	  if(!decfock)
	    sol.J=scr.calcJ(sol.P,intthr);
	  else {
	    arma::mat Phlp=decconv*sol.P*arma::trans(decconv);
	    arma::mat Jhlp=scr.calcJ(Phlp,intthr);
	    sol.J=arma::trans(decconv)*Jhlp*decconv;
	  }
	  if(verbose) {
	    printf("done (%s)\n",t.elapsed().c_str());
	    fflush(stdout);
	  }
	}

	if(kshort!=0.0) {
	  if(verbose) {
	    printf("Computing screened exchange matrix.\nScreening integrals with tolerance %.3e ... ",intthr);
	    fflush(stdout);
	    t.set();
	  }
	  if(sol.Pa_im.n_rows == sol.Pa.n_rows && sol.Pa_im.n_cols == sol.Pa.n_cols) {
	    
	    arma::cx_mat cPa(sol.Pa*std::complex<double>(1.0,0.0) + sol.Pa_im*std::complex<double>(0.0,1.0));
	    arma::cx_mat cPb(sol.Pb*std::complex<double>(1.0,0.0) + sol.Pb_im*std::complex<double>(0.0,1.0));
	    arma::cx_mat cKa, cKb;
	    if(!decfock) {
	      scr_rs.calcK(cPa,cPb,cKa,cKb,intthr);
	    } else {
	      
	      arma::cx_mat Pahlp=decconv*cPa*arma::trans(decconv);
	      arma::cx_mat Pbhlp=decconv*cPb*arma::trans(decconv);
	      
	      scr_rs.calcK(Pahlp,Pbhlp,cKa,cKb,intthr);
	      
	      cKa=arma::trans(decconv)*cKa*decconv;
	      cKb=arma::trans(decconv)*cKb*decconv;
	    }
	    Kashort=arma::real(cKa);
	    Kashort_im=arma::imag(cKa);
	    Kbshort=arma::real(cKb);
	    Kbshort_im=arma::imag(cKb);
	  } else {
	    if(!decfock) {
	      scr_rs.calcK(sol.Pa,sol.Pb,Kashort,Kbshort,intthr);
	    } else {
	      
	      arma::mat Pahlp=decconv*sol.Pa*arma::trans(decconv);
	      arma::mat Pbhlp=decconv*sol.Pb*arma::trans(decconv);
	      
	      scr_rs.calcK(Pahlp,Pbhlp,Kashort,Kbshort,intthr);
	      
	      Kashort=arma::trans(decconv)*Kashort*decconv;
	      Kbshort=arma::trans(decconv)*Kbshort*decconv;
	    }
	  }
	  if(verbose) {
	    printf("done (%s)\n",t.elapsed().c_str());
	    fflush(stdout);
	  }
	}

      } else {
	
	if(verbose) {
	  printf("Forming HF Coulomb matrix ... ");
	  fflush(stdout);
	  t.set();
	}
	sol.J=tab.calcJ(sol.P);
	if(verbose) {
	  printf("done (%s)\n",t.elapsed().c_str());
	  fflush(stdout);
	}

	if(kfull!=0.0) {

	  if(verbose) {
	    printf("Forming HF exchange matrix ... ");
	    fflush(stdout);
	    t.set();
	  }
//# 592 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
	  if(sol.Pa_im.n_rows == sol.Pa.n_rows && sol.Pa_im.n_cols == sol.Pa.n_cols) {
	    
	    arma::cx_mat cPa(sol.Pa*std::complex<double>(1.0,0.0) + sol.Pa_im*std::complex<double>(0.0,1.0));
	    arma::cx_mat cKa(tab.calcK(cPa));
	    Kafull=arma::real(cKa);
	    Kafull_im=arma::imag(cKa);
	  } else
	    Kafull=tab.calcK(sol.Pa);
	  if(sol.Pb_im.n_rows == sol.Pb.n_rows && sol.Pb_im.n_cols == sol.Pb.n_cols) {
	    
	    arma::cx_mat cPb(sol.Pb*std::complex<double>(1.0,0.0) + sol.Pb_im*std::complex<double>(0.0,1.0));
	    arma::cx_mat cKb(tab.calcK(cPb));
	    Kbfull=arma::real(cKb);
	    Kbfull_im=arma::imag(cKb);
	  } else
	    Kbfull=tab.calcK(sol.Pb);

	  if(verbose) {
	    printf("done (%s)\n",t.elapsed().c_str());
	    fflush(stdout);
	  }
	}

	if(kshort!=0.0) {
	  if(verbose) {
	    printf("Forming short-range HF exchange matrix ... ");
	    fflush(stdout);
	    t.set();
	  }
//# 638 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
	  if(sol.Pa_im.n_rows == sol.Pa.n_rows && sol.Pa_im.n_cols == sol.Pa.n_cols) {
	    
	    arma::cx_mat cPa(sol.Pa*std::complex<double>(1.0,0.0) + sol.Pa_im*std::complex<double>(0.0,1.0));
	    arma::cx_mat cKa(tab_rs.calcK(cPa));
	    Kashort=arma::real(cKa);
	    Kashort_im=arma::imag(cKa);
	  } else
	    Kashort=tab_rs.calcK(sol.Pa);
	  if(sol.Pb_im.n_rows == sol.Pb.n_rows && sol.Pb_im.n_cols == sol.Pb.n_cols) {
	    
	    arma::cx_mat cPb(sol.Pb*std::complex<double>(1.0,0.0) + sol.Pb_im*std::complex<double>(0.0,1.0));
	    arma::cx_mat cKb(tab_rs.calcK(cPb));
	    Kbshort=arma::real(cKb);
	    Kbshort_im=arma::imag(cKb);
	  } else
	    Kbshort=tab_rs.calcK(sol.Pb);
	  if(verbose) {
	    printf("done (%s)\n",t.elapsed().c_str());
	    fflush(stdout);
	  }
	}
      }
    }
  }

  
//# 676 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
  sol.Ka=kfull*Kafull + kshort*Kashort;
  if(sol.Pa.n_rows == sol.Pa_im.n_rows && sol.Pa.n_cols == sol.Pa_im.n_cols)
    sol.Ka_im=kfull*Kafull_im + kshort*Kashort_im;
  else
    sol.Ka_im.clear();
  sol.Kb=kfull*Kbfull + kshort*Kbshort;
  if(sol.Pb.n_rows == sol.Pb_im.n_rows && sol.Pb.n_cols == sol.Pb_im.n_cols)
    sol.Kb_im=kfull*Kbfull_im + kshort*Kbshort_im;
  else
    sol.Kb_im.clear();
//# 697 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"


 */


  std::cout << "\nBEFORE sol.en.Exc=0.0;\n";

  sol.en.Exc=0.0;
  if(dft.x_func>0 || dft.c_func>0) {
    if(verbose) {
      printf("Computing DFT exchange and correlation ... ");
      fflush(stdout);
      t.set();
      std::cout << "\nAFTER t.set;\n";
    }
    double Nelnum;

    std::cout << "\nBEFORE double nelgrid = grid.compute_Nel(sol.Pa,sol.Pb);\n";
    //grid.compute_Nel(const arma::mat & P);
    //double nelgrid = grid.compute_Nel(sol.Pa,sol.Pb);
    //printf("\nnumber of electrons from grid: nelgrid = %.5f\n",nelgrid);
    std::cout << "\nAFTER double nelgrid = grid.compute_Nel(sol.Pa,sol.Pb);\n";

    grid.eval_Fxc(dft.x_func,dft.c_func,sol.Pa,sol.Pb,sol.XCa,sol.XCb,sol.en.Exc,Nelnum);



    double rel_diff=(Nelnum-Nel)*100.0/Nel;

    if(verbose) {
      printf("done (%s)\n",t.elapsed().c_str());
      printf("Numerically integrated density is %.5f (%+.4f %%).\n",Nelnum,rel_diff);
    }
    if(fabs(rel_diff)>1e-2) {
      std::ostringstream oss;
      
      oss << "Warning - numerically integrated density seems inaccurate.\n";
      if(verbose)
	std::cout << oss.str();
      
    }
    std::cout << "\nAFTER if(fabs(rel_diff)>1e-2)\n";
  }
  
  sol.en.Esic=0.0;
  
  if(dft.nl) {
    if(verbose) {
      printf("Computing non-local correlation ... ");
      fflush(stdout);
      t.set();
    }


    arma::mat XC(sol.XCa);
    XC.zeros();
    grid.eval_VV10(nlgrid,dft.vv10_b,dft.vv10_C,sol.P,XC,sol.en.Enl);
    sol.XCa+=XC;
    sol.XCb+=XC;

    if(verbose) {
      printf("done (%s)\n",t.elapsed().c_str());
    }
  }

  
//# 785 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"

  sol.Ha=Hcore+sol.J+sol.XCa;
  sol.Hb=Hcore+sol.J+sol.XCb;
  if(kshort!=0.0 || kfull!=0.0) {
    sol.Ha-=sol.Ka;
    sol.Hb-=sol.Kb;
    sol.en.Exc-=0.5*(arma::trace(sol.Pa*sol.Ka)+arma::trace(sol.Pb*sol.Kb));
    if(sol.Pa.n_rows == sol.Pa_im.n_rows && sol.Pa.n_cols == sol.Pa_im.n_cols)
      sol.en.Exc+=0.5*(arma::trace(sol.Pa_im*sol.Ka_im)+arma::trace(sol.Pb_im*sol.Kb_im));
  }
//# 811 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"

//# 822 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
  
  if(freeze.size()>0) {

    freeze_orbs(freeze,sol.Ca,S,sol.Ha,verbose);
    freeze_orbs(freeze,sol.Cb,S,sol.Hb,verbose);
  }
  fflush(stdout);
  if(dimcalc) {
    
    std::vector<nucleus_t> nuclei(basisp->get_nuclei());
    for(size_t i=0;i<nuclei.size();i++)
      if(nuclei[i].r.x!=0.0 || nuclei[i].r.y!=0.0)
        throw std::logic_error("Nuclei must be on z axis for dimer calculation!\n");
    
    arma::ivec mvals(basisp->get_m_values());
    

    sol.Ha=block_m(sol.Ha,mvals);
    sol.Hb=block_m(sol.Hb,mvals);
  }
  
  sol.en.Ekin=arma::trace(sol.P*T);
  sol.en.Enuca=arma::trace(sol.P*Vnuc);
  sol.en.Enucr=Enuc;
  sol.en.Eone=arma::trace(sol.P*Hcore);
  sol.en.Ecoul=0.5*arma::trace(sol.P*sol.J);
  
  sol.en.Eel=sol.en.Ecoul+sol.en.Exc+sol.en.Eone+sol.en.Enl;
  sol.en.E=sol.en.Eel+sol.en.Enucr;
  
//# 876 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
  if(!arma::is_finite(sol.Ha)) {
//# 887 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
    throw std::runtime_error("Alpha Fock operator is not finite.\n");
  }
  if(!arma::is_finite(sol.Hb)) {
//# 900 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
    throw std::runtime_error("Beta Fock operator is not finite.\n");
  }

  if(!std::isfinite(sol.en.E)) {
//# 935 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
    std::ostringstream oss;
    oss << "\nSomething wrong with total energy " << sol.en.E <<"?\nEnding program.\n";
    throw std::runtime_error(oss.str());
  }

  std::cout << "\nEXITING SCF::Fock_UDFT\n";
}


