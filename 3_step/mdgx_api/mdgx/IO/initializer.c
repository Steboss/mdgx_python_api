#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "../mdgx.h"

//-----------------------------------------------------------------------------
// GetEwaldNamelist: this function reads in Ewald parameters, as would be
//                   found in the &ewald namelist in PMEMD or SANDER command
//                   input files.
//
// Arguments:
//   dcinp:   direct space command information
//   rcinp:   reciprocal space command information
//   inp:     the input file (mdin)  --> for the moment we are not interested in using it
//-----------------------------------------------------------------------------

void dircon_reccon(dircon *dcinp, reccon *rcinp)
{ 
  int i, collect, genordr;

  // Default settings for SPME
  //initialize the reccon
  rcinp->ng = (int*)malloc(3*sizeof(int));
  rcinp->ng[0] = -1;
  rcinp->ng[1] = -1;
  rcinp->ng[2] = -1;
  rcinp->ordr[0] = 4;
  rcinp->ordr[1] = 4;
  rcinp->ordr[2] = 4;
  rcinp->S = -1.0;
  // Default settings for MLE
  rcinp->ggordr = 8;
  rcinp->nlev = 1;
  rcinp->nslab = 1;
  rcinp->nstrip = 1;
  for (i = 0; i < 4; i++) {
    rcinp->cfac[i] = -1.0;
    rcinp->PadYZ[i] = -1;
  }
  rcinp->cfac[0] = 1.0;
  rcinp->PadYZ[0] = 1;
  //initialize dircon

  double Ecut;       // The electrostatic non-bonded cutoff
  double Vcut;       // The van-der Waals non-bonded cutoff
  double Mcut;       // The maximum of Ecut and Vcut
  double MaxDens;    // The maximum expected density of atoms (determines the
                     //   storage size in direct space decomposition cells)
  double invMcut;    // The inverse of Mcut
  double invEcut;    // The inverse of Ecut
  double sigma;      // The Gaussian width for spreading charges in preparation
                     //   for Ewald calculations

  ///////
  dcinp->Dtol = 6.0e-6;
  dcinp->lkpspc = 0.0625;
  dcinp->ewcoeff = -1.0;
  dcinp->LRvdw = 1;
  // Set the electrostatic cutoff to zero so that direct space interactions
  // are not counted.  This will generate other issues for exclusions, but
  // those will be dealt with in the routines for making lookup tables and
  // adjusting bonded interactions.
  dcinp->Ecut = 0.0;


  // Compute the Ewald coefficient and then the Gaussian spread
  if (rcinp->S <= 0.0 && dcinp->ewcoeff <= 0.0) {
    dcinp->ewcoeff = EwaldCoefficient(dcinp->Ecut, dcinp->Dtol);
    rcinp->S = 0.5/dcinp->ewcoeff;
  }
  else {
    if (rcinp->S > 0.0) {
      dcinp->ewcoeff = 0.5/rcinp->S;
    }
    else {
      rcinp->S = 0.5/dcinp->ewcoeff;
    }
    dcinp->Dtol = (1.0 - erf(dcinp->ewcoeff*dcinp->Ecut))/dcinp->Ecut;
  }

}

/*
//-----------------------------------------------------------------------------
// main
//-----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  int i, isys; //miscellanea
  reccon rcinp; //reciprocal space control
  dircon dcinp; //direct space control
  //coord* crd;  //these can be passed with parmed.py 
  cellgrid* CG;
  bckit* PPk;
  //prmtop* tp; //this is the topology, parmed.py
  trajcon tj;  //trajectory control system
  ipqcon ipqinp; //input for ipolq charg
  configs cfsinp; 
  spdata spelist;
  FrcTab Etab, EHtab;
  Energy* sysUV;
  execon etimers;
  cdftrj* Acdf;
  fset myfit;
  prmset myparms;

  //initialize dircon and reccon 
  dircon_reccon(dcinp,rcinp);
}
*/