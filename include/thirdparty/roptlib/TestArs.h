/*
This is the test file for the Brocokett problem defined in StieBrockett.h and StieBrockett.cpp.

---- WH
*/

// from   ROPTLIB-beta/test/TestProdStieSumBrockett.h

#ifndef TESTPRODSTIESUMBROCKETT_H
#define TESTPRODSTIESUMBROCKETT_H

/*Output to console*/
#include <iostream>

/*Generate random number*/
#include "others_randgen.h"

/*Computational time*/
#include <ctime>

/*If this test file is called from Matlab, then functions in DriverMexProb.h are used.*/
#include "test_DriverMexProb.h"

/*Problem related classes*/
#include "problems_Problem.h"
#include "problems_StieBrockett.h"
#include "problems_SphereTxRQ.h"
#include "problems_ProdStieSumBrockett.h"

/*Manifold related classes*/
#include "manifolds_Manifold.h"
#include "manifolds_MultiManifolds.h"
//#include "manifolds_Stiefel/StieVariable.h"
#include "manifolds_Stiefel.h"
#include "manifolds_SphereTx.h"

/*Linesearch based solvers*/
#include "solvers_RSD.h"
#include "solvers_RNewton.h"
#include "solvers_RCG.h"
#include "solvers_RBroydenFamily.h"
#include "solvers_RWRBFGS.h"
#include "solvers_RBFGS.h"
#include "solvers_LRBFGS.h"
#include "solvers_RGS.h"
#include "solvers_LRBFGSSub.h"
//#include "solvers_RBFGSLPSub.h"

/*Trust-region based solvers*/
#include "solvers_SolversSMTR.h"
#include "solvers_RTRSD.h"
#include "solvers_RTRNewton.h"
#include "solvers_RTRSR1.h"
#include "solvers_LRTRSR1.h"

/*The global head file*/
#include "others_def.h"

using namespace ROPTLIB;

/*The main test function*/
void testProdStieSumBrockett(void);

#endif // end of TESTPRODSTIESUMBROCKETT_H
