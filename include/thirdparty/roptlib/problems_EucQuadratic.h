/*
This file defines the class for the problem min_{x \in R^d} 0.5 * x^T A x, where A is a d by d symmetric positive definite matrix

Problem --> EucQuadratic

---- WH
*/

#ifndef EUCQUADRATIC_H
#define EUCQUADRATIC_H

#include "manifolds_Euclidean.h"
#include "problems_Problem.h"
#include "others_def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class EucQuadratic : public Problem{
	public:
		EucQuadratic(Vector M);
		virtual ~EucQuadratic(void);
		virtual realdp f(const Variable &x) const;
		virtual Vector &EucGrad(const Variable &x, Vector *result) const;
		virtual Vector &EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const;

        Vector A;
	};
}; /*end of ROPTLIB namespace*/
#endif /* end of EUCQUADRATIC_H */
