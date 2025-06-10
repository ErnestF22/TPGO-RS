function rhess = manopt_stiefel_ehess2rhess(X, egrad, ehess, Xdot)
%Useful for EHESS to RHESS CONVERSION (egrad also needed)
%Thanks to Manopt's stiefelfactory
XtG = multiprod(multitransp(X), egrad);
symXtG = multisym(XtG);
HsymXtG = multiprod(Xdot, symXtG);
rhess = stiefel_tangentProj(X, ehess - HsymXtG);
end