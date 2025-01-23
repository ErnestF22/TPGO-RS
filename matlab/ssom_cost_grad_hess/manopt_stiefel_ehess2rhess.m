function rhess = manopt_stiefel_ehess2rhess(X, egrad, ehess, H)
%Useful for EHESS to RHESS CONVERSION (egrad also needed)
%Thanks to Manopt's stiefelfactory
XtG = multiprod(multitransp(X), egrad);
symXtG = multisym(XtG);
HsymXtG = multiprod(H, symXtG);
rhess = stiefel_tangentProj(X, ehess - HsymXtG);
end