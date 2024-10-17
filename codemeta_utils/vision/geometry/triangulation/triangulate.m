function XEst=triangulate(x,P)
XEstLin=triangulate_lin(x,P);
XEst=triangulate_nonlin(x,P,XEstLin);
