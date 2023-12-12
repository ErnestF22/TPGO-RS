%Compute gradient of bearing+range cost
function gradPhi=bearingCostGeneralRangesGradient(y,yg,ny,nyg,funs)
c=bearingComputeCosine(y,yg);
e=ny.*c-nyg;
d=size(y,1);
gradPhi=-sum((ones(d,1)*funs.df(e)).*yg,2);
