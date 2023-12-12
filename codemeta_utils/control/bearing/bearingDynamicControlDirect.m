%Compute the acceleration control from the bearing potential and derivatives
function u=bearingDynamicControlDirect(y,dy,yg,funs,alpha)
gPhi=bearingCostGeneral_gradient(y,yg,funs);
u=-alpha*gPhi+sum(dy,2);
