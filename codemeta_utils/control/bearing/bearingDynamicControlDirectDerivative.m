%Compute derivative of bearingDynamicControlDirect using bearing derivatives
function du=bearingDynamicControlDirectDerivative(y,dy,ddy,yg,funs,alpha)
du1=bearingControlDirectDerivative(y,dy,yg,funs);
du2=sum(ddy,2);
du=alpha*du1+du2;
