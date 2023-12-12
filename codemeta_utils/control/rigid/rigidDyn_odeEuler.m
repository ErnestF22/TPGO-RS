%Integrate second order ODE on the space of rigid poses using Euler's method
%function [t,x]=rigidDyn_odeEuler(dx,tspan,x0,options)
function [t,x]=rigidDyn_odeEuler(dx,tspan,x0,options)
[t,x]=generic_odeEuler(dx,tspan,x0,@funTransition,options);

function x=funTransition(x,dt,dx)
[dw,dv]=rigidDyn_inputUnpack(dx);
[R,w,T,v]=rigidDyn_stateUnpackRT(x);
R=rot_exp(R,rot_hat(R,dt*w));
w=w+dt*dw;
T=T+dt*v;
v=v+dt*dv;
x=rigidDyn_statePackRT(R,w,T,v);
