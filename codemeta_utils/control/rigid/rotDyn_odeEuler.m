%Integrate second order ODE on the space of rotations using Euler's method
%function [t,x]=rotDyn_odeEuler(dx,tspan,x0,options)
function [t,x]=rotDyn_odeEuler(dx,tspan,x0,options)
[t,x]=generic_odeEuler(dx,tspan,x0,@funTransition,options);

function x=funTransition(x,dt,dx)
dw=rotDyn_inputUnpack(dx);
[R,w]=rotDyn_stateUnpack(x);
R=rot_exp(R,rot_hat(R,dt*w));
w=w+dt*dw;
x=rotDyn_statePack(R,w);
