function rotDyn_controlPD_cost_test
global RReference
resetRands()
RReference=rot_randn();
[R,~,~,~,w,~,dw]=rot_geodFun([],[],'speed','cubic');

funCheckDer(@(t) funAndDer(R(t),w(t),dw(t)))


function [phi,dphi]=funAndDer(R,w,dw)
global RReference
J=diag([5;2;1]);
method='packed';
switch lower(method)
    case 'localcoordinates'
        [phi,gradphi]=rotDyn_controlPD_cost(R,w,'inertiaMatrix',J,'RReference',RReference,'methodGradient','localCoordinates');
        dphi=gradphi'*[w;dw];
    case 'packed'
        [phi,gradphi]=rotDyn_controlPD_cost(R,w,'inertiaMatrix',J,'RReference',RReference,'methodGradient','packed');
        dphi=rotDyn_metricPacked(gradphi,rotDyn_statePack(R*hat3(w),dw));
end

