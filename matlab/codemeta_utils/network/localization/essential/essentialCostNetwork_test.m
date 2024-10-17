function essentialCostNetwork_test

t_node=testNetworkBuildTestNetwork();%buildTestTagNetwork();
[Ri0,Ti0]=testNetworkGetRotTransl(t_node);
%[~,Qij]=testNetworkGetRelativeEssential(t_node);
E=testNetworkGetEdges(t_node);
Qij=t_node.QEijtruth;

% [Rit,~,~,~,vi]=rot_randGeodFun(Ri0);
% [Tit,~,~,dTi]=real_randGeodFun(Ti0);
GiInit=RT2G(Ri0,Ti0);
[Git,~,~,~,vi]=rot3r3_randGeodFun(GiInit);
f=@(t) costDer(Git(t),Qij,E,vi);

%check_der(f,'function',linspace(-5,5,30))
plotfun(f,linspace(-1e-3,1e-3,30))

function [c,dc]=costDer(Gi,Qij,E,vi)
[c,gradc]=essentialCostNetwork(G2R(Gi),G2T(Gi),Qij,E);
dc=vi(:)'*gradc(:);
