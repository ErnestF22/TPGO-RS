function rotLocFrobCost_test
t_node=testNetworkBuildTestNetwork('N',7);
[R,v]=rot_randGeodFun(t_node.Ritruth);
RRel=t_node.Rijtruth;
E=t_node.E;

G=rotLocFrobCostMatrix(E,RRel);

funCheckDer(@(t) costAndDer(E,R(t),RRel,v(t)))
funCompare(@(t) rotLocFrobCost(E,R(t),RRel), @(t) costFromMatrix(R(t),G))

function [c,dc]=costAndDer(E,R,RRel,v)
flagSeparateGradient=true;
if flagSeparateGradient
    c=rotLocFrobCost(E,R,RRel);
    gradc=rotLocFrobCost_grad(E,R,RRel);
else
    [c,gradc]=rotLocFrobCost(E,R,RRel);
end
    
dc=sum(multitrace(multiprod(multitransp(v),gradc)));

function c=costFromMatrix(R,G)
RVec=matStack(multitransp(R));
%c=trace(RVec'*G*RVec)/2;
c=trace(RVec(:)'*kron(eye(3),G)*RVec(:))/2;
