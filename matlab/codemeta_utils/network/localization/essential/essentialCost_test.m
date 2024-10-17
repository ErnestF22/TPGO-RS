function essentialCost_test
resetRands()
[gi,~,~,~,vgiVec]=rot3r3_randGeodFun();
[gj,~,~,~,vgjVec]=rot3r3_randGeodFun();
Qij=essential_fromG(gi(0),gj(0));

f=@(t) costDer(gi(t),gj(t),Qij,[vgiVec;vgjVec]);
df=@(t) derDder(gi(t),gj(t),Qij,[vgiVec;vgjVec]);

%check_der(f,'function',linspace(-5,5,200))
check_der(df,'function',linspace(-5,5,200))

function [c,dc]=costDer(gi,gj,Qij,vij)
[c,gradc]=essentialCost(G2R(gi),G2T(gi),G2R(gj),G2T(gj),Qij);
dc=vij'*gradc;

function [dc,ddc]=derDder(gi,gj,Qij,vij)
[~,gradc,Hessc]=essentialCost(G2R(gi),G2T(gi),G2R(gj),G2T(gj),Qij);
dc=vij'*gradc;
ddc=vij'*Hessc*vij;
