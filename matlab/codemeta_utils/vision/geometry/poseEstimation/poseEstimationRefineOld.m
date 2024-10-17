function [REst,TEst]=poseEstimationRefine(R,T,X,x)
G0=RT2G(R,T);
f=@(G) cost(G,X,x);
df=@(G) derCost(G,X,x);
[GEst,errorsGiEst]=lie_minimize(rot3r3_funs(),f,df,G0,'maxIt',2000,'stepsize',2); %,'showCost'
REst=G2R(GEst);
TEst=G2T(GEst);


function c=cost(G,X,x)
c=sum(poseEstimationCost_GFormat(G,X,x));

function d=derCost(G,X,x)
[~,gradc]=poseEstimationCost_GFormat(G,X,x);
d=sum(gradc,3);
