%Compute the differential of the gradient of the cost w.r.t. measurements
%function DygradPhi=bearingCostGeneral_DgradientMeasurements(YEval,YGoal,funs)
function DygradPhi=bearingCostGeneral_DgradientMeasurements(YEval,YGoal,funs)
df=funs.df;
ddf=funs.ddf;
[dimY,NY]=size(YEval);
DygradPhi=zeros(dimY,dimY,NY);

c=min(1,max(-1,sum(YGoal.*YEval)));

for iX=1:NY
    Yei=YEval(:,iX);
    Ygi=YGoal(:,iX);
    ci=c(iX);
    dfci=df(ci);
    ddfci=ddf(ci);
    DygradPhi(:,:,iX)=bearingCostGeneral_termsDygrad(Yei,Ygi,dfci,ddfci);
end
