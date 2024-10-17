function POCAffineSfm
load('sampleDynamicSfMDataset_interpolation')

F=size(xb,3);
P=size(xb,2);
W=matStack(xb);

c=mean(W,2);
WCentered=W-c*ones(1,P);
%rank(WCentered)

[UWcentered,SWCentered,VWCentered]=svd(WCentered,'econ');
MAffine=UWcentered(:,1:3);
SAffine=SWCentered(1:3,1:3)*VWCentered(:,1:3)';
%[MMetric,KMetric]=sfm_rawRotationsAdjust(MAffine,'method','linear');

[RReducedEst,KEst]=affineSfMMetricUpgrade(matUnstack(MAffine,2));
SEst=KEst\SAffine;

REst=zeros(3,3,F);
for f=1:F
    REst(:,:,f)=orthCompleteBasis(RReducedEst(:,:,f)')';
end

flagFlip=true;
if flagFlip
    KFlip=diag([1 1 -1]);
    REst=multiprod(KFlip,multiprod(REst,KFlip));
    SEst=KFlip'*SEst;
end
RFirstFrame=RbsTruth(:,:,1);
KAlign=rot_proj(REst(:,:,1))'*RFirstFrame;
RAlign=multiprod(REst,KAlign);
SAlign=KAlign\SEst;

TbEst=[reshape(c,2,F); ones(1,F)];
TAlign=-multiprodMatVec(RAlign,TbEst);
plotPoints(SAlign,'x');
hold on
plotPoints(TAlign,'-')
hold off
keyboard