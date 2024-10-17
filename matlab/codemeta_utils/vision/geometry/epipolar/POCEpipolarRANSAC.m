function POCEpipolarRANSAC
resetRands()
load POCEpipolarRANSAC_dataset
ransacThreshold=1e-3;
NOutliers=30;

imgSize=size(Im1);
K(3,3)=1;
halfImgSize=imgSize/2;
for k=1:2
    K(k,k)=halfImgSize(k);
    K(k,3)=halfImgSize(k);
end

x1=[x1 flipud(sfm_rawTransformFeature(K,2*rand(2,NOutliers)-1))];
x2=[x2 flipud(sfm_rawTransformFeature(K,2*rand(2,NOutliers)-1))];

x1Norm=sfm_rawTransformFeature(K,x1,'invert');
x2Norm=sfm_rawTransformFeature(K,x2,'invert');


figure(1)
sfm_rawDisplayMatch(Im1,Im2,x1,x2);

[E,residuals,flagInlier,allResiduals,ESamples] = sfm_rawEssentialRansac8pt(x1Norm,x2Norm,300,ransacThreshold);

[e,flagInlier]=sfm_rawMatchFilterWithEssential(x1Norm,x2Norm,E,'threshold',ransacThreshold);
x1Inlier=x1(:,flagInlier);
x2Inlier=x2(:,flagInlier);
figure(2)
sfm_rawDisplayMatch(Im1,Im2,x1Inlier,x2Inlier);

[e,flagInlierAuto]=sfm_rawMatchFilterWithEssential(x1Norm,x2Norm,E);
x1InlierAuto=x1(:,flagInlierAuto);
x2InlierAuto=x2(:,flagInlierAuto);
figure(3)
sfm_rawDisplayMatch(Im1,Im2,x1InlierAuto,x2InlierAuto);

Q=essential_fromE(ESamples);
[QMax,d,f,idxMax]=essential_kernelmode(Q,0.01,'methodDistance','sparse');

EMax=essential_getE(QMax);
[e,flagInlierMax]=sfm_rawMatchFilterWithEssential(x1Norm,x2Norm,EMax);
x1InlierMax=x1(:,flagInlierMax);
x2InlierMax=x2(:,flagInlierMax);
figure(4)
sfm_rawDisplayMatch(Im1,Im2,x1InlierMax,x2InlierMax);

save([mfilename '_data'])
