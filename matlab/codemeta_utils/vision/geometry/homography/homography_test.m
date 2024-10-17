function homography_test
flagDisplay=true;
flagCheckPlane=true;
flagCheckHomographyMap=true;
flagCheckTriangulation=true;

resetRands(1)
N=sphere_randn([1;0;0],0.1);
d=0;
NVec=planeNdToNVec(N,d);
NPoints=30;
X=planeGeneratePoints(NVec,NPoints);

methodAbsolutePoses='reference';
R1=rot([0;pi/2;0])*rot_randn(eye(3),0.1);
T1=-R1*[0;0;5];
NVec1=rigidTransform(R1,T1,NVec,'planes','methodAbsolutePoses',methodAbsolutePoses);
[N1,d1]=planeNVecToNd(NVec1);
R2=rot([0;pi/2;0])*rot_randn(eye(3),1);
T2=-R2*[0;0;5];

%normalize scene so that d1=1
T1=T1/abs(d1);
T2=T2/abs(d1);
X=X/abs(d1);
N1=-N1;
d1=1;
NVec1=planeNdToNVec(N1,d1);


G1=RT2G(R1,T1);
G2=RT2G(R2,T2);
X1=rigidTransformG(G1,X,'methodAbsolutePoses',methodAbsolutePoses);
X2=rigidTransformG(G2,X,'methodAbsolutePoses',methodAbsolutePoses);

if flagCheckPlane
    disp('N1''*X1-d1')
    disp(N1'*X1-d1)
end

G=cat(3,G1,G2);
G21=computeRelativePoseFromG(G2,G1,'methodAbsolutePoses',methodAbsolutePoses);
%[R21,T21]=G2RT(G21);

x=projectFromG(G,X,'methodAbsolutePoses',methodAbsolutePoses);
x1=x(:,:,1);
x2=x(:,:,2);


if flagDisplay
    figure(1)
    draw3dcameraFromG(G1,'methodAbsolutePoses',methodAbsolutePoses,'scale',0.2)
    hold on
    draw3dcameraFromG(G2,'methodAbsolutePoses',methodAbsolutePoses,'scale',0.2)
    plotPoints(X)
    hold off
    axis equal

    figure(2)
    plotPoints(x(:,:,1),'r');
    hold on
    plotPoints(x(:,:,2),'b');
    hold off
end

H=homographyFromG(G21,NVec1);
HEst=homographyNormalize(homographyEstimateLinear(x1,x2),x1,x2);

if flagCheckHomographyMap
    x2Est=homographyMap(H,x(:,:,1));
    disp('max(abs(x(:,:,2)-x2Est),[],2)')
    disp(max(abs(x(:,:,2)-x2Est),[],2))
end

if flagCheckTriangulation
    lambdaRatio=homographyTriangulateDepths(H,x1,x2);
    disp('max(abs(X2(3,:)./X1(3,:)-lambdaRatio(2,:)),[],2)')
    disp(max(abs(X2(3,:)./X1(3,:)-lambdaRatio(2,:)),[],2))
end



x1a=homogeneous(x1(:,1),3);
x2a=homogeneous(x2(:,1),3);
x2aOrth=computeOrthogonalComplement(x2a);
%disp(x2aOrth'*H*x1a)
%disp(kron(x1a',x2aOrth')*H(:))
%disp(kron(x1a,x2aOrth)'*H(:))
disp('[H HEst]')
disp([H HEst])
[R21,T21]=G2RT(G21);
[R21Est,T21Est,N1Est,lambda,tp]=homographyToRT(HEst,x1,x2);
disp('[R21 R21Est]')
disp([R21 R21Est])
disp('[T21/d1 T21Est]')
disp([T21/d1 T21Est])
disp('[N1 N1Est]')
disp([N1 N1Est])
G21Est=RT2G(R21Est, T21Est);

G21EstB=poseEstimationHomographyG(x1,x2);

disp('cnormalize([G21(1:3,:) G21Est(1:3,:) G21EstB(1:3,:)])');
disp(cnormalize([G21(1:3,:) G21Est(1:3,:) G21EstB(1:3,:)]));

save([mfilename '_data'])
