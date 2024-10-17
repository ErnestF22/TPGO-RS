function POCDynamicSfM
load('sampleDynamicSfMDataset_original')
wb=wIMU;
dwb=J\(Gammab-rotDyn_gyroscopicTerm(wb,'inertiaMatrix',J));
hwb=hat3(wb);
bpp=permute(alphaIMU,[1 3 2]);

NFrames=size(wb,2);
NX=size(xb,2);

A=repmat(P,[1 1 NFrames]);
%dispError(matStack(xb)-matStack(multiprod(A,MTruth))*STruth+matStack(multiprod(multiprod(A,TTruth),ones(1,NX))))
%dispError(matStack(xb)-[matStack(multiprod(A,MTruth)) -vec(multiprod(A,TTruth))]*homogeneous(STruth,4))
%dispError(xb-multiprod(multiprod(A,MTruth),STruth)+multiprod(multiprod(A,TTruth),ones(1,NX)))
%dispError(multivec(xb+multiprod(multiprod(A,TTruth),ones(1,NX)))-multiprodMatVec(multikron(STruth',A),multivec(MTruth)))
bM=multivec(xb+multiprod(multiprod(A,TTruth),ones(1,NX)));
AM=multikron(STruth',A);
%dispError(bM-multiprodMatVec(AM,multivec(MTruth)))


Ap=-multiprod(P,hwb);
Bp=-A;
%dispError(matStack(dxb)-matStack(multiprod(Ap,MTruth))*STruth+vec(multiprod(Ap,TTruth)-multiprod(Bp,TpTruth))*ones(1,NX))
%dispError(matStack(dxb)-[matStack(multiprod(Ap,MTruth)) -vec(multiprod(Ap,TTruth)-multiprod(Bp,TpTruth))]*homogeneous(STruth,4))
%dispError(dxb-multiprod(multiprod(Ap,MTruth),STruth)+multiprod(multiprod(Ap,TTruth)-multiprod(Bp,TpTruth),ones(1,NX)))
%dispError(multivec(dxb+multiprod(multiprod(Ap,TTruth)-multiprod(Bp,TpTruth),ones(1,NX)))-multiprodMatVec(multikron(STruth',Ap),multivec(MTruth)))
bM=[bM;multivec(dxb+multiprod(multiprod(Ap,TTruth)-multiprod(Bp,TpTruth),ones(1,NX)))];
AM=[AM;multikron(STruth',Ap)];
%dispError(bM-multiprodMatVec(AM,multivec(MTruth)))

hwbSq=multiprod(hwb,hwb);
hdwb=hat3(dwb);
App=multiprod(P,hwbSq-hdwb);
Bpp=-2*Ap;
Cpp=-A;
%dispError(matStack(ddxb)-matStack(multiprod(App,MTruth))*STruth+vec(multiprod(App,TTruth)-multiprod(Bpp,TpTruth)-multiprod(Cpp,bpp)+multiprod(Cpp,multiprod(MTruth,g)))*ones(1,NX))
%dispError(ddxb-multiprod(multiprod(App,MTruth),STruth)+multiprod(multiprod(App,TTruth)-multiprod(Bpp,TpTruth)-multiprod(Cpp,bpp)+multiprod(Cpp,multiprod(MTruth,g)),ones(1,NX)))
%dispError(ddxb-multiprod(multiprod(App,MTruth),STruth)+multiprod(multiprod(App,TTruth)-multiprod(Bpp,TpTruth)-multiprod(Cpp,bpp),ones(1,NX))+multiprod(multiprod(Cpp,multiprod(MTruth,g)),ones(1,NX)))
%dispError(multivec(ddxb+multiprod(multiprod(App,TTruth)-multiprod(Bpp,TpTruth)-multiprod(Cpp,bpp),ones(1,NX)))-multiprodMatVec(multikron(STruth',App),multivec(MTruth))+multivec(multiprod(Cpp,multiprod(MTruth,g*ones(1,NX)))))
G=g*ones(1,NX);
%dispError(multivec(ddxb+multiprod(multiprod(App,TTruth)-multiprod(Bpp,TpTruth)-multiprod(Cpp,bpp),ones(1,NX)))-multiprodMatVec(multikron(STruth',App)-multikron(G',Cpp),multivec(MTruth)))
bM=[bM;multivec(ddxb+multiprod(multiprod(App,TTruth)-multiprod(Bpp,TpTruth)-multiprod(Cpp,bpp),ones(1,NX)))];
AM=[AM;multikron(STruth',App)-multikron(G',Cpp)];
r=zeros(1,size(AM,3));
for ix=1:size(AM,3)
    r(ix)=rank(AM(:,:,ix));
end
if ~all(r==9)
    warning('Problem ill-posed for some time instants')
end
%dispError(bM-multiprodMatVec(AM,multivec(MTruth)))
MEst=reshape(multisolve(AM,bM),3,3,[]);
%dispError(MEst-MTruth)
%MEst=rot_proj([MEstPart;zeros(1,3,size(MEstPart,3))]);
W=[matStack(xb);matStack(dxb);matStack(ddxb)];
[UW,SW,VW]=svd(W,'econ');
MEstPre=UW(:,1:4);
SEstPre=SW(1:4,1:4,:)*VW(:,1:4)';
%subspace(SEstPre',STruth')
x=SEstPre'\ones(NX,1);
STransform=[eye(3) zeros(3,1);x'];
SEstSym=STransform*SEstPre;
MEstSym=MEstPre/STransform;
STransformCenter=[eye(3) -mean(SEstSym(1:3,:),2);zeros(1,3) 1];
SEstSymCenter=STransformCenter*SEstSym;
MEstSymCenter=MEstSym/STransformCenter;
disp(subspace(SEstSymCenter(1:3,:)',STruth'))
%dispError(MEstSymCenter*SEstSymCenter-W)

NT=size(xb,3);
AM=[num2blkdiag(A); num2blkdiag(Ap); num2blkdiag(App)];
MEstSymReduced=AM\MEstSymCenter(:,1:3);
SEstSymReduced=SEstSymCenter(1:3,:);
%dispError(AM*MEstSymReduced-MEstSymCenter(:,1:3))
[MEstSymReducedAdjusted,STransformMRotations]=sfm_rawRotationsAdjust(MEstSymReduced,'method','linear');
SEstSymReducedAdjusted=STransformMRotations\SEstSymReduced;
%disp(subspace(MEstSymReducedAdjusted,matStack(MTruth)))
STransformMTruth=MEstSymReducedAdjusted\matStack(MTruth);
STransformM=MEstSymReducedAdjusted(1:3,1:3)'*Rbc';
%dispError(matStack(MTruth)-MEstSymReducedAdjusted*STransformMTruth)
dispError(STruth-STransformM\SEstSymReducedAdjusted)

%% Translation component using a single linear system
TAllTruth=[-vec(multiprod(A,TTruth));
    -vec(multiprod(Ap,TTruth)-multiprod(Bp,TpTruth));
    -vec(multiprod(App,TTruth)-multiprod(Bpp,TpTruth)-multiprod(Cpp,bpp)+multiprod(Cpp,multiprod(MTruth,g)))
    ];
%dispError(MEstSymCenter(:,4)-TAllTruth)
AMat=[num2blkdiag(-A) sparse(2*NT,3*NT)];
ApMat=[num2blkdiag(-Ap) num2blkdiag(Bp)];
AppMat=[num2blkdiag(-App) num2blkdiag(Bpp)];
xTTruth=[vec(TTruth);vec(TpTruth)];
% TAllB=[AMat*xTTruth;
%     ApMat*xTTruth;
%     AppMat*xTTruth-vec(-multiprod(Cpp,bpp))+norm(g)*num2blkdiag(Cpp)*MEstSymReducedAdjusted*STransformMTruth(:,3)
%     ];
ATVec=[AMat;ApMat;AppMat];
bTVec=-[sparse(4*NT,1);vec(-multiprod(Cpp,bpp))-norm(g)*num2blkdiag(Cpp)*MEstSymReducedAdjusted*STransformMTruth(:,3)];
TAllB=ATVec*xTTruth+bTVec;
dispError(TAllB-TAllTruth)
WRed=W-AM*MEstSymReduced*SEstSymReduced;
disp(norm(WRed-bTVec*ones(1,24)-ATVec*xTTruth*ones(1,24),'fro'));
disp(norm(xTTruth-ATVec\(mean(WRed,2)-bTVec)))

%% Translation component using multiple, frame-wise systems
% Recall that TTruth contains taub, while TpTruth contains nub
xTTruthC=[squeeze(TTruth);squeeze(TpTruth)];
% TAllC=[-multiprodMatVec([A zeros(size(A))],xTTruthC);
%     -multiprodMatVec([Ap -Bp],xTTruthC);
%     -multiprodMatVec([App -Bpp],xTTruthC)+multiprodMatVec(Cpp,squeeze(bpp))-multiprodMatVec(Cpp,multiprodMatVec(MTruth,g))
%     ];
% TAllC=-multiprodMatVec([A zeros(size(A));Ap -Bp;App -Bpp],xTTruthC)...
%     +[zeros(4,size(A,3));multiprodMatVec(Cpp,squeeze(bpp))-multiprodMatVec(Cpp,multiprodMatVec(MTruth,g))];
% TAllC=-multiprodMatVec([A zeros(size(A));Ap -Bp;App -Bpp],xTTruthC)...
%     +[zeros(4,size(A,3));multiprodMatVec(Cpp,squeeze(bpp))+norm(g)*multiprodMatVec(Cpp,reshape(MEstSymReducedAdjusted*STransformMTruth(:,3),3,[]))];
AT=-[A zeros(size(A));Ap -Bp;App -Bpp];
bT=@(s) -[zeros(4,size(A,3));multiprodMatVec(Cpp,squeeze(bpp))+norm(g)*multiprodMatVec(Cpp,reshape(MEstSymReducedAdjusted*s,3,[]))];
N=size(MEstSymCenter,1);
ord=matStack(reshape(1:N,2,N/6,[]));
%TAllC=multiprodMatVec(AT,xTTruthC)-bT;
%dispError(TAllC-reshape(MEstSymCenter(ord,4),size(ord)))
%dispError(reshape(MEstSymCenter(ord,4),size(ord))+bT(STransformMTruth(:,3))-multiprodMatVec(AT,xTTruthC))
xT=@(s) multisolve(AT,bT(s)+reshape(MEstSymCenter(ord,4),size(ord)));
%dispError(xT(STransformM(:,3))-xTTruthC)
%s=sphere_randn();
%dispError(multiprodMatVec(AT,xT(s))-bT(s)-reshape(MEstSymCenter(ord,4),size(ord)))
s=STransformM(:,3);
xTs=xT(s);

%% First derivative and smoothness constraints on translations
filterOrder=4;
filterWindow=21;
filterHalfWindow=floor(filterWindow/2);
TSampling=mean(diff(t));
%matrix representation of numerical derivative based on S-G filter
[~,g] = sgolay(filterOrder,filterWindow);   
MDer=-spconvmtx(sparse(g(:,2)),NFrames)/TSampling;
MDer=MDer(1+filterHalfWindow:end-filterHalfWindow,:);
idxValid=1+filterHalfWindow:NFrames-filterHalfWindow;
% TsbEst=multiprodMatVec(invR(Rbs),reshape(xTTruth(1:NFrames*3),3,[]));%xTs(1:3,:));
% dTsbEst=multiprodMatVec(invR(Rbs),reshape(xTTruth((NFrames*3+1):end),3,[]));%xTs(4:6,:));
TsbEst=reshape(num2blkdiag(invR(Rbs))*xTTruth(1:NFrames*3),3,[]);
dTsbEstNumVec=kron(MDer,speye(3))*num2blkdiag(invR(Rbs))*xTTruth(1:NFrames*3);
dTsbEst=reshape(num2blkdiag(invR(Rbs))*xTTruth((NFrames*3+1):end),3,[]);
%funCheckDerInterpInterp(t,TsbEst,dTsbEst)
%dTsbEstNum=TsbEst*MDer';
%ADerxT=reshape(dTsbEstNumVec-num2blkdiag(invR(Rbs))*xTTruth((NFrames*3+1):end),3,[]);
ADerxTVec=[kron(MDer,speye(3))*num2blkdiag(invR(Rbs)) -num2blkdiag(invR(Rbs))]*xTTruth;
ADerxT=reshape(ADerxTVec,3,[]);
%disp(rmse(ADerxT(:,idxValid)))
%plot(t(idxValid),ADerxT(:,idxValid)')
idxValidMat=subMatrix(reshape(1:3*NFrames,3,[]),1:3,idxValid);
ADerxTVecValid=ADerxTVec(idxValidMat,:);
disp(rmse(reshape(ADerxTVecValid,3,[])))

%plot(t(idxValid),dTsbEstNum(:,idxValid)',t,dTsbEst,'x')
%dispError(dTsbEstNum-dTsbEstNumVec)
%ADer=[kron(MDer,speye(3))*num2blkdiag(invR(Rbs)) -num2blkdiag(invR(Rbs))];
%ADerxT=reshape(ADer*xTTruth,3,[]);
%dispError(ADerxT(:,idxValid))

%% Second derivative and smoothness constraints on translations
%funCheckDerInterpInterp(t,dTsbEst,multiprodMatVec(invR(RbsTruth),alphaIMU-multiprodMatVec(RbsTruth,gravityVector)))
%rmse(multiprodMatVec(RbsTruth,gravityVector),...
%    -norm(gravityVector)*reshape(MEstSymReducedAdjusted*STransformM(:,3),3,[]))
%funCheckDerInterpInterp(t,dTsbEst,...
%    multiprodMatVec(invR(RbsTruth),alphaIMU+norm(gravityVector)*reshape(MEstSymReducedAdjusted*STransformM(:,3),3,[])))
bDDer=vec(multiprodMatVec(invR(RbsTruth),alphaIMU+norm(gravityVector)*reshape(MEstSymReducedAdjusted*STransformM(:,3),3,[])));
ADder=[sparse(3*NFrames,3*NFrames) kron(MDer,speye(3))*num2blkdiag(invR(Rbs))];
%rmse(dTsbEst,reshape(ADder*xTTruth,3,[]))
ADderxTVec=ADder*xTTruth-bDDer;
dispError(ADderxTVec(idxValidMat))

%% Derivative and smoothness constraints on rotations
Rwb=rot_exp(eye(3),rot_hat(eye(3),-TSampling*wb(:,1:end-1)));
ARDer=dynamicSfM_rotationSparseConvolutionMatrix(Rwb);
dispError(ARDer*MEstSymReducedAdjusted)


function D=num2blkdiag(A)
ACell=squeeze(num2cell(A,[1 2]));
D=spblkdiag(ACell{:});


function dispError(x)
disp(norm(vec(x),'inf'))

function xVec=multivec(x)
sz=size(x);
xVec=reshape(x,[],sz(end));

function x=multisolve(AM,bM)
Nx=size(AM,3);
d=size(AM,2);
x=zeros(d,Nx);
for ix=1:Nx
    x(:,ix)=AM(:,:,ix)\bM(:,ix);
end
