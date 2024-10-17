function [XEst,RbsEst,taubEst,nubEst,output]=dynamicSfm(xb,dxb,ddxb,Gammab,wIMU,alphaIMU,varargin)
J=eye(3);
P=[eye(2) zeros(2,1)];
g=gravityVector();
flagProjectMFactorRotations=false;
flagRetriangulateStructure=false;
flagCollect=nargout>4;
flagShowDiagnosticInfo=false;
RFirstFrame=eye(3);
methodTranslationFactor='original';
methodTranslationExtraction='multisystem';
flagDerivativeConstraint=false;
flagSecondDerivativeConstraint=false;
derivativeConstraintsWeight=1;
secondDerivativeConstraintsWeight=2;
flagDebugSubstituteRotations=false;
flagRotationDerivativeConstraint=false;
flagEstimateAlignedGravity=true;
rotationDerivativeConstraintWeight=1;
rotationExtractionWeights=[1 1 1];
%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'inertiamatrix'
            ivarargin=ivarargin+1;
            J=varargin{ivarargin};
        case 'rfirstframe'
            ivarargin=ivarargin+1;
            RFirstFrame=varargin{ivarargin};
        case 'showdiagnosticinfo'
            flagShowDiagnosticInfo=true;
        case 'projectmfactorrotations'
            flagProjectMFactorRotations=true;
        case 'retriangulatestructure'
            flagRetriangulateStructure=true;
        case 'rotationextractionweights'
            ivarargin=ivarargin+1;
            rotationExtractionWeights=varargin{ivarargin};
        case 't'
            ivarargin=ivarargin+1;
            t=varargin{ivarargin};
        case 'methodtranslationfactor'
            ivarargin=ivarargin+1;
            methodTranslationFactor=lower(varargin{ivarargin});
        case 'methodtranslationextraction'
            ivarargin=ivarargin+1;
            methodTranslationExtraction=lower(varargin{ivarargin});
        case 'rotationderivativeconstraint'
            flagRotationDerivativeConstraint=true;
            ivarargin=ivarargin+1;
            TSampling=lower(varargin{ivarargin});
        case 'rotationderivativeconstraintweight'
            ivarargin=ivarargin+1;
            rotationDerivativeConstraintWeight=lower(varargin{ivarargin});
        case 'translationderivativeconstraint'
            methodTranslationExtraction='singlesystem';
            flagDerivativeConstraint=true;
            [filterCoefficients,ivarargin]=parseFilterCoefficients(varargin,ivarargin);
        case 'flagtranslationderivativeconstraint'
            methodTranslationExtraction='singlesystem';
            ivarargin=ivarargin+1;
            flagDerivativeConstraint=lower(varargin{ivarargin});
        case 'translationderivativeconstraintweight'
            ivarargin=ivarargin+1;
            derivativeConstraintsWeight=lower(varargin{ivarargin});
        case 'translationsecondderivativeconstraint'
            methodTranslationExtraction='singlesystem';
            flagSecondDerivativeConstraint=true;
            [secondFilterCoefficients,ivarargin]=parseFilterCoefficients(varargin,ivarargin);
        case 'flagtranslationsecondderivativeconstraint'
            methodTranslationExtraction='singlesystem';
            ivarargin=ivarargin+1;
            flagSecondDerivativeConstraint=lower(varargin{ivarargin});
        case 'flagestimatealignedgravity'
            ivarargin=ivarargin+1;
            flagEstimateAlignedGravity=lower(varargin{ivarargin});
        case 'translationsecondderivativeconstraintweight'
            ivarargin=ivarargin+1;
            secondDerivativeConstraintsWeight=lower(varargin{ivarargin});
        case 'debugsubstituterotations'
            ivarargin=ivarargin+1;
            flagDebugSubstituteRotations=true;
            debugSubstituteRotations=lower(varargin{ivarargin});
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

% detect some conflicting configurations
if flagEstimateAlignedGravity && ~(flagDerivativeConstraint && flagSecondDerivativeConstraint)
    error('In order to estimate the gravity direction, both derivative constraints on the translations should be present')
end

if (flagDerivativeConstraint || flagSecondDerivativeConstraint) && ~strcmpi(methodTranslationExtraction,'singleSystem')
    error('Derivative constraints on the translations need ''methodTranslationExtraction'',''singleSystem''')
end

wb=wIMU;
dwb=J\(Gammab-rotDyn_gyroscopicTerm(wb,'inertiaMatrix',J));
hwb=hat3(wb);
hwbSq=multiprod(hwb,hwb);
hdwb=hat3(dwb);

NFrames=size(wb,2);
NX=size(xb,2);

%Prepare various coefficient matrices
A=repmat(P,[1 1 NFrames]);
Ap=-multiprod(P,hwb);
Bp=-A;

App=multiprod(P,hwbSq-hdwb);
Bpp=-2*Ap;
Cpp=-A;
bpp=permute(alphaIMU,[1 3 2]);

if flagCollect
    output.A=A;
    output.Ap=Ap;
    output.Bp=Bp;
    output.App=App;
    output.Bpp=Bpp;
    output.Cpp=Cpp;
    output.bpp=bpp;
end

%Matrix with all measurements
W=[matStack(xb);matStack(dxb);matStack(ddxb)];
if flagCollect
    output.W=W;
end
if flagShowDiagnosticInfo
    disp('Sum of svals over 4:');
    s=svd(W);
    disp(sum(s(5,:)));
end

%Factorization (reconstruction up to projectivity)
[UW,SW,VW]=svd(W,'econ');
MEstProjective=UW(:,1:4);
SEstProjective=SW(1:4,1:4,:)*VW(:,1:4)';
if flagCollect
    output.MEstProjective=MEstProjective;
    output.SEstProjective=SEstProjective;
end    
if flagShowDiagnosticInfo
    disp('Error of rank-4 projective factorization')
    disp(norm(W-MEstProjective*SEstProjective,'fro'))
end

%Recover normal to plane at infinity by imposing that the last row of the
%estimated structure should contain all ones
v=SEstProjective'\ones(NX,1);
KProjective=[eye(3) ones(3,1);v'];

%Apply transformation to obtain reconstruction up to an affine transformation
SEstAffine=KProjective*SEstProjective;
MEstAffine=MEstProjective/KProjective;
if flagCollect
    output.SEstAffine=SEstAffine;
    output.MEstSym=MEstAffine;
    output.KProjective=KProjective;
end    
if flagShowDiagnosticInfo
    disp('RMSE on fourth elements of SEstAffine being equal to one')
    disp(rmse(SEstAffine(4,:),1))
    disp('Error of rank-4 affine factorization')
    disp(norm(W-MEstAffine*SEstAffine,'fro'))
end

%Impose the constraint that the world reference frame is at the centroid of
%the points
KAffine=[eye(3) -mean(SEstAffine(1:3,:),2);zeros(1,3) 1];

%Apply transformation to obtain reconstruction up to a similarity
SEstSimilarity=KAffine*SEstAffine;
MEstSimilarity=MEstAffine/KAffine;
if flagCollect
    output.SEstSimilarity=SEstSimilarity;
    output.MEstSimilarity=MEstSimilarity;
    output.KAffine=KAffine;
end
if flagShowDiagnosticInfo
    disp('Error of rank-4 similarity factorization')
    disp(norm(W-MEstSimilarity*SEstSimilarity,'fro'))
end    

%Extract reduced motion matrix from the first factor
%with a reweight of the equations
AM=[num2blkdiag(A);
    num2blkdiag(Ap);
    num2blkdiag(App)];
MEstSimilarityReweighted=MEstSimilarity;
WReweighted=W;
idxSections=reshape(1:NFrames*3*2,NFrames*2,[]);
for iSection=1:3
    MEstSimilarityReweighted(idxSections(iSection),:)=rotationExtractionWeights(iSection)*MEstSimilarityReweighted(idxSections(iSection),:);
    AM(idxSections(iSection),:)=rotationExtractionWeights(iSection)*AM(idxSections(iSection),:);
    WReweighted(idxSections(iSection),:)=rotationExtractionWeights(iSection)*WReweighted(idxSections(iSection),:);
end
if flagRotationDerivativeConstraint
    Rwb=rot_exp(eye(3),rot_hat(eye(3),-TSampling*wb(:,1:end-1)));
    ARDer=dynamicSfM_rotationSparseConvolutionMatrix(Rwb);
    AM=[AM; rotationDerivativeConstraintWeight*ARDer];
    MEstSimilarityReweighted=[MEstSimilarityReweighted; zeros(3*(NFrames-1),4)];
end
MEstSimilarityReduced=AM\MEstSimilarityReweighted(:,1:3);
if flagShowDiagnosticInfo
    disp('Error on LS problem for extracting reduced motion matrix')
    disp(norm(AM*MEstSimilarityReduced-MEstSimilarityReweighted(:,1:3),'fro'))
end

%Structure in non-homogeneous coordinates
SEstSimilarityReduced=SEstSimilarity(1:3,:);
if flagCollect
    output.MEstSimilarityReduced=MEstSimilarityReduced;
    output.SEstSimilarityReduced=SEstSimilarityReduced;
    output.AM=AM;
end    
if flagShowDiagnosticInfo
    disp('Error on rank-4 factorization with reduced motion and structure matrices')
    %disp(norm(W-[AM(1:6*NFrames,:)*MEstSimilarityReduced MEstSimilarity(:,4)]*SEstSimilarity,'fro'))
    disp(norm(WReweighted-[AM(1:6*NFrames,:)*MEstSimilarityReduced MEstSimilarity(:,4)]*homogeneous(SEstSimilarityReduced,4),'fro'))
    disp('  Note: this will not be zero with ideal data if the rotationExtractionWeights are not all the same')
end

%Metric upgrade
[MEstMetric,KMetric]=sfm_rawRotationsAdjust(MEstSimilarityReduced,'method','linear');
SEstMetric=KMetric\SEstSimilarityReduced;

%correct sign
s=median(sign(multidet(matUnstack(MEstMetric))));
if s<0
    MEstMetric=-MEstMetric;
    SEstMetric=-SEstMetric;
    KMetric=-KMetric;
end

if flagProjectMFactorRotations
    if flagCollect
        output.MEstMetricPreProjection=MEstMetric;
    end
    MEstMetric=matStack(rot_proj(matUnstack(MEstMetric)));
end

if flagShowDiagnosticInfo
    disp('Error on rank-4 factorization with metric motion and structure matrices')
    disp(norm(WReweighted-AM(1:6*NFrames,:)*MEstMetric*SEstMetric-MEstSimilarityReweighted(1:6*NFrames,4)*ones(1,NX),'fro'))
    disp('  Note: this will not be zero with ideal data if the rotationExtractionWeights are not all the same')
end

if flagRetriangulateStructure
    if flagCollect
        output.SEstMetricPreRetriangulation=SEstMetric;
    end
    SEstMetric=(AM(1:6*NFrames,:)*MEstMetric)\(WReweighted-MEstSimilarityReweighted(1:6*NFrames,4)*ones(1,NX));
end

if flagCollect
    output.MEstMetric=MEstMetric;
    output.SEstMetric=SEstMetric;
    output.STransformMetric=KMetric;
end


%Align reconstruction with a rotation so that the first frame is
%(approximately) RFirstFrame
KAlign=rot_proj(MEstMetric(1:3,1:3))'*RFirstFrame;
MEstAligned=MEstMetric*KAlign;
SEstAligned=KAlign\SEstMetric;
if flagCollect
    output.MEstAligned=MEstAligned;
    output.SEstAligned=SEstAligned;
    output.KAlign=KAlign;
end

%Assign outputs for rotations and structure
RbsEst=rot_proj(matUnstack(MEstAligned));
XEst=SEstAligned;

if flagDebugSubstituteRotations
    RbsEst=debugSubstituteRotations;
    MEstAligned=matStack(RbsEst);
end

%Recover translations and linear velocities
switch lower(methodTranslationFactor)
    case 'original'
        TFactor=MEstSimilarityReweighted(1:6*NFrames,4);
    case 'mean'
        TFactor=mean(W,2);
    case 'reduced'
        WReduced=WReweighted-AM(1:6*NFrames,:)*MEstAligned*SEstAligned;
        TFactor=mean(WReduced,2);
    otherwise
        error('Method for translation factor not recognized')
end
if flagCollect
    output.methodTranslationFactor=methodTranslationFactor;
    output.TFactor=TFactor;
end

if flagShowDiagnosticInfo
    fprintf('Method for translation extraction: %s\n\n',methodTranslationExtraction);
end

switch lower(methodTranslationExtraction)
    case 'multisystem'
        AT=-[A zeros(size(A));Ap -Bp;App -Bpp];
        bT=-[zeros(4,size(A,3));multiprodMatVec(Cpp,squeeze(bpp))+norm(g)*multiprodMatVec(Cpp,reshape(MEstAligned(:,3),3,[]))];
        N=size(MEstSimilarity,1);
        ord=matStack(reshape(1:N,2,N/6,[]));
        xT=multisolve(AT,bT+reshape(TFactor(ord),size(ord)));
        if flagCollect
            output.AT=AT;
            output.bT=bT;
            output.ord=ord;
            output.xT=xT;
        end
    case 'singlesystem'
        AMat=[num2blkdiag(-A) sparse(2*NFrames,3*NFrames)];
        ApMat=[num2blkdiag(-Ap) num2blkdiag(Bp)];
        AppMat=[num2blkdiag(-App) num2blkdiag(Bpp)];
        %coefficient matrix from projection and derivatives
        ATVec=[AMat;ApMat;AppMat];
        %coefficient matrix for global gravity direction
        AgVec=[sparse(4*NFrames,3);norm(g)*num2blkdiag(Cpp)*MEstAligned];
        %known term from projection and measurements
        bTVec=TFactor+[sparse(4*NFrames,1);vec(-multiprod(Cpp,bpp))];
        if flagDerivativeConstraint
            filterWindow=length(filterCoefficients);
            filterHalfWindow=floor(filterWindow/2);
            MDer1=spconvmtx(filterCoefficients,NFrames,'same');
            ATDer=[kron(MDer1,speye(3))*num2blkdiag(invR(RbsEst)) -num2blkdiag(invR(RbsEst))];
            idxValid=1+filterHalfWindow:NFrames-filterHalfWindow;
            idxValidMat=subMatrix(reshape(1:3*NFrames,3,[]),1:3,idxValid);
            ATDderValid=ATDer(idxValidMat,:);
            ATVec=[ATVec; derivativeConstraintsWeight*ATDderValid];
            AgVec=[AgVec; sparse(size(ATDderValid,1),3)];
            bTVec=[bTVec;zeros(size(ATDderValid,1),1)];
            if flagCollect
                output.filterCoefficients=filterCoefficients;
            end
        end
        
        if flagSecondDerivativeConstraint
            filterWindow=length(secondFilterCoefficients);
            filterHalfWindow=floor(filterWindow/2);
            MDer2=spconvmtx(secondFilterCoefficients,NFrames,'same');
            ATDder=[sparse(3*NFrames,3*NFrames) kron(MDer2,speye(3))*num2blkdiag(invR(RbsEst))];
            AgDder=-num2blkdiag(invR(RbsEst))*norm(gravityVector)*MEstAligned;
            bTDder=vec(multiprodMatVec(invR(RbsEst),alphaIMU));
            idxValid=1+filterHalfWindow:NFrames-filterHalfWindow;
            idxValidMat=subMatrix(reshape(1:3*NFrames,3,[]),1:3,idxValid);
            ATDderValid=ATDder(idxValidMat,:);
            AgDderValid=AgDder(idxValidMat,:);
            bTDderValid=bTDder(idxValidMat,:);
            ATVec=[ATVec; secondDerivativeConstraintsWeight*ATDderValid];
            AgVec=[AgVec; secondDerivativeConstraintsWeight*AgDderValid];
            bTVec=[bTVec; secondDerivativeConstraintsWeight*bTDderValid];
            if flagCollect
                output.secondFilterCoefficients=secondFilterCoefficients;
            end
        end
        if ~flagEstimateAlignedGravity
            gAligned=[0;0;1];
            xTVec=ATVec\(bTVec-AgVec*gAligned);
        else
            xTVecgAligned=[ATVec AgVec]\bTVec;
            xTVec=xTVecgAligned(1:end-3);
            gAligned=xTVecgAligned(end-2:end);
        end
        ordVec=matStack(reshape(1:length(xTVec),3,[],2));
        xT=xTVec(ordVec);
        if flagCollect
            output.ATVec=ATVec;
            output.bTVec=bTVec;
            output.ordVec=ordVec;
            output.xTVec=xTVec;
            output.xT=xT;
            output.gAligned=gAligned;
        end
        
    otherwise
        error('methodTranslationExtraction not recognized.')
end
if flagCollect
    output.methodTranslationExtraction=methodTranslationExtraction;
end

%Assign outputs
taubEst=xT(1:3,:);
nubEst=xT(4:6,:);

function xT=unscrambleTranslation(xTVec)
ordVec=matStack(reshape(1:length(xTVec),3,[],2));
xT=xTVec(ordVec);

function D=num2blkdiag(A)
ACell=squeeze(num2cell(A,[1 2]));
D=spblkdiag(ACell{:});

function x=multisolve(AM,bM)
Nx=size(AM,3);
d=size(AM,2);
x=zeros(d,Nx);
for ix=1:Nx
    x(:,ix)=AM(:,:,ix)\bM(:,ix);
end

function [filterCoefficients,ivarargin]=parseFilterCoefficients(args,ivarargin)
ivarargin=ivarargin+1;
filterCoefficients=lower(args{ivarargin});
if ischar(filterCoefficients)
    switch filterCoefficients
        case 'sgolay'
            ivarargin=ivarargin+1;
            filterOrder=lower(args{ivarargin});
            ivarargin=ivarargin+1;
            filterWindow=lower(args{ivarargin});
            ivarargin=ivarargin+1;
            TSampling=lower(args{ivarargin});
            [~,filterCoefficients] = sgolay(filterOrder,filterWindow);   
            filterCoefficients=-filterCoefficients(:,2)/TSampling;
        otherwise
            error('Filter option for numerical derivative constraint not recognized');
    end
end
