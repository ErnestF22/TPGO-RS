function [t_node,errors]=localization_MLE_rigid(t_node,varargin)
flagInitRotationTree=true;          %discart init already present in t_node and use a tree
flagRescaleDispersionMatrices=true; %multiply all the dispersion matrices by a constant so that they have unitary evals in average
flagSaveIntermediateData=true;      %save state with intermediate and final results
flagShowMessages=false;
flagCollectErrors=false;
flagInit=true;
flagHasTruth=isfield(t_node,'gitruth');
flagHasEType=false;
%optsLieMinimize={'disablePhase12'};
optsLieMinimize={};
rotationMaxIt=5000;
totalMaxIt=50000;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'flaginit'
            ivarargin=ivarargin+1;
            flagInit=varargin{ivarargin};
        case 'noinit'
            flagInit=false;
        case 'displayit'
            optsLieMinimize=[optsLieMinimize 'displayIt'];
        case 'optslieminimize'
            ivarargin=ivarargin+1;
            optsLieMinimize=[optsLieMinimize varargin{ivarargin}];
        case 'showmessages'
            flagShowMessages=true;
        case 'flagshowmessages'
            ivarargin=ivarargin+1;
            flagShowMessages=varargin{ivarargin};
        case 'rotationmaxit'
            ivarargin=ivarargin+1;
            rotationMaxIt=varargin{ivarargin};
        case 'totalmaxit'
            ivarargin=ivarargin+1;
            totalMaxIt=varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

if nargout>1
    flagCollectErrors=true;
end

NNodes=testNetworkGetNumberOfNodes(t_node);
if isfield(t_node,'gi')
    [Ri0,Ti0]=testNetworkGetRotTransl(t_node);
else
    Ri0=repmat(eye(3),[1 1 NNodes]);
    Ti0=zeros(3,NNodes);
end
[Rij,Tij]=testNetworkGetRelativeRotTranslScales(t_node);
[GammaijR,GammaijT]=testNetworkGetDispersionMatricesRotationTranslation(t_node);
Gammaij=testNetworkGetDispersionMatricesRotationTranslation(t_node);
if flagRescaleDispersionMatrices
    meanEval=mean(squeeze(cat(3,...
        Gammaij(1,1,:), Gammaij(2,2,:), Gammaij(3,3,:),...
        Gammaij(4,4,:), Gammaij(5,5,:), Gammaij(6,6,:))));
    GammaijR=GammaijR/meanEval;
    GammaijT=GammaijT/meanEval;
    Gammaij=Gammaij/meanEval;
end
E=testNetworkGetEdges(t_node);

if isfield(t_node,'EType')
    flagHasEType=true;
    EType=t_node.EType;
else
    EType=ones(size(E,1),1);
end

ticFlag(flagShowMessages);

if flagHasTruth
    [RiTruth,TiTruth]=G2RT(t_node.gitruth);
end
if flagInit
    if flagHasEType
        error('Initialization with non-default edge types not implemented yet')
    end
    if flagInitRotationTree
        %init R with tree
        fprintfFlag(flagShowMessages,'Initialization of rotations with random spanning tree...')
        ETree=testNetworkGenerateRandomTree(t_node);
        t_node=testNetworkLocalizeTree(t_node,ETree,'References');
        tocFlag(flagShowMessages);
    end

    [RiInit,TiInit]=testNetworkGetRotTransl(t_node);

    fRotations=@(R) rotationLogLikelihoodNetwork(R,Rij,GammaijR,E);
    %gradfRotations=@(R) rot_hat(R,gradRotationLogLikelihoodNetwork(R,Rij,GammaijR,E));

    fprintfFlag(flagShowMessages,'Optimization of rotations alone...');
    %[RiEst,errorsRiEst]=lie_minimize(rot_funs(),fRotations,gradfRotations,RiInit,'stepsize',0.1,'displayIt');
    [RiEst,errorsRiEst]=lie_minimizeGradNewton(rot_funs(),fRotations,RiInit,...
        'maxIt',rotationMaxIt,'epsilon',0.1,optsLieMinimize{:});%,'showCost'
    tocFlag(flagShowMessages);

    if flagShowMessages
        if flagHasTruth
            disp('Cost rotation only RiInit/RiEst/RiTruth')
            disp([fRotations(RiInit) fRotations(RiEst) fRotations(RiTruth)])
        else
            disp('Cost rotation only RiInit/RiEst')
            disp([fRotations(RiInit) fRotations(RiEst)])
        end
    end

    fTranslations=@(R,T) translationLogLikelihoodNetwork(R,T,Tij,GammaijT,E);

    %Linear initialization of the translations
    fprintfFlag(flagShowMessages,'Linear initialization of translations...');
    [A,b]=translationOnlyLogLikelihoodNetworkMatrix(RiEst,Tij,GammaijT,E);
    %A is singular by construction, so we need to use the SVD
    [U,S,~]=svd(A);
    TiInit=reshape(-U(:,1:end-3)*(S(1:end-3,1:end-3)\(U(:,1:end-3)'*b)),3,[]);
    tocFlag(flagShowMessages);

    if flagShowMessages
        if flagHasTruth
            disp('Cost translation only at RiEst TiInit/TiTruth')
            disp([fTranslations(RiEst,TiInit) fTranslations(RiEst,TiTruth)])
        else
            disp('Cost translation only at RiEst TiInit')
            disp(fTranslations(RiEst,TiInit))
        end
    end

    if flagSaveIntermediateData
        save([mfilename '_data'])
    end
    %load([mfilename '_data'])
else
    RiEst=Ri0;
    TiInit=Ti0;
end

%Optimization of the total cost function
% fRotationsTranslations=@(G) logLikelihoodNetwork(G2R(G),G2T(G),Rij,Tij,Gammaij,E);
% gradfRotationsTranslations=@(G) gradLogLikelihoodNetwork_GFormat(G,Rij,Tij,Gammaij,E);
% GiInit=RT2G(RiEst, TiInit);
% 
% fprintfFlag(flagShowMessages,'Optimization of the total cost function...');
% [GiEst,errorsGiEst]=lie_minimize(rot3r3_funs(),fRotationsTranslations,gradfRotationsTranslations,GiInit,...
%     'stepsize',0.01,'maxit',10000); %'showCost',
fRotationsTranslations=@(G) logLikelihoodNetwork(G2R(G),G2T(G),Rij,Tij,Gammaij,E,EType);
GiInit=RT2G(RiEst, TiInit);

fprintfFlag(flagShowMessages,'Optimization of the total cost function...');
[GiEst,errorsGiEst]=lie_minimizeGradNewton(rot3r3_funs(),fRotationsTranslations,GiInit,...
    'epsilon',0.1,'MaxIt',totalMaxIt,optsLieMinimize{:}); %,'showCost'
tocFlag(flagShowMessages);

RiEst2=G2R(GiEst);
TiEst=G2T(GiEst);
if flagShowMessages
    if flagHasTruth
        disp('Cost total (RiEst,TiInit)/(RiEst,TiEst)/(RiTruth/TiTruth)')
        disp([fRotationsTranslations(RT2G(RiEst,TiInit))...
            fRotationsTranslations(RT2G(RiEst2,TiEst))...
            fRotationsTranslations(RT2G(RiTruth,TiTruth))])
    else
        disp('Cost total (RiEst,TiInit)/(RiEst,TiEst)')
        disp([fRotationsTranslations(RT2G(RiEst,TiInit))...
            fRotationsTranslations(RT2G(RiEst2,TiEst))])
    end
end

t_node.gi=GiEst;

if flagCollectErrors
    errors.errorsGiEst=errorsGiEst;
    if flagInit
        errors.errorsRiEst=errorsRiEst;
        errors.RiInit=RiInit;
    end
    errors.RiEst=RiEst;
    errors.TiInit=TiInit;
end

if flagSaveIntermediateData
    save([mfilename '_data2'])
end

% function v=computeGradfRotationsTranslations(R,T,Rij,Tij,GammaijR,GammaijT,E)
% vR1=-rot_hat(R,gradRotationLogLikelihoodNetwork(R,Rij,GammaijR,E));
% [vR2,vT]=gradTranslationLogLikelihoodNetwork(R,T,Tij,GammaijT,E);
% vR2=-rot_hat(R,vR2);
% v=[vR1+vR2 permute(vT,[1 3 2]); zeros(1,4,size(vR2,3))];
