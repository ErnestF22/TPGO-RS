%Infer outlier label from cycle closure errors
%function flagE=sfm_rawRotationsCycleErrors2EdgeLabels(C,e,varargin)
%This is the convex optimization plus branch and bound method described in
%   Zach, Klopschitz, Pollefeys
%   "Disambiguating Visual Relations Using Loop Constraints"
%   CVPR 2010
%Inputs
%   C   matrix with cycle vectors
%   e   cycle errors
%Optional Inputs
%   'weightingPrior',s  Method to weight the parameters for the prior
%       'uniform'           Standard unitary weights
%       'count'             Multiply by number of times edge appears in cycles
%       'countNormalized'   As previous but divide by total (weights are
%                           between zero and one)
function flagE=sfm_rawRotationsCycleErrors2EdgeLabels(C,e,varargin)
methodWeightingPrior='uniform';
tolRounding=0.001;
epsLog=eps; %this is a regularization constant to avoid computing log(0)

%probabilities for cycle errors given cycle labels
sigmaInliners=2*pi/180;     %variance for distribution of all-inliers loop
pErrorInliers=@(d,l) exp(-d/sigmaInliners)/(sigmaInliners*(1-exp(-pi/sigmaInliners)));
pErrorOutliers=@(d,l) (1-exp(-d/sigmaInliners))/(pi-sigmaInliners*(1-exp(-pi/sigmaInliners)));

%prior probabilities for each label 
NEdges=size(C,1);
pLabelsPrior=ones(NEdges,1)*[0.9 0.1];

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'perrorinliers'
            ivarargin=ivarargin+1;
            pErrorInliers=varargin{ivarargin};
        case 'perroroutliers'
            ivarargin=ivarargin+1;
            pErrorOutliers=varargin{ivarargin};
        case 'plabelsprior'
            ivarargin=ivarargin+1;
            pLabelsPrior=varargin{ivarargin};
        case 'tolrounding'
            ivarargin=ivarargin+1;
            tolRounding=varargin{ivarargin};
        case 'weightingprior'
            ivarargin=ivarargin+1;
            methodWeightingPrior=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

%set up cost parameters and get first solution
rhoe=-log(pLabelsPrior+epsLog)*[-1;1];
switch lower(methodWeightingPrior)
    case 'uniform'
        %nothing to do
    case 'count'
        rhoe=rhoe.*sum(abs(C),2);
    case 'countnormalized'
        rhoe=rhoe.*sum(abs(C),2)/sum(abs(C(:)));
    otherwise
        error('methodWeightingPrior not recognized')
end
l=sum(abs(C));  
rhol=(-log(pErrorOutliers(e,l)+epsLog))-(-log(pErrorInliers(e,l)+epsLog));
xe=solveLinProgram(rhoe,rhol,C);
xe=roundSolution(xe,tolRounding);

if ~isconverged(xe,tolRounding)
    xe=solveBranch(xe,C,rhoe,rhol,tolRounding);
end

flagE=logical(xe);

%%% Supporting functions %%%
%given a non-converged solution, do branch and bound by fixing a
%non-converged variable at a time
function xe=solveBranch(xe,C,rhoe,rhol,tolSol)
idxFix=getIdxToFix(xe,tolSol);
xeNew=cell(2,1);
vNew=NaN(2,1);

%iterate on possible values for valFix
for valFix=1:2
    xe(idxFix)=valFix-1;
    if isconverged(xe,tolSol)
        %xe(idxFix) was the last unsure variable
        xeNew{valFix}=xe;
    else
        %solve the reduced problem
        xeNew{valFix}=solveReduced(xe,C,rhoe,rhol,tolSol);
        %if is not converged yet, then do the recursive call
        if ~isconverged(xeNew{valFix},tolSol)
            xeNew{valFix}=solveBranch(xe,C,rhoe,rhol,tolSol);
        end
    end
    %evaluate the cost
    xl=computexl(C,xeNew{valFix});
    vNew(valFix)=xeNew{valFix}'*rhoe+xl'*rhol;
end

%pick the lower cost solution
if vNew(1)<vNew(2)
    xe=xeNew{1};
else
    xe=xeNew{2};
end

%solve linear problem only on variables that have not converged yet
function xe=solveReduced(xe,C,rhoe,rhol,tolSol)
NCycles=size(C,2);
dConverged=min(1-xe,xe);

flagFixedEdges=dConverged<tolSol;

%remove cycles which have either
% - all of their edges already fixed 
% - at least one edge set to one
%In both these cases, the value of the non-fixed xe would not affect the
%value of the corresponding xl
flagFixedCycles=any(C & ~flagFixedEdges*ones(1,NCycles)) ...
    | any(C(flagFixedEdges,:) & xe(flagFixedEdges)*ones(1,NCycles));

CReduced=C(~flagFixedEdges,~flagFixedCycles);
rholReduced=rhol(~flagFixedCycles);
rhoeReduced=rhoe(~flagFixedEdges);

xeReduced=solveLinProgram(rhoeReduced,rholReduced,CReduced);

xe(~flagFixedEdges)=xeReduced;
xe=roundSolution(xe,tolSol);

function xe=roundSolution(xe,tolSol)
dConverged=min(1-xe,xe);
flagConverged=dConverged<tolSol;
xe(flagConverged)=round(xe(flagConverged));

function flag=isconverged(xe,tolSol)
flag=all((min(1-xe,xe))<tolSol);

function c=countToFix(xe,tolSol)
c=sum((min(1-xe,xe))>tolSol);

function idx=getIdxToFix(xe,tolSol)
dConverged=min(1-xe,xe);
idxNotConverged=find(dConverged>tolSol);
[~,idxIdx]=min(dConverged(idxNotConverged));
idx=idxNotConverged(idxIdx);

function xl=computexl(C,xe)
NCycles=size(C,2);
xl=zeros(NCycles,1);
for il=1:NCycles
    xl(il)=max(xe(C(:,il)>0));
end

function xe=solveLinProgram(rhoe,rhol,C)
cvx_quiet(true)
cvx_begin
    variable xe(length(rhoe))
    variable xl(length(rhol))
    minimize (xe'*rhoe+xl'*rhol)
    subject to
        xe>=0;
        xe<=1;
        xl>=0;
        xl<=1;
        for il=1:size(C,2)
            c=C(:,il)>0;
            xl(il)>=xe(c);
            xl(il)<=sum(xe(c));
        end
cvx_end


