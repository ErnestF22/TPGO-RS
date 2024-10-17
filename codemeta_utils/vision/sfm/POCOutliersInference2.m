function POCOutliersInference2
%Part 2: implement CVPR10 paper
%"Disambiguating Visual Relations Using Loop Constraints"

resetRands()

%dataset
load sfm_test_data_synthetic_clean.mat
E=sfm_matchEdges(data,'memberMatch','matchFiltered');
w=sfm_matchCounts(data,'memberMatch','matchFiltered');
w=max(w)-w;
NEdges=size(E,1);
R=G2R(data.matchPoseEstimated);
xeTruth=zeros(NEdges,1);

%generate outliers
NOutliers=5;
v=[45 90]*pi/180;
idxOutliers=randperm(NEdges,NOutliers);
xeTruth(idxOutliers)=1;
R(:,:,idxOutliers)=rot_randNotch(R(:,:,idxOutliers),v,[],'U',diag([1 1 0]));

%sample cycles
[C,e]=sfm_rawRotationsSampleCycleErrors(R,[E w],'triangles','NSamples',100);
NCycles=size(C,2);
xlTruth=zeros(NCycles,1);
for iCycles=1:NCycles
    c=C(:,iCycles)>0;
    xlTruth(iCycles)=max(xeTruth(c));
end

%inference
tolSol=0.001;
flagDisplayDistributions=true;
%This is part 2: implementation of the CVPR10 paper
%"Disambiguating Visual Relations Using Loop Constraints"
sigmaInliners=1.5*pi/180;     %variance for distribution of all-inliers loop
pErrorInliers=@(d,l) exp(-d/sigmaInliners)/(sigmaInliners*(1-exp(-pi/sigmaInliners)));
pErrorOutliers=@(d,l) (1-exp(-d/sigmaInliners))/(pi-sigmaInliners*(1-exp(-pi/sigmaInliners)));
%Check that distributions integrate to one
fIn=@(x) pErrorInliers(x*pi/180,3);
fOut=@(x) pErrorOutliers(x*pi/180,3);
disp(integral(fIn,0,pi))
disp(integral(fOut,0,pi))
if flagDisplayDistributions
    funPlot(fIn,linspace(0,180))
    hold on
    funPlot(fOut,linspace(0,180))
    hold off
end

NEdges=size(C,1);

%prior probabilities for each label 
pLabelsPrior=ones(NEdges,1)*[0.9 0.1];

%set up cost parameters and get first solution
rhoe=-log(pLabelsPrior)*[-1;1];
l=sum(abs(C));  
rhol=(-log(pErrorOutliers(e,l)))-(-log(pErrorInliers(e,l)));
xe=solveLinProgram(rhoe,rhol,C);
xeInit=xe;
if ~isconverged(xe,tolSol)
    xe=solveBranch(xe,C,rhoe,rhol,tolSol);
end

disp([xeInit xe xeTruth])

%given a non-converged solution, do branch and bound by fixing a
%non-converged variable at a time
function xe=solveBranch(xe,C,rhoe,rhol,tolSol)
idxFix=getIdxToFix(xe,tolSol);
disp(['solveBranch: ' num2str(countToFix(xe,tolSol)) ' values to fix'])
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
