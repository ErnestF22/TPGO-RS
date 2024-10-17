%function [rotErr,translErr,scaleRatios]=testNetworkComputeErrors(t_node,varargin)
%Compute errors on relative rotations, translations and ratios of scale by comparing
%the relative poses obtained from the current state with the ground truth.
%IF THE GROUND TRUTH IS NOT AVAILABLE, compute the errors with respect the
%measurements. The same can be achieved with the option 'toMeasurements'
%Optional arguments
%   'Measurements'      compute errors using the pairwise measurements instead of
%                       the estimated poses
%   'toMeasurements'    compute errors by comparing to the measurements
%                       instead of the ground truth (it is the default if
%                       the ground truth is not available)
%
%The rotation and translation errors are given in radians
%NOTE: if both optional arguments 'Measurements' and 'toMeasurements' are
%active, then, obviously, all the computed errors will be zero.

%%AUTORIGHTS%%

function [rotErr,translErr,scaleRatios,translErrNorm]=testNetworkComputeErrors(t_node,varargin)
flagFromMeasurements=false;
flagToMeasurements=~testNetworkHasGroundTruth(t_node);
methodAbsolutePoses='reference';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'measurements'
            flagFromMeasurements=true;
        case 'tomeasurements'
            flagToMeasurements=true;
        case 'references'
            methodAbsolutePoses='reference';
        case 'poses'
            methodAbsolutePoses='pose';
        case 'methodabsoluteposes'
            ivarargin=ivarargin+1;
            methodAbsolutePoses=varargin{ivarargin};
       otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

structType=testNetworkDetectStructType(t_node);
N=testNetworkGetNumberOfNodes(t_node);

%allocate and initialize
switch structType
    case 'single'
        NEdges=t_node.NEdges;
        t_node.rotErr=NaN(NEdges,1);
        t_node.translErr=NaN(NEdges,1);
        t_node.translErrNorm=NaN(NEdges,1);
        t_node.scaleRatios=NaN(NEdges,1);

        t_node=testNetworkExecFunOnEdges(t_node,...
            @(t,i) edgeErrorSingle(t,i,flagFromMeasurements,flagToMeasurements,methodAbsolutePoses));

        rotErr=t_node.rotErr;
        translErr=t_node.translErr;
        translErrNorm=t_node.translErrNorm;
        scaleRatios=t_node.scaleRatios;
        
    case 'array'
        [t_node.rotErr]=deal(NaN(1,N));
        [t_node.translErr]=deal(NaN(1,N));
        [t_node.translErrNorm]=deal(NaN(1,N));
        [t_node.scaleRatios]=deal(NaN(1,N));
        
        E=testNetworkGetEdges(t_node);
        t_node=testNetworkExecFunOnEdges(t_node,...
            @(t_node,iNode,jNode) edgeErrorArray(t_node,iNode,jNode,flagFromMeasurements,flagToMeasurements,methodAbsolutePoses),E);

        A=cat(1,t_node.aij);
        rotErrMat=cat(1,t_node.rotErr);
        translErrMat=cat(1,t_node.translErr);
        translErrNormMat=cat(1,t_node.translErrNorm);
        scaleRatiosMat=cat(1,t_node.scaleRatios);

        idx=find(A~=0);
        rotErr=rotErrMat(idx);
        translErr=translErrMat(idx);
        translErrNorm=translErrNormMat(idx);
        scaleRatios=scaleRatiosMat(idx);
end


function t_node=edgeErrorSingle(t_node,iEdge,flagFromMeasurements,flagToMeasurements,methodAbsolutePoses)
if flagFromMeasurements
    gij=t_node.gij(:,:,iEdge);
else
    iNode=t_node.E(iEdge,1);
    jNode=t_node.E(iEdge,2);
    
    gi=t_node.gi(:,:,iNode);
    gj=t_node.gi(:,:,jNode);
    gij=computeRelativePoseFromG(gi,gj,'methodAbsolutePoses',methodAbsolutePoses);
end
if ~flagToMeasurements
    gijtruth=t_node.gijtruth(:,:,iEdge);
else
    gijtruth=t_node.gij(:,:,iEdge);
end

[t_node.rotErr(iEdge),...
    t_node.translErr(iEdge),...
    t_node.scaleRatios(iEdge),...
    t_node.translErrNorm(iEdge)]=computeEdgeErr(gij,gijtruth);


function t_node=edgeErrorArray(t_node,iNode,jNode,flagFromMeasurements,flagToMeasurements,methodAbsolutePoses)
if flagFromMeasurements
    gij=t_node(iNode).gij(:,:,jNode);
else
    gi=t_node(iNode).gi;
    gj=t_node(jNode).gi;
    gij=computeRelativePoseFromG(gi,gj,'methodAbsolutePoses',methodAbsolutePoses);
end
if ~flagToMeasurements
    gijtruth=t_node(iNode).gijtruth(:,:,jNode);
else
    gijtruth=t_node(iNode).gij(:,:,jNode);
end

[t_node(iNode).rotErr(jNode),...
    t_node(iNode).translErr(jNode),...
    t_node(iNode).scaleRatios(jNode),...
    t_node(iNode).translErrNorm(jNode)]=computeEdgeErr(gij,gijtruth);


% function gij=computeRelPose(gi,gj,methodAbsolutePoses)
% switch methodAbsolutePoses
%     case 'reference'
%         gij=invg(gi)*gj;
%     case 'pose'
%         gij=gi*invg(gj);
%     otherwise
%         error('Absolute pose method specification not valid')
% end

function [rotErr,translErr,scaleRatio,translErrNorm]=computeEdgeErr(gij,gijtruth)
Rij=gij(1:3,1:3);
Rijtruth=gijtruth(1:3,1:3);
[tij,lambdaij]=cnormalize(gij(1:3,4));
[tijtruth,lambdaijtruth]=cnormalize(gijtruth(1:3,4));

rotErr=rot_dist(Rij,Rijtruth);
translErr=acos(max(-1,min(1,tij'*tijtruth)));
scaleRatio=lambdaij/lambdaijtruth;
translErrNorm=norm(tij-tijtruth);
