% Script that runs ADMM
% =========================================================================

function admm_averagingKmeans()
flagDebug=false;
resetRands(1)
%%
%initializing all variables:
rho = 2;
% array of edges
E=[1 2];%;2 3;3 4];%;5 3; 5 4];
E=[E;fliplr(E)];
%Generalize to each node having a different number of points
d=2;
xiTemplate=[rand(d,10) rand(d,10)+[1;0] rand(d,10)+1];
xi0={xiTemplate,xiTemplate+[3;0],xiTemplate+[0;3],xiTemplate+3};


%Implement also a version where you have K copies of mu at every node, i.e.
%nodes(iNode).mu is a [d x K] matrix
K=3;


d=size(xi0{1},1);
nbNodes = max(vec(E));
if length(xi0)<nbNodes
    error('Number of initial values does not match maximum index in edges')
end
xi0=xi0(1:nbNodes);
nbIterations = 60;
nodes=repmat(struct('xi',[],'mu',[],'z',[],'lambda',[],'neighbors',[]),nbNodes,1);

%initialization loop:
for iNode=1:nbNodes
    nodes(iNode).xi=xi0{iNode};
    nodes(iNode).neighbors= E(E(:,1)==iNode,2);
    nbNeighbors = length(nodes(iNode).neighbors);
    nodes(iNode).mu = repmat(mean(nodes(iNode).xi,2),1,K);
    nodes(iNode).mu = permute(nodes(iNode).mu,[1 3 2]);
    nodes(iNode).z = zeros(d,nbNeighbors,K); 
    nodes(iNode).lambda = zeros(d,nbNeighbors,2,K); 
end
disp(length(xi0));
allMu=NaN([nbIterations d nbNodes K]);
allMu(1,:,:,:)=[nodes.mu]; 

for it=1:nbIterations
    %Loop to update mu's
    for iNode=1:nbNodes
        [ziNode,lambdaiNode]=muData(nodes,iNode);
        %At the first iteration, ziNode should be empty so that the nodes
        %get initialized with pure k-means
        if it==1
            ziNode=NaN(d,0);
        end
        muInit=squeeze(nodes(iNode).mu);
        muUpdated = kmeans_augmented2(nodes(iNode).xi,K,lambdaiNode,ziNode,rho,muInit); %update mu
        nodes(iNode).mu=permute(muUpdated,[1 3 2]);
        %nodes(iNode).mu(:,2)=muUpdate(nodes(iNode).xi,ziNode(:,:,2),lambdaiNode(:,:,2),rho);
        if flagDebug
            admm_averagingKmeans_plot(nodes,xi0)
            %fprintf('mu_{%d}: ziNode=[%s], lambdaiNode=[%s], muNew=[%s]\n',...
            %    iNode,vecString(ziNode),vecString(lambdaiNode),vecString(nodes(iNode).mu))
        end
    end
    
    %Loop to update z's
    for iNode=1:nbNodes
        nbNeighbors=length(nodes(iNode).neighbors);
        for idxJNode=1:nbNeighbors
            jNode=nodes(iNode).neighbors(idxJNode);
            [muiNode,lambdaiNode]=zData(nodes,iNode,jNode);
            nodes(iNode).z(:,idxJNode,:)=zUpdate(muiNode,lambdaiNode,rho); 
            %nodes(iNode).z(:,idxJNode,2)=zUpdate(muiNode(:,3:4),lambdaiNode(:,:,:,2),rho);
            if flagDebug
                fprintf('z_{%d,%d}: muiNode=[%s], lambdaiNode=[%s], zNew=[%s]\n',...
                    iNode,jNode,vecString(muiNode),vecString(lambdaiNode),vecString(nodes(iNode).z(:,idxJNode)))
            end
        end
    end    
    %Loop to update lambda
    for iNode=1:nbNodes
        nbNeighbors=length(nodes(iNode).neighbors);
        for idxJNode=1:nbNeighbors
            jNode=nodes(iNode).neighbors(idxJNode);
            [muiNode,ziNode,lambdaiNode]=lambdaData(nodes,iNode,jNode);
            %at this point:
            %   muiNode [d x 1 x K]
            %   ziNode  [d x direction x K]
            %   lambdaUpdated [d x 1 (neighbor) x directions x K]
            lambdaUpdated=lambdaUpdate(muiNode, ziNode,lambdaiNode, rho);%update lambda
            nodes(iNode).lambda(:,idxJNode,:,:)=lambdaUpdated;
            if flagDebug
                admm_averagingKmeans_plot(nodes,xi0,allMu)
            end
       end
    end
    allMu(it+1,:,:,:) = shiftdim([nodes.mu]);
end
plotPoints([xi0{:}])
hold on
switch d
    case 2
        plot(squeeze(allMu(:,1,:)),squeeze(allMu(:,2,:)))
    case 3
        plot3(squeeze(allMu(:,1,:)),squeeze(allMu(:,2,:)),squeeze(allMu(:,3,:)))
    otherwise
        error('Cannot display this number of dimensions')
end
for iNode=1:nbNodes
    plotPoints(nodes(iNode).mu,'x')
end
hold off
disp('Final mu''s')
disp(squeeze([nodes.mu]));
end
%%Get local information for a node's update
function [ziNode,lambdaiNode]=muData(nodes,iNode)
% ziNode [d x 2*nbNeighbors x K]
% lambdaiNode [d x 2*nbNeighbors x K]
neighbors=nodes(iNode).neighbors;
nbNeighbors=length(neighbors);
[d,K]=size(nodes(iNode).mu);
ziNode=zeros(d,nbNeighbors,K); 
for idxJNode=1:nbNeighbors
    jNode=neighbors(idxJNode);
    %Find index of iNode in the neighbors of jNode
    flagINodeInJ = nodes(jNode).neighbors==iNode;
    ziNode(:,idxJNode,:,:) = nodes(jNode).z(:,flagINodeInJ,:,:);
end
ziNode=[ziNode nodes(iNode).z];
lambdaiNode=reshape(nodes(iNode).lambda,d,[],K);
end

function [muiNode,lambdaiNode]=zData(nodes,iNode,jNode)
%muiNode contains the mu's for an edge
%Expected dimensions: d x 2 (Directions) x K
muiNode=[nodes(iNode).mu,...
    nodes(jNode).mu];
%Find index of iNode in the neighbors of jNode
flagJNodeInI = nodes(iNode).neighbors==jNode;
flagINodeInJ = nodes(jNode).neighbors==iNode;
%lambdaiNode contains the lambda's necessary for updating z from nodes i and j
%Expected dimensions: d x 2 x K
lambdaiNode=[permute(nodes(iNode).lambda(:,flagJNodeInI,1,:),[1 2 4 3])...
    permute(nodes(jNode).lambda(:,flagINodeInJ,2,:),[1 2 4 3])];
end

function [muiNode,ziNode,lambdaiNode]=lambdaData(nodes,iNode,jNode)
muiNode=nodes(iNode).mu;
%Find index of iNode in the neighbors of jNode
flagJNodeInI = nodes(iNode).neighbors==jNode;
flagINodeInJ = nodes(jNode).neighbors==iNode;
ziNode=[nodes(iNode).z(:,flagJNodeInI,:) nodes(jNode).z(:,flagINodeInJ,:)];
lambdaiNode=nodes(iNode).lambda;
end

%%
% % mu equation function:
% function mu=muUpdate(xi, z, lambda, rho)
% nbZ=size(z,2);
% nbXi=size(xi,2); %M_i in the paper
% mu=(2*sum(xi,2)+ sum(lambda,2) + 2*rho*sum(z,2))/(2*nbZ*rho + 2*nbXi);
% end
%%
% z equation function:
function z=zUpdate(mu, lambda, rho)
%update z's for a single edge.
%Expected size: d x 1 x K
z =-sum(lambda,2)/(4*rho)+1/2*(sum(mu,2));
end
%%
% lambda equation function:
function lambda=lambdaUpdate(mu, z, lambda, rho)
%at this point:
%   mu [d x 1 x K]
%   z  [d x direction x K]
%   lambda [d x 1 (neighbor) x directions x K]
z=permute(z,[1 4 2 3]);
mu=permute(repmat(mu,[1 2 1]),[1 4 2 3]);
lambda= lambda + (1/rho)*(reshape(z,size(lambda)) - mu); %mu(:,:,1:2)
end

function [Xv] = vec(X)
[a,b] = size(X);
Xv = reshape(X,a*b,1);
end

%%
function s=vecString(v)
s=num2str(vec(v)',' %d');
end


