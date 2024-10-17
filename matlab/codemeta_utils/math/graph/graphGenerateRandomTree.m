%function E=graphGenerateRandomTree(E)
%Extract a random tree subgraph from the graph given by E
%NOTES:
% - the given graph is assumed to be connected and it is assumed that
%   the edges are not repeated
% - the edges in ETree might be reversed
% - the order of the edges is visit compatible, i.e., the starting endpoint
%   always appears after being a arriving endpoint before (except for the
%   first one)
function [ETree,VTree]=graphGenerateRandomTree(E)
E=[E;fliplr(E)];

V=unique(E);

VTree=V(randi(length(V)));

ETree=[];
VNoTree=setdiff(V,VTree);
EBridge=E(ismember(E(:,1),VTree) & ismember(E(:,2),VNoTree),:);

while ~isempty(EBridge)
    iNewEdge=randi(size(EBridge,1));
    vnew=EBridge(iNewEdge,2);
    
    VTree=[VTree; vnew];
    ETree=[ETree; EBridge(iNewEdge,:)];
    VNoTree(VNoTree==vnew)=[];

    EBridge=E(ismember(E(:,1),VTree) & ismember(E(:,2),VNoTree),:);
end

function c=randi(m)
c=floor(rand*m)+1;
