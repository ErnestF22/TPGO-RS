function quickshift_breakTreeMerge_testDistance
treeEdges=[1 1 2 1 4 4 4 7 7];
treeDistances=[0 5 1 3 1 1 3 2 1];
[~,edgesSorted]=sort(treeDistances,'ascend');
rho=1.1;
strategy='max';
fMerge=@(cmpData1,cmpData2,treeData,~) fMergeDistance(cmpData1,cmpData2,treeData,rho,strategy);

treeData=num2cell(treeDistances);

[treeEdgesNew,info]=quickshift_breakTreeMerge(treeEdges,treeData,edgesSorted,fMerge);

disp(edgesSorted)
NEdges=length(treeEdges);
for iEdge=1:NEdges
    fprintf('%2d --(d:%3.1f)--> O:%2d/N:%2d (c:%2d)\n',iEdge,treeDistances(iEdge),treeEdges(iEdge),treeEdgesNew(iEdge),info.treeComponents(iEdge));
end

function [flag, componentDataMerged]=fMergeDistance(cmpDistance1,cmpDistance2,edgeDistance,rho,strategy)
if isempty(cmpDistance1) && isempty(cmpDistance2)
    %if both are singleton components, just join them
    flag=true;
    componentDataMerged=edgeDistance;
else
    %note that the max operation takes care of the case where one of the
    %two components is a singleton
    flag=edgeDistance<rho*max([cmpDistance1 cmpDistance2]);
    if flag
        switch strategy
            case 'max'
                componentDataMerged=max([cmpDistance1,cmpDistance2,edgeDistance]);
            case 'min'
                componentDataMerged=min([cmpDistance1,cmpDistance2,edgeDistance]);
        end
    else
        componentDataMerged=NaN;
    end
end
