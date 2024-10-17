function [constraintList,constraintParametersList,E]=locSegSampleConstraintsGallery(sampleNum,E,x,R,constraintList,constraintParametersList)
if ~exist('constraintList','var')
    constraintList={};
end
if ~exist('constraintParametersList','var')
    constraintParametersList={};
end
NNodes=size(x,2);
NEdges=size(E,1);
disp('Constraints from gallery:')
switch sampleNum
    case 1
        disp('  all fixed bearing')
        [constraintList{end+1:end+NEdges}]=deal('fixedBearings');
        E=[E; (1:NNodes)'*[1 1]];
        [constraintList{end+1:end+NNodes}]=deal('fixedRotations');
        [constraintParametersList{end+1:end+NEdges+NNodes}]=deal({});
    case 2
        disp('  all relative bearings')
        [constraintList{end+1:end+NEdges}]=deal('relativeBearings');
        t=generateRelativeBearingMeasurements(x,R,E);
        for iEdge=1:NEdges
            iNode=E(iEdge,1);
            constraintParametersList{end+1}={R(:,:,iNode), t(:,iEdge)};
        end
    case 3
        disp('  all distances')
        [constraintList{end+1:end+NEdges}]=deal('distances');
        [constraintParametersList{end+1:end+NEdges}]=deal({});
    case 4
        disp('  all relative rotations')
        [constraintList{end+1:end+NEdges}]=deal('relativeRotations');
        for iEdge=1:NEdges
            iNode=E(iEdge,1);
            jNode=E(iEdge,2);
            Ri=R(:,:,iNode);
            Rj=R(:,:,jNode);
            constraintParametersList{end+1}={Ri,Rj,Ri*Rj'};
        end
    case 5
        disp('  all relative translations')
        [constraintList{end+1:end+NEdges}]=deal('relativeTranslations');
        T=generateRelativeTranslationMeasurements(x,R,E);
        for iEdge=1:NEdges
            iNode=E(iEdge,1);
            constraintParametersList{end+1}={R(:,:,iNode), T(:,iEdge)};
        end
    case 6
        disp('  all fixed translations')
        [constraintList{end+1:end+NEdges}]=deal('fixedTranslations');
        [constraintParametersList{end+1:end+NEdges}]=deal({});
        
end

function T=generateRelativeTranslationMeasurements(x,R,E)
dimData=size(x,1);
NEdges=size(E,1);
T=zeros(dimData,NEdges);
for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    xi=x(:,iNode);
    xj=x(:,jNode);
    Tij=xi-xj;
    Ri=R(:,:,iNode);
    T(:,iEdge)=Ri'*Tij;
end

function t=generateRelativeBearingMeasurements(x,R,E)
t=cnormalize(generateRelativeTranslationMeasurements(x,R,E));
