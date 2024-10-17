function POCOutliersInference
%This is part 1: sampling of the cycles and evaluation of the errors
load sfm_test_data_fountain
data=sfm_matchPoseTruth(data,'memberMatch','matchFiltered');
NSamples=500;
%A=sfm_matchAdjMatrix(data,'member','matchFiltered');
%[E,wInit]=adj2edges(A,'undirected');
E=[data.matchFiltered.idxImg]';
wInit=cellfun(@(x) size(x,2),{data.matchFiltered.idxMatch})';
wInit=max(wInit)-wInit;
C=grCycleSamples([E wInit],'NSamples',NSamples);
NSamples=size(C,2);
R=G2R(data.matchPoseEstimated);
RTruth=G2R(data.matchPoseTruth);
l=zeros(NSamples,1);
e=zeros(NSamples,1);
for iSample=1:NSamples
    c=C(:,iSample);
    [idxc,signc]=edgesSortIndex(E,c);
    RCombined=rotationOrientedCombine(R,idxc,signc);
    l(iSample)=sum(abs(c));
    e(iSample)=rot_dist(RCombined,eye(3))*180/pi;
end
disp(rot_dist(R,RTruth,'vector')'*180/pi);
[lUnique,meanE]=indexedMean(l,e);
plot(l,e,'.',lUnique,meanE,'-')

%keyboard

function [mUnique,meanX]=indexedMean(m,x)
mUnique=unique(m);
meanX=zeros(size(mUnique));
for im=1:length(mUnique)
    meanX(im)=mean(x(m==mUnique(im)));
end

function RCombined=rotationOrientedCombine(R,idxC,signC)
RCombined=eye(3);
for ic=1:length(idxC)
    if signC(ic)>0
        RCombined=RCombined*R(:,:,idxC(ic));
    else
        RCombined=RCombined*R(:,:,idxC(ic))';
    end
end
