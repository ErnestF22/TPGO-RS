function tagsDatasetBuildSynthetic
resetRands()
flagDisplay=true;
flagSave=true;

L=10;   %size of square where tags are positioned
hMin=5;   %height of cameras
hMax=10;

NPoses=10;
NTags=6;

RCamera=[
    1  0  0;
    0 -1  0;
    0  0 -1
    ];
RTags=eye(3);

TTags=[
    L*(rand(2,NTags)-0.5);
    zeros(1,NTags)
    ];

N=NTags+NPoses;
idxCameras=1:NPoses;
idxTags=NPoses+(1:NTags);

methodAbsolutePoses='reference';

G=zeros(4,4,N);
G(4,4,:)=1;
for iCamera=1:NPoses
    G(1:3,1:3,idxCameras(iCamera))=RCamera*rot_randn(eye(3),0.5);
    G(1:3,4,idxCameras(iCamera))=[0;0;0]+[2*randn;2*randn;hMin+(hMax-hMin)*rand];
end
for iTag=1:NTags
    G(1:3,1:3,idxTags(iTag))=RTags;
    G(1:3,4,idxTags(iTag))=TTags(:,iTag);
end

tagSize=2;
tagXFrame=tagSize/2*[1 -1 -1 1; 1 1 -1 -1; 0 0 0 0];
tagX=[];
for iTag=1:NTags
    tagX=cat(3,tagX,rigidTransformG(G(:,:,idxTags(iTag)),tagXFrame,'methodAbsolutePoses',methodAbsolutePoses,'cw'));
end

for iCamera=1:NPoses
    for iTag=1:NTags
        posex{iCamera}(iTag).id=iTag-1;
        posex{iCamera}(iTag).x=projectFromG(G(:,:,idxCameras(iCamera)),tagX(:,:,iTag),'methodAbsolutePoses',methodAbsolutePoses);
    end
end

tagG=G(:,:,idxTags);
poseG=G(:,:,idxCameras);
poseNTags=NTags*ones(NPoses,1);

if flagSave
    save('tagDatasetSynthetic','posex','poseG','NPoses','NTags','poseNTags','tagX','tagG','tagXFrame');
end

if flagDisplay
    figure(1)
    testNetworkDisplay(G,'references')
    hold on
    plotPoints(tagX)
    hold off
    figure(2)
    for iCamera=1:NPoses
        plotPoints([posex{iCamera}.x]);
        hold on
    end
    axis([-1 1 -1 1])
end