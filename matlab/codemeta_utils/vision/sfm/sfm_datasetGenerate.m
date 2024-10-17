function data=sfm_datasetGenerate(varargin)
sigmaNoise=0;
NCameras=7;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'sigmanoise'
            ivarargin=ivarargin+1;
            sigmaNoise=varargin{ivarargin};
        case 'ncameras'
            ivarargin=ivarargin+1;
            NCameras=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

E=adj2edges(adjgallery(NCameras,'kneigh',2),'undirected');

[G,X]=testNetworkCreateAbsolutePoses(NCameras,'references','NPoints',40);

NPoints=size(X,2);
x=projectFromG(G,X,'references');
if sigmaNoise>0
    x=x+sigmaNoise*randn(size(x));
end

data.imgIdName=char2cell(num2str((1:NCameras)','%04d'),[],true);
data.imgSize=2*ones(2,NCameras);
data.resizeFactor=1;
data.imgFileName=cell(1,NCameras);

feature=repmat(struct('locationNormalized',[]),1,NCameras);
for iCamera=1:NCameras
    feature(iCamera).locationNormalized=x(:,:,iCamera);
    feature(iCamera).location=feature(iCamera).locationNormalized;
end

data.feature=feature;
data.matchIdxImg=E';

NEdges=size(E,1);
match=repmat(struct('idxImg',[],'idxMatch',[1;1]*(1:NPoints),'scores',ones(1,NPoints)),1,NEdges);
for iEdge=1:NEdges
    match(iEdge).idxImg=E(iEdge,:)';
end
data.match=match;

data=sfm_poseAdd(data,G);
