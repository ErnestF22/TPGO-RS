%Generate a synthetic sfm dataset from the poses G and the 3-D points X
function data=sfm_datasetGenerateFromGX(G,X,varargin)
NPointsOverlap=15;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'npointsoverlap'
            ivarargin=ivarargin+1;
            NPointsOverlap=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

NCameras=size(G,3);
data.imgIdName=char2cell(num2str((1:NCameras)','%04d'),[],true);
data.imgSize=2*ones(2,NCameras);
data.resizeFactor=1;
data.imgFileName=cell(1,NCameras);

x=projectFromG(G,X,'references');
z=projectGetDepthsFromG(G,X,'references');

flagXVisible=squeeze(x(1,:,:)<1 & x(1,:,:)>-1 ...
    &  x(2,:,:)<1 & x(2,:,:)>-1)' ...
    & z>0;
idxXVisible=zeros(size(flagXVisible));
for iCamera=1:NCameras
    flagCamera=flagXVisible(iCamera,:);
    idxXVisible(iCamera,flagCamera)=1:sum(flagCamera);
end

feature=repmat(struct('locationNormalized',[]),1,NCameras);
for iCamera=1:NCameras
    feature(iCamera).locationNormalized=x(:,flagXVisible(iCamera,:),iCamera);
end

match=[];
for iCamera=1:NCameras
    for jCamera=iCamera+1:NCameras
        flagOverlap=flagXVisible(iCamera,:) & flagXVisible(jCamera,:);
        if sum(flagOverlap)>NPointsOverlap
            idxImg=[iCamera;jCamera];
            idxMatch=[idxXVisible(iCamera,flagOverlap); idxXVisible(jCamera,flagOverlap)];
            scores=ones(1,sum(flagOverlap));
            match=[match; struct('idxImg',idxImg,'idxMatch',idxMatch,'scores',scores)];
        end
    end
end

data.feature=feature;
data.matchIdxImg=[match.idxImg];
data.match=match;

data=sfm_poseAdd(data,G);
