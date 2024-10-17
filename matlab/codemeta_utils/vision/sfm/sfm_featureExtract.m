%Extracts SIFT feature and adds them to the data structure
%function data=sfm_featureExtract(data)
%Get list of images from imgFileName and use vl_sift to extract features.
%This function heuristically adjusts the first octave parameters of vl_sift
%to adapt to the resolution of the image.
%
%Input
%   data    data structure with file information for NImages
%Optional inputs
%   'showStats'     Display information as each image is processed
%   'optsSIFT'      Cell array of additional parameters to pass to vl_sift
%   
%Output
%   data    the same data structure with the new field imgFeature containing
%           an NImages long structure array with the following fields
%       featLocation    [2 x NFeature] array with pixel location
%       featScale       [1 x NFeature] array with scales
%       featAngle       [1 x NFeature] array with angles (in radians)
%       featDescriptor  [128 x NFeature] array with SIFT descriptors
%       featIdName      IdName of the image
%       featIdNumber    number of the feature in the image
function data=sfm_featureExtract(data,varargin)
flagShowStats=false;
flagAutoFirstOctave=true;
NImages=length(data.imgFileName);
optsSift={};

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'showstats'
            flagShowStats=true;
        case 'optssift'
            ivarargin=ivarargin+1;
            optsSift=[optsSift{:} varargin{ivarargin}];
        case 'flagautofirstoctave'
            ivarargin=ivarargin+1;
            flagAutoFirstOctave=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

feature=repmat(struct('location',[],'scale',[],'angle',[],'descriptor',[]),1,NImages);
for iImage=1:NImages
    if flagShowStats
        fprintf('Extracting features image %d\n',iImage);
    end
    im=single(rgb2gray(sfm_getImageById(data,iImage)));
    if flagAutoFirstOctave
        imWidth=data.imgSize(1,iImage);
        siftFirstOctave=round(1.32*log2(imWidth)-12);
        if flagShowStats
            fprintf('\tFirst octave: %d\n',siftFirstOctave);
        end
        optsSiftOctave={'FirstOctave' siftFirstOctave};
    else
        optsSiftOctave={};
    end
    [f,d]=vl_sift(im,optsSiftOctave{:},optsSift{:});
    feature(iImage).location=f(1:2,:);
    feature(iImage).scale=f(3,:);
    feature(iImage).angle=f(4,:);
    feature(iImage).descriptor=d;
end
data.feature=feature;
