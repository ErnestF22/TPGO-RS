%function t_node=testNetworkAddGeoTag(t_node,geoTagIdx, varargin)
%Add geotags, i.e., the projection of the camera centers on the plane Z=0
%for the cameras with indeces geoTagIdx
%
%Optional arguments
%   'random'        Instead of a set of indeces, geoTagIdx must be a
%                   scalar, and a set of geoTagIdx indeces is generated at
%                   random
%   'references'
%   'poses'
%   'MethodAbsolutePoses', method
%   'given',geoTagT     Use the given geotags
function t_node=testNetworkAddGeoTag(t_node,geoTagIdx, varargin)
methodAbsolutePoses='reference';
flagGenerateRandomIdx=false;
flagGiven=false;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'references'
            methodAbsolutePoses='reference';
        case 'poses'
            methodAbsolutePoses='pose';
        case 'methodabsoluteposes'
            ivarargin=ivarargin+1;
            methodAbsolutePoses=lower(varargin{ivarargin});
        case 'random'
            flagGenerateRandomIdx=true;
            if ~isscalar(geoTagIdx)
                error('MATLAB:wrongDimension', 'geoTagIdx must be a scalar when the optional parameter ''random'' is given')
            end
        case 'given'
            flagGiven=true;
            ivarargin=ivarargin+1;
            givenGeoTagT=lower(varargin{ivarargin});
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

structType=testNetworkDetectStructType(t_node);
NNodes=testNetworkGetNumberOfNodes(t_node);

if flagGenerateRandomIdx
    randIdx=randperm(NNodes);
    geoTagIdx=randIdx(1:geoTagIdx);
end

if ~flagGiven
    switch methodAbsolutePoses
        case 'pose'
            %NOTE: since geotags are in world coordinates, we need to invert
            %the matrix gi so that we get the camera center in world
            %coordinates too
            [Ri,geoTagT]=testNetworkGetRotTransl(t_node,'fieldName','gitruth','flagInvertG',true);
        case 'reference'
            [Ri,geoTagT]=testNetworkGetRotTransl(t_node,'fieldName','gitruth','flagInvertG',false);
        otherwise
            error('MATLAB:methodPosesInvalid','methodAbsolutePoses invalid');
    end
else
    geoTagT=zeros(size(givenGeoTagT,1),NNodes);
    geoTagT(:,geoTagIdx)=givenGeoTagT;
end

%geotags have zero Z coordinate
geoTagT(3,:)=0;

switch structType
    case 'single'
        t_node.flagHasGeoTag=false(1,NNodes);
        t_node.flagHasGeoTag(geoTagIdx)=true;
        %set to zero geotag coordinates if we do not have geotags
        geoTagT(:,~t_node.flagHasGeoTag)=0;
        t_node.geoTag=geoTagT;
    case 'array'
        [t_node.flagHasGeoTag]=deal(false);
        [t_node(geoTagIdx).flagHasGeoTag]=deal(true);
        [t_node.geoTag]=deal(zeros(3,1));
        for iNode=geoTagIdx
            t_node(iNode).geoTag=geoTagT(:,iNode);
        end
end
