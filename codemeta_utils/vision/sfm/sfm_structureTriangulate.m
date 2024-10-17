function data=sfm_structureTriangulate(data,varargin)
flagDisplayStats=false;
%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'displaystats'
            flagDisplayStats=true;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

NPoints=length(data.structure.featureCount);
location=zeros(3,NPoints);
residual=zeros(1,NPoints);
if flagDisplayStats
    w=getTextWaitBar(NPoints);
    w(0)
end

for iPoint=1:NPoints
    locationImg=sfm_getFeatureLocationByStructureId(data,iPoint);
    projection=sfm_getProjectionFromStructureId(data,iPoint);
    [location(:,iPoint),residual(iPoint)]=triangulate_nonlin(permute(locationImg,[1 3 2]),projection);
    if flagDisplayStats
        w(iPoint)
    end
end
data.structure.location=location;
data.structure.residual=residual;
