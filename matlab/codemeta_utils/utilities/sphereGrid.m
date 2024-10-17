%Generate spherical grid
%function [p,NRadii,NLatitudes,NLongitudes]=sphereGrid(varargin)
%Generate a grid on concentric spheres.
%Optional inputs
%   'NRadii',n          number of concentric spheres
%   'NLatutudes',n      number of latitudes on each sphere
%   'NLongitudes',n     number of longitudes on each sphere
%   'Unique'            remove duplicate points. The output is a 2D matrix.
%   'Sorted'            as 'unique' plus sort points using the z coordinate
%Output
%   p   [3 x NRadii x NLatitudes x NLongitudes] array with the 3-D
%       coordinates. If 'Unique' is selected, then is a [3 x NPoints]
%       array.
function [p,NRadii,NLatitudes,NLongitudes]=sphereGrid(varargin)
NRadii=4;
NLatitudes=5;
NLongitudes=8;
stepRadii=1;
flagUnique=false;
flagSort=false;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'nradii'
            ivarargin=ivarargin+1;
            NRadii=varargin{ivarargin};
        case 'nlatitudes'
            ivarargin=ivarargin+1;
            NLatitudes=varargin{ivarargin};
        case 'nlongitudes'
            ivarargin=ivarargin+1;
            NLongitudes=varargin{ivarargin};
        case 'stepradii'
            ivarargin=ivarargin+1;
            stepRadii=varargin{ivarargin};
        case 'unique'
            flagUnique=true;
        case 'sorted'
            flagUnique=true;
            flagSort=true;
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

NRadii=NRadii+1;

p=zeros(3,NRadii,NLatitudes,NLongitudes);
for iRadii=1:NRadii
    r=stepRadii*(iRadii-1);
    for iLongitudes=1:NLongitudes
        theta=2*pi*(iLongitudes-1)/NLongitudes;
        for iLatitudes=1:NLatitudes
            phi=pi*(iLatitudes-1)/(NLatitudes-1);
            p(1,iRadii,iLatitudes,iLongitudes)=r*sin(phi)*cos(theta);
            p(2,iRadii,iLatitudes,iLongitudes)=r*sin(phi)*sin(theta);
            p(3,iRadii,iLatitudes,iLongitudes)=r*cos(phi);
        end
    end
end

if flagUnique
    p=reshape(p,size(p,1),[]);
    p=uniqueApprox(p);
    if flagSort
        [~,idxSort]=sort(p(3,:),'descend');
        p=p(:,idxSort);
    end        
end
