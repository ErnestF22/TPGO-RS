%Normalize the homogeneous vectors defining planes
%function n=planeNormalize(N,varargin)
%Input
%   N   [D+1 x NPlanes] array containing the homogeneous vectors defining
%       the NPlanes in R^D
%Optional Input
%   'unitNormal'    (default) or
%   'unitDistance'  select the method used for normalization
%Output
%   n   [D+1 x NPlanes] array obtained by dividing each vector in N by one
%       of the following:
%       'unitNormal'    norm of N(1:D,:)
%       'unitDistance'  N(D+1,:)
function n=planeNormalize(N,varargin)
method='unitNormal';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'unitnormal'
            method='unitaryNormal';
        case 'unitdistance'
            method='unitDistance';
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

D = size(N,1)-1;
switch lower(method)
    case 'unitnormal'
        normN = sqrt(sum(N(1:D,:).^2,1));
    case 'unitdistance'
        normN = N(D+1,:);
end

idxNotZero=abs(normN)>1e-14;
n=zeros(size(N));
if sum(idxNotZero)>0
    n(:,idxNotZero) =  N(:,idxNotZero) ./ (ones(D+1,1)*normN(idxNotZero));
end
