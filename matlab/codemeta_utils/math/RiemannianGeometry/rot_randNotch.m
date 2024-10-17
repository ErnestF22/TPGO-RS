%Generates random rotations with distances in a given interval
%function rot_randNotch(R,v,N)
%Inputs
%   R   base rotation
%   v   vector with start and end distances
%   N   number of rotations to generate
%Optional inputs
%   'U',U   square root of covariance matrix to generate vectors in the
%           tangent space
function RNoise=rot_randNotch(R,v,N,varargin)
U=eye(size(R(:,:,1)));
if ~exist('R','var') || isempty(R)
    R=eye(3);
end
if ~exist('v','var') || isempty(v)
    v=[45 180]*pi/180;
end
if ~exist('N','var') || isempty(N)
    N=size(R,3);
end

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'u'
            ivarargin=ivarargin+1;
            U=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

d=rot_dim(R);
u=(ones(d,1)*(v(1)+rand(1,N)*(v(2)-v(1)))).*cnormalize(U*randn(d,N));
RNoise=rot_exp(R,rot_hat(R,u));
