%Generates points that are regularly spaced in polar coordinates
%function [x,r,theta]=polarGrid(Nr,Ntheta,varargin)
%Inputs
%   Nr      Number of directions for the radial component
%   Ntheta  Number of directions for the angular component
%Optional inputs
%   'rMin',r        Set minimum value for radial component
%   'rMax',r        Set maximum value for radial component
%   'thetaMin',r    Set minimum value for angular component
%   'thetaMax',r    Set maximum value for angular component

function [x,r,theta]=polarGrid(Nr,Ntheta,varargin)
rMin=1;
rMax=5;
thetaMin=-pi;
thetaMax=(Ntheta-2)/Ntheta*pi;
%flatten x into a [2 x Nr*Ntheta] array, 
%otherwise it is a [Nr x Ntheta x 2] array
flagFlatten=true;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'rmin'
            ivarargin=ivarargin+1;
            rMin=varargin{ivarargin};
        case 'rmax'
            ivarargin=ivarargin+1;
            rMax=varargin{ivarargin};
        case 'thetamin'
            ivarargin=ivarargin+1;
            thetaMin=varargin{ivarargin};
        case 'thetamax'
            ivarargin=ivarargin+1;
            thetaMax=varargin{ivarargin};
        case 'noflatten'
            flagFlatten=false;
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

r=linspace(rMin,rMax,Nr);
theta=linspace(thetaMin,thetaMax,Ntheta);
x(:,:,1)=r'*cos(theta);
x(:,:,2)=r'*sin(theta);

if flagFlatten
    x=reshape(permute(x,[3 1 2]),2,Nr*Ntheta);
end