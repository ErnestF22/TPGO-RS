%Compute the velocity control from the bearing potential
%function dx=bearingControlDirect(y,yg,funs)
%Inputs
%   y       Current bearing vectors
%   yg      Goal's bearing vectors
%   funs    Reshaping function for the generalized cost
%
%For simulink: hardcode funs, comment out optional arguments, and remove funs
%and varargin inputs
function dx=bearingControlDirect(y,yg,funs,varargin)
flagPrecondition=false;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'precondition'
            flagPrecondition=true;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

if ~flagPrecondition
    dx=-bearingCostGeneral_gradient(y,yg,funs);
else
    dx=bearingCostGeneral_gradient(y,yg,funs);
    [~,Hdx]=bearingCostGeneral_gradient(yg,yg,funs);
    %Hdx=Hdx/det(Hdx);
    dx=-Hdx\dx;
end
    