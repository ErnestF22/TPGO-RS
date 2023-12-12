%Returns the gyroscopic term in Euler's rotation equations
%function dw=rotDyn_gyroscopicTerm(w,varargin)
%dw=hat(w)*J*w

function dw=rotDyn_gyroscopicTerm(w,varargin)
J=eye(3);

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'inertiamatrix'
            ivarargin=ivarargin+1;
            J=varargin{ivarargin};
            flagCancelDynamics=true;
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

dw=multiprodMatVec(hat3(w),J*w);
