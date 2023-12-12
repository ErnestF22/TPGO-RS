function [kd,kv,kp,beta] = TR3xR3_fminsearch(e_A, e_B, varargin)
% Use derivative free method Nelder and Mead Algorithm and find optimal
% gains kd, kv, kp. This method is a builtin as the fminsearch matlab
% function.
% INPUTS:
%   e_A := eigenvalue range for the A matrix corresponding to rho(x,xd)
%      [min;max]
%   e_B := eigenvalue range for the B matrix correspondding to rho(xd,0)
%       [min;max]
% OUTPUTS:
%   kd,kv,kp := combo of gains with best result
%   beta := best convergence rate

% Default params
gains0 = [1;1;1];
options = optimset();
flagNonAugSys = false; % Flag to solve for the non-augmented system

ivarargin=1;
len_varargin = length(varargin);
while ivarargin<=len_varargin
    switch lower(varargin{ivarargin})
        case 'init_gains'
            % Initial gains must be sent as a 3x1 vector
            ivarargin=ivarargin+1;
            gains0 = varargin{ivarargin};
        case 'showprogress'
            options.PlotFcns = @optimplotfval;
        case 'maxiter'
            % Requires a positive integer
            ivarargin=ivarargin+1;
            options.MaxIter = varargin{ivarargin};
        case 'nonaugmentedsystem'
            % Solve the LMI for the standard non-augmented system
            % The path to the LMI function must be added, 
            % IE run the commend: addpath('../../pointmass_3D/');
            flagNonAugSys = true;
        otherwise
            % ignore the entry and increase index until next entry is
            % a string
            while ivarargin+1 <= len_varargin && ~ischar(varargin{ivarargin+1})
                ivarargin=ivarargin+1;
            end
    end
    ivarargin=ivarargin+1;
end

% Define bisection search as a function. Note: since fminsearch tries to
% minimize, but we want to maximize, the sign needs to be changed
% gains = [kd;kv;kp]
if flagNonAugSys
    bisectionSearch = @(gains) -TR3xR3_betaBisection(e_A,e_B,gains(2),gains(1),0,varargin{:});
    gains0 = gains0(1:2); % Clip the last gain
else
    bisectionSearch = @(gains) -TR3xR3_betaBisection(e_A,e_B,gains(2),gains(1),gains(3),varargin{:});
end

% Solve
gains = fminsearch(bisectionSearch,gains0,options);

% Return results
kd = gains(1); kv = gains(2); 
if flagNonAugSys
    kp = 0;
else
    kp = gains(3);
end
beta = -bisectionSearch(gains);

end

