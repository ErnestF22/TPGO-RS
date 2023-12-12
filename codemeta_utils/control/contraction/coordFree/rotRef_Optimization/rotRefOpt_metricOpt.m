function [betaOut] = rotRefOpt_metricOpt(x,varargin)
% Last Edited: Nov. 9 2020 by Bee Vang
% Solve for a new metric given the current state and gains. Goal is to find
% a metric that gives the best convergence rate bound \beta by bisection
% search. Calls the function rotRef_betaBisection, for setup see rotRef_contractionTest.m
% INPUTS:
%   x := [24x1] current state x (packed). The first 12 elements are on
%       TSO(3), last 9 are on SO(3), and remaining 3 are the gains
%       [kd;kv;kref]
% OUTPUTS:
%   betaOut := positive scalar mini guaranteed exp. conv. rate
%   M_nnOut := [3x3] nonnatural metric that satisfies contraction condidtion


%% Setup parameters (mag_W, mag_R, mag_RRef)
[R,w,RReference] = rotDyn_stateUnpack(x,'augmentedsystem');
kd = x(22);
kv = x(23);
kref = x(24);

mag_W = norm(w);
mag_R = norm(rot_log(R,RReference));
mag_RRef = norm(rot_log(RReference,eye(3)));

%% Find best beta and return (function handles must be included for the LMI opt problem)
betaOut = rotRef_betaBisection(mag_R,mag_RRef,kd,kv,kref,mag_W,varargin{:});

end

