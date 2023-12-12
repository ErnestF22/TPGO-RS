function [kd,kv,beta] = rotBundle_gainSearch(mag_R, mag_W, varargin)
% Use derivative free method Nelder and Mead Algorithm and find optimal
% gains kd and kv. This method is a builtin as the fminsearch matlab
% function.
% INPUTS:
%   mag_R := positive scalar representing max distance from R to eye(3)
%   magW := positive scalar representing maximum velocity magnitude
% OUTPUTS:
%   kd,kv := combo of gains with best result
%   beta := best convergence rate

%% Default params
warning('off','all')
gains0 = [1;1];
options = optimset();
saveDirectory = 'rotBundle_contractionResults';

%% Optional params
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
        case 'gainsearch_maxiter'
            % Requires a positive integer
            ivarargin=ivarargin+1;
            options.MaxIter = varargin{ivarargin};
        case 'savedir'
            ivarargin=ivarargin+1;
            saveDirectory = varargin{ivarargin};
        otherwise
            % ignore the entry and increase index until next entry is
            % a string
            while ivarargin+1 <= len_varargin && ~ischar(varargin{ivarargin+1})
                ivarargin=ivarargin+1;
            end
    end
    ivarargin=ivarargin+1;
end

%% Nelder Meld
% Define bisection search as a function. Note: since fminsearch tries to
% minimize, but we want to maximize, the sign needs to be changed
% gains = [kd;kv;kp]
bisectionSearch = @(gains) -rotBundle_betaBisection_new(mag_R,gains(2),gains(1),mag_W,varargin{:});

% Solve
[gains,beta] = fminsearch(bisectionSearch,gains0,options);
% gains = fminunc(bisectionSearch,gains0,options);

% Return results
kd = gains(1); kv = gains(2);
beta = -beta; % Invert sign of bisectionSearch

if ~exist(saveDirectory,'dir')   
    mkdir(saveDirectory);
end

% save to file
clear bisectionSearch
fileSaveName = [saveDirectory '/Results_' datestr(datetime,'yyyymmdd_HHMMss') '.mat'];
while exist(fileSaveName,'file')
    % Wait 1 second before generating a filename
    pause(1);
    fileSaveName = [saveDirectory '/Results_' datestr(datetime,'yyyymmdd_HHMMss') '.mat'];
end
% Save once the filename does not exist
save(fileSaveName);
end