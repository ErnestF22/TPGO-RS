function [] = POC_EvalTermsContraction(kd,kv_list,maxD,magW,varargin)
% Compute the terms that are multipled by m2 and m3 assuming beta = 0. Find
% which inequality breaks for a given kd. Each element of kv_list will
% generate a plot that assumes kv is fixed.
% INPUTS:
%   kd := scalar gain
%   kv_list := list of scalar gains to test
%   maxD := max distance error
%   magW := max speed
% OUTPUTS:
%   results := a function of (kd,kv) [8x3] coeffs., column 1 is multiply by m2, colm. 2 is mult.
%       by m3, the sum of the terms should be less than or equal to colm. 3

% Define Parameters
syms m2 m3;
Lambda = [1;maxD/2*cot(maxD/2)];
minM2 = 0.013;
maxM2 = 0.017;
minM3 = 0;
maxM3 = 0.005;
maxSteps = 200;
beta = 0;
lineStyle={'r-','r:','g-','g:','m-','m:','k-','k:'};

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'minm2'
            ivarargin=ivarargin+1;
            minM2=varargin{ivarargin};
        case 'maxm2'
            ivarargin=ivarargin+1;
            maxM2=varargin{ivarargin};
        case 'minm3'
            ivarargin=ivarargin+1;
            minM3=varargin{ivarargin};
        case 'maxm3'    
            ivarargin=ivarargin+1;
            maxM3=varargin{ivarargin};
        case 'maxsteps'
            ivarargin=ivarargin+1;
            maxSteps=varargin{ivarargin};
        case 'beta'
            ivarargin=ivarargin+1;
            beta=varargin{ivarargin};
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

% Define Coefficients for the constraints [m2*coeff1 + m3*coeff2 + const]

coeff = @(kv,beta) [ ...
    [-kd*Lambda+beta-kv/2-magW/4,...
    -kd*Lambda/2+magW^2/8+magW/4*kv,...
    beta+1/2*[1;1] ];...
    % Gerg Disc 3 and 4
    [-kd*Lambda-beta+kv/2-magW/4,...
    kd*Lambda/2+magW^2/8+magW/4*kv,...
    beta-1/2*[1;1] ];...
    % Gerg Disc 5 and 6
    [1+beta-kv/2-magW/4*[1;1],...
    beta-kv-kd*Lambda/2+magW^2/8+magW/4*kv,...
    +1/2*[1;1] ];...
    % Gerg Disc 7 and 8
    [1-beta+kv/2-magW/4*[1;1],...
    beta-kv+kd*Lambda/2+magW^2/8+magW/4*kv,...
    -1/2*[1;1] ]];

mTermsMat = [m2;m3;1];

% For each kv in kv_list, plot the feasible region (assuming fixed kv)
figure;
iRows = ceil(sqrt(length(kv_list))); % and number of columns
for i = 1:length(kv_list)
    kv = kv_list(i);
    constants = coeff(kv,beta);
    constraints = constants*mTermsMat<=0;
    constraints(end+1,:) = [m3>m2^2];
%     constraints = [m3>m2^2];
    
    % Plot the feasible regions
    ax(i) = subplot(iRows,iRows,i);
    funImageConstraints(constraints,m2,m3,{},'xgrid',linspace(minM2,maxM2,maxSteps),'ygrid',linspace(minM3,maxM3,maxSteps));
    hold on
    for j = 1:size(constants,1)
        plotline(linspace(0,1),linspace(0,1),constants(j,1:2)',constants(j,3),char(lineStyle(j)));
    end
    
    title(['kd = ' num2str(kd) ', kv = ' num2str(kv) ', beta = ' num2str(beta)]);
    xlabel('m2');
    ylabel('m3');
    xlim([minM2 maxM2]);
    ylim([minM3 maxM3]);
end
% linkaxes(ax,'xy');

