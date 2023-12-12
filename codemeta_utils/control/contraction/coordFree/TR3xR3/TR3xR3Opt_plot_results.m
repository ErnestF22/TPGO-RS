function [] = TR3xR3Opt_plot_results(x,t,varargin)
% Last Edited: Oct. 27 2020 by Bee Vang
% Create informative plots for dynamics evolving on TR^3xR^3
% INPUTS:
%   x := state vector 12 x 1 with [xpos;xvel;xref;kd;kv;kref]
%   t := simulation time array
% OUTPUTS:
%   N/A -- plots

%% Set default parameters
flagEnableTitle = false;
flagPlotControl = false;
flagSaveFigures = false;
flagPlotGains = false;
%% Load optional parameters
ivarargin=1;
len_varargin = length(varargin);
while ivarargin<=len_varargin
    switch lower(varargin{ivarargin})
        case 'enabletitle'
            flagEnableTitle = true;
        case 'plot_control'
            % Plot the control vector in R^3. Requires function handle to
            % the control function and state vector x from rotRef_simulation
            flagPlotControl = true;
            ivarargin = ivarargin + 1;
            control_func = varargin{ivarargin};
        case 'plot_gains'
            flagPlotGains = true;
        otherwise
            % ignore the entry and increase index until next entry is
            % a string
            while ivarargin+1 <= len_varargin && ~ischar(varargin{ivarargin+1})
                ivarargin=ivarargin+1;
            end
    end
    ivarargin=ivarargin+1;
end

%% Extract the States
xpos = x(:,1:3); xvel = x(:,4:6); xref = x(:,7:9);
if size(x,2) >= 12
    kd = x(:,10); kv = x(:,11); kref = x(:,12);
end

%% Generate Plots
% Plot distance from xpos to 0
figure
norm_xpos = vecnorm(xpos,2,2);
plot(t,norm_xpos,'LineWidth',3);
set(gca,'TickLabelInterpreter','latex')
ylabel('$\|x_{pos}\|$','Interpreter','latex')
if flagEnableTitle
    title('Distance from x_{pos} to 0')
end

% Plot distance from xpos to xref
figure
norm_xpos_xref = vecnorm(xpos-xref,2,2);
plot(t,norm_xpos_xref,'LineWidth',3);
set(gca,'TickLabelInterpreter','latex')
ylabel('$\|x_{pos}-x_{ref}\|$','Interpreter','latex')
if flagEnableTitle
    title('Distance from x_{pos} to x_{ref}')
end

% Plot speed
figure
norm_velocity = vecnorm(xvel,2,2);
plot(t,norm_velocity,'LineWidth',3);
set(gca,'TickLabelInterpreter','latex')
ylabel('$\|x_{vel}\|$','Interpreter','latex')
if flagEnableTitle
    title('Velocity Error')
end

% Plot the control effort
if flagPlotControl
    figure
    controlVec = zeros(3,length(t));
    for ii = 1:length(t)
        tempControl = control_func(t(ii),x(ii,:)');
        controlVec(:,ii) = tempControl(1:3);
    end
    plot(t,controlVec,'LineWidth',3);
    set(gca,'TickLabelInterpreter','latex')
    ylabel('$\Gamma$','Interpreter','latex')
    if flagEnableTitle
        title('Control Inputs')
    end
    if flagSaveFigures
        customSaveFig("figControl"+addSaveFigName);
    end
end

% Plot the gains
if flagPlotGains
    figure
    plot(t,[kd,kv,kref],'LineWidth',3);
    legend('k_{d}','k_{v}','k_{ref}')
    set(gca,'TickLabelInterpreter','latex')
    ylabel('Gains','Interpreter','latex')
    if flagEnableTitle
        title('Gains Over Time')
    end
end

end

function [] = customSaveFig(saveFilename)
savefig(saveFilename);
stdFigSize = [3 1.25];
% savefigure(char(saveFilename),'epsc',stdFigSize);
pause(1);
end