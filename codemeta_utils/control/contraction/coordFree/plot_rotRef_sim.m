function [] = plot_rotRef_sim(R,w,RRef,t,varargin)
% Create informative plots for dynamics evolving on TSO(3)xSO(3)
% INPUTS:
%   R := [3x3] Matrix array of the current state 
%   w := [3x1] Vector array of the current angular velocity
%   RRef := [3x3] Matrix array of the reference trajectory
%   t := simulation time array
% OUTPUTS:
%   N/A -- just plots
% Set figure formatting
setFigFontSize(8);
axisLimit = [-0.5 1];
stdFigSize = [3 1.25];
% Default parameters
flagSaveFigures = false;
flagPlotContractionMat = false;
flagPlotControl = false;
flagEnableTitle = false;
flagPlotGreshBounds = false;
flagPlotLyapunov = false;
flagPlotLyapunov_krasovski = false;
% Load optional parameters
ivarargin=1;
len_varargin = length(varargin);
while ivarargin<=len_varargin
    switch lower(varargin{ivarargin})
        case 'savefigure'
            % Set the initial angular velocity
            flagSaveFigures = true;
        case 'plot_contraction_mat'
            % Plot the eigenvalue of the contraction matrix along the
            % trajectory. Requires function handle to contraction matrix.
            flagPlotContractionMat = true;
            ivarargin = ivarargin + 1;
            M_contract = varargin{ivarargin};
        case 'plot_control'
            % Plot the control vector in R^3. Requires function handle to
            % the control function and state vector x from rotRef_simulation
            flagPlotControl = true;
            ivarargin = ivarargin + 1;
            control_func = varargin{ivarargin};
            ivarargin = ivarargin + 1;
            x = varargin{ivarargin};
        case 'enabletitle'
            flagEnableTitle = true;
        case 'plot_lyapunov'
            flagPlotLyapunov = true;
            ivarargin = ivarargin + 1;
            gains = varargin{ivarargin};
            kd = gains(1); kv = gains(2); kref = gains(3);
        case 'plot_lyapunov_krasovski'
            flagPlotLyapunov_krasovski = true;
            ivarargin = ivarargin + 1;
            gains = varargin{ivarargin};
            kd = gains(1); kv = gains(2); kref = gains(3);
            ivarargin = ivarargin + 1;
            M_nn = varargin{ivarargin};
        case 'plot_greshbounds'
            % Plot the maximum bounds given by the greshgorin disc along
            % the trajectory
            % Requires: [kd,kv,kref]; beta; M_nonnatural
            flagPlotGreshBounds = true;
            ivarargin = ivarargin + 1;
            gains = varargin{ivarargin};
            kd = gains(1); kv = gains(2); kref = gains(3);
            ivarargin = ivarargin + 1;
            beta = varargin{ivarargin};
            ivarargin = ivarargin + 1;
            M_nn = varargin{ivarargin};
        otherwise
            % ignore the entry and increase index until next entry is
            % a string
            while ivarargin+1 <= len_varargin && ~ischar(varargin{ivarargin+1})
                ivarargin=ivarargin+1;
            end
    end
    ivarargin=ivarargin+1;
end

% Plot the distance from R to identity
figure
normR_I = zeros(length(t),1);
for ii = 1:length(t)
    normR_I(ii) = norm(rot_log(R(:,:,ii),eye(3)));
end
plot(t,normR_I,'LineWidth',3);
set(gca,'TickLabelInterpreter','latex')
ylabel('$\|(\log_{I_{3}}R)^{\vee}\|$','Interpreter','latex')
if flagEnableTitle
    title('Distance from R to I')
end
if flagSaveFigures
    savefigure('figNormR_I','epsc',stdFigSize);
    pause(1);
end

% Plot the distance from R to RRef
figure
normR_Ref = zeros(length(t),1);
for ii = 1:length(t)
    normR_Ref(ii) = norm(rot_log(R(:,:,ii),RRef(:,:,ii)));
end
plot(t,normR_Ref,'LineWidth',3);
set(gca,'TickLabelInterpreter','latex')
ylabel('$\|(\log_{R}R_{ref})^{\vee}\|$','Interpreter','latex')
if flagEnableTitle
    title('Distance from R to RRef')
end
if flagSaveFigures
    savefigure('figNormR_RRef','epsc',stdFigSize);
    pause(1);
end

% Plot the distance from RRef to idenity
figure
normRRef_I = zeros(length(t),1);
for ii = 1:length(t)
    normRRef_I(ii) = norm(rot_log(RRef(:,:,ii),eye(3)));
end
plot(t,normRRef_I,'LineWidth',3);
set(gca,'TickLabelInterpreter','latex')
ylabel('$\|(\log_{I_{3}}R_{ref})^{\vee}\|$','Interpreter','latex')
if flagEnableTitle
    title('Distance from RRef to I')
end
if flagSaveFigures
    savefigure('figNormRRef_I','epsc',stdFigSize);
    pause(1);
end

% Plot the angular speed
figure
normW = vecnorm(w);
plot(t,normW,'LineWidth',3);
set(gca,'TickLabelInterpreter','latex')
ylabel('$\|\omega\|$','Interpreter','latex')
if flagEnableTitle
    title('Angular Speed')
end
if flagSaveFigures
    savefigure('figNormW','epsc',stdFigSize);
    pause(1);
end

% Plot max eigenvalues of contraction matrix
if flagPlotContractionMat
    figure
    maxEval = zeros(length(t),1);
    for ii = 1:length(t)
        M = M_contract(R(:,:,ii),w(:,ii),RRef(:,:,ii));
        M=(M+M')/2;
        maxEval(ii) = max(eig(M));
    end
    plot(t,maxEval,'LineWidth',3);
    set(gca,'TickLabelInterpreter','latex')
    ylabel('$\lambda_{max}(\mathcal{M})$','Interpreter','latex')
    if flagEnableTitle
        title('Max Contraction Eigenvalue')
    end
    if flagSaveFigures
        savefigure('figEig','epsc',stdFigSize);
        pause(1);
    end
end

% Plot the control effort
if flagPlotControl
    figure
    controlVec = zeros(3,length(t));
    for ii = 1:length(t)
        controlVec(:,ii) = control_func(t(ii),x(:,ii));
    end
    plot(t,controlVec,'LineWidth',3);
    set(gca,'TickLabelInterpreter','latex')
    ylabel('$\Gamma$','Interpreter','latex')
    if flagEnableTitle
        title('Control Inputs')
    end
    if flagSaveFigures
        savefigure('figControl','epsc',stdFigSize);
        pause(1);
    end
end

% Plot gershgorin bounds
if flagPlotGreshBounds
    figure
    greshBounds = zeros(3,length(t));
    for ii = 1:length(t)
        greshBounds(:,ii) = rotRef_contractionMatrix_greshBound(normR_Ref(ii),normRRef_I(ii),kd,kv,kref,normW(ii),beta,M_nn);
    end
    plot(t,greshBounds,'LineWidth',3);
    set(gca,'TickLabelInterpreter','latex')
    ylabel('$Greshgorin Bounds$','Interpreter','latex')
    if flagEnableTitle
        title('Gersgorin Bounds')
    end
    if flagSaveFigures
        savefigure('figGreshBounds','epsc',stdFigSize);
        pause(1);
    end
end

% Plot Lyapunov Function
% V(R,w,RRef) = 1/2*norm(w)^2 + kd*1/2*norm(rot_log(R,RRef))^2 + 1/2*norm(rot_log(RRef,eye(3))
if flagPlotLyapunov
    figure
    LyapunovF = zeros(1,length(t));
    for ii = 1:length(t)
        LyapunovF(:,ii) = 1/(2*kd)*normW(ii)^2 + 1/2*normR_Ref(ii)^2 + 1/2*normRRef_I(ii)^2;
    end
    plot(t,LyapunovF,'LineWidth',3);
    hold on
    plot(t,LyapunovF(1)*exp(-.4022*t),'LineWidth',3);
    set(gca,'TickLabelInterpreter','latex')
    ylabel('$Lyapunov$','Interpreter','latex')
    if flagEnableTitle
        title('Lyapunov Function')
    end
    if flagSaveFigures
        savefigure('figLyapunov','epsc',stdFigSize);
        pause(1);
    end
    
    % Plot the derivatives
    figure
    LyapunovF_der = zeros(3,length(t));
    for ii = 1:length(t)
        LyapunovF_der(:,ii) = [-kv/kd*normW(ii)^2;...
            1/2*kref*trace(-rot_log(RRef(:,:,ii),R(:,:,ii))'*rot_log(RRef(:,:,ii),eye(3)));...
            -1/2*kref*trace(rot_log(RRef(:,:,ii),eye(3))'*rot_log(RRef(:,:,ii),eye(3)))];
    end
    LyapunovF_der(4,:) = sum(LyapunovF_der,1);
    plot(t,LyapunovF_der,'LineWidth',3);
    legend('omega','R','Ref','Total');
    set(gca,'TickLabelInterpreter','latex')
    ylabel('$Lyapunov Der$','Interpreter','latex')
    if flagEnableTitle
        title('Lyapunov Function Der')
    end
    if flagSaveFigures
        savefigure('figLyapunov_der','epsc',stdFigSize);
        pause(1);
    end
end

% Plot Krasovski Lyapunov function
% V(R,w,RRef) = rotRef_metric_nonNatural2(Z,X,X,M_nn)^2
if flagPlotLyapunov_krasovski
    figure
    LyapunovF = zeros(1,length(t));
    X = @(R,U,RRef) [U;kd*rot_log(R,RRef)-kv*U;kref*rot_log(RRef,eye(3))];
    for ii = 1:length(t)
        R0=R(:,:,ii);
        U0=R0*hat3(w(:,ii));
        RRef0=RRef(:,:,ii);
        LyapunovF(ii) = rotRef_metric_nonNatural2([R0;U0;RRef0],...
            X(R0,U0,RRef0),X(R0,U0,RRef0),M_nn);
    end
    plot(t,LyapunovF,'LineWidth',3);
    hold on
    plot(t,LyapunovF(1)*exp(-.4022*t),'LineWidth',3);
    set(gca,'TickLabelInterpreter','latex')
    ylabel('$Lyapunov (Krasovski)$','Interpreter','latex')
    if flagEnableTitle
        title('Lyapunov (Krasovski) Function')
    end
    if flagSaveFigures
        savefigure('figLyapunov_krasovski','epsc',stdFigSize);
        pause(1);
    end
end
    
end