function [] = rotRefOptAll_plotSim(datafile,varargin)
% LAST EDITED: Nov. 14, 2020 by Bee Vang
% Plot results from simulations for optimal controller with dynamic gains, 
% metric parameters, and convergence rate
% Create informative plots for dynamics evolving on TSO(3)xSO(3)
% INPUTS:
%   x := state vector 24x1 of [R;\omega;Rref;kd;kv;kref]
%   t := simulation time array
% OUTPUTS:
%   N/A -- just plots

%% Default parameters
close all; % Make sure we're starting with clean figures
% Add the file into a cell array
add_dataname{1} = datafile;
add_labelname{1} = 'First';
iAddFilesCounter = 1;
LINEWIDTH = 1;
DISCONT_TOL = 20;
% Set figure formatting
setFigFontSize(8);
flagSaveFigures = false;
flagPlotContractionMat = false;
flagPlotControl = false;
flagEnableTitle = false;
flagPlotGreshBounds = false;
flagPlotGains = false;
flagPlotLyapunov_krasovski = false;
flagPlotConvergenceRate = false;
flagPlotMetric = false;
flagStaticParamsSolveConvRate = false;
%% Load optional parameters
ivarargin=1;
len_varargin = length(varargin);
while ivarargin<=len_varargin
    switch lower(varargin{ivarargin})
        case 'save_figure'
            % Set the initial angular velocity
            flagSaveFigures = true;
            ivarargin = ivarargin + 1;
            prefixSaveFigName = varargin{ivarargin};
        case 'plot_contraction_mat'
            % Plot the eigenvalue of the contraction matrix along the
            % trajectory. Requires function handle to contraction matrix.
            flagPlotContractionMat = true;
        case 'plot_control'
            % Plot the control vector in R^3. Requires function handle to
            % the control function and state vector x from rotRef_simulation
            flagPlotControl = true;
        case 'enable_title'
            flagEnableTitle = true;
            % Requires a pre-title name
            ivarargin = ivarargin + 1;
            preTitle = varargin{ivarargin};
        case 'plot_lyapunov_krasovski'
            flagPlotLyapunov_krasovski = true;
            ivarargin = ivarargin + 1;
            gains = varargin{ivarargin};
            kd = gains(1); kv = gains(2); kref = gains(3);
            ivarargin = ivarargin + 1;
            M_nn = varargin{ivarargin};
        case 'plot_gains'
            % Plot the gains over time (assumes x has atleast 24 rows)
            flagPlotGains = true;
        case 'plot_convergence_rate'
            % Plot the optimal convergence rate from opt problem
            flagPlotConvergenceRate = true;
        case 'plot_metric'
            % Plot the min eigenvalue of the metric matrix M_nn
            flagPlotMetric = true;
        case 'static_parameters_solve_convergence_rate'
            % Skip optimization process, used for simulating unchanged
            % parameter system
            flagStaticParamsSolveConvRate = true;
        case 'addfile'
            % send additional files to compare in plots
            iAddFilesCounter = iAddFilesCounter + 1;
            % requires filename
            ivarargin = ivarargin + 1;
            add_dataname{iAddFilesCounter} = varargin{ivarargin};
            % requires a data label name
            ivarargin = ivarargin + 1;
            add_labelname{iAddFilesCounter} = varargin{ivarargin};
        case 'datalabel'
            % Update the data label for the first data file
            ivarargin = ivarargin + 1;
            add_labelname{1} = varargin{ivarargin};
        otherwise
            % ignore the entry and increase index until next entry is
            % a string
            while ivarargin+1 <= len_varargin && ~ischar(varargin{ivarargin+1})
                ivarargin=ivarargin+1;
            end
    end
    ivarargin=ivarargin+1;
end

%% Extract States
for aa = 1:iAddFilesCounter
    DataStruct = load(add_dataname{aa});
    x = DataStruct.x;
    [R,w,RRef,kd,kv,kref,M_nn,b]=rotDynOptAll_stateUnpack(x);
    t = DataStruct.t;

%% Plot the distance from R to Identity
    figure(1)
    normR_I = zeros(length(t),1);
    for ii = 1:length(t)
        normR_I(ii) = norm(rot_log(R(:,:,ii),eye(3)));
    end
    plot(t,normR_I,'LineWidth',LINEWIDTH);
    hold on
    legend_1{aa} = add_labelname{aa};

%% Plot the distance from R to RRef
    figure(2)
    normR_Ref = zeros(length(t),1);
    for ii = 1:length(t)
        normR_Ref(ii) = norm(rot_log(R(:,:,ii),RRef(:,:,ii)));
    end
    plot(t,normR_Ref,'LineWidth',LINEWIDTH);
    hold on
    legend_2{aa} = add_labelname{aa};
    
%% Plot the distance from RRef to idenity
    figure(3)
    normRRef_I = zeros(length(t),1);
    for ii = 1:length(t)
        normRRef_I(ii) = norm(rot_log(RRef(:,:,ii),eye(3)));
    end
    plot(t,normRRef_I,'LineWidth',LINEWIDTH);
    hold on
    legend_3{aa} = add_labelname{aa};

%% Plot the angular speed
    figure(4)
    normW = vecnorm(w);
    plot(t,normW,'LineWidth',LINEWIDTH);
    hold on
    legend_4{aa} = add_labelname{aa};

%% Plot the control effort
    if flagPlotControl
        control_func = DataStruct.control;
        controlVec = zeros(3,length(t));
        mod_controlVec = controlVec;
        flag_DISCONT_CTRL = false;
        for ii = 1:length(t)
            tempControl = control_func(t(ii),x(:,ii));
            controlVec(:,ii) = tempControl(1:3);
            % For plotting with discontinous controls
            mod_controlVec(:,ii) = tempControl(1:3);
            if ii > 1
                tempDiff = abs(tempControl(1:3)-controlVec(:,ii-1));
                idxDiscont = find(tempDiff>DISCONT_TOL);
                if ~isempty(idxDiscont)
                    mod_controlVec(idxDiscont,ii)=inf;
                    flag_DISCONT_CTRL= true;
                end
            end
        end
        
        % Plot all control of each individual data set as one plot
        pause(1) % Seems that the first the gcf doesnt get updated
        figure
        if flag_DISCONT_CTRL
            plot(t,controlVec,'--','LineWidth',LINEWIDTH);
            hold on
            set(gca,'ColorOrderIndex',1)
            plot(t,mod_controlVec,'LineWidth',LINEWIDTH);
        else
            plot(t,controlVec,'LineWidth',LINEWIDTH);
        end
        set(gca,'TickLabelInterpreter','latex')
        ylabel('$\Gamma$','Interpreter','latex')
        xlabel('$\mathrm{Time } (s)$','Interpreter','latex')
        if flagEnableTitle
            title([preTitle '(' add_labelname{aa} '): Gains over Time'])
        end
        if flagSaveFigures
            customSaveFig(prefixSaveFigName+"_"+add_labelname{aa}+"_Control");
        end
        pause(1)
        close % close this figure after completion incase another plot is using this particular figure number
        
        % Plot to show norm(u) and sum(norm(u))
        figure(5)
        subplot(2,1,1)
        normU = vecnorm(controlVec);
        plot(t,normU,'LineWidth',LINEWIDTH);
        hold on
        legend_5_2{aa} = add_labelname{aa};
        
        subplot(2,1,2)
        t_diff = diff(t); t_diff(end+1) = mean(t_diff);
        sumControlEffort = cumsum(normU'.*t_diff);
    %     sumControlEffort = cumsum(normU);
        plot(t,sumControlEffort,'LineWidth',LINEWIDTH);
        hold on
        legend_5_4{aa} = add_labelname{aa};

        % Plot distance to eye(3) vs cummulative control norm
        figure(6)
        plot(sumControlEffort,normR_I,'LineWidth',LINEWIDTH);
        hold on
        legend_6{aa} = add_labelname{aa};        
    end
    
% %     figure (TODO: complete, will have to resolve optimization problem)
% %     % plot the rate of change of norm(u)^2
% %     controlRate = zeros(length(t),1);
% %     for ii = 1:length(t)
% %         controlRate(ii) = 2*[kd*rot_vee(R,rot_log(R,RRef))'*rot_vee(R,rot_log(R,RRef)) - kv*w'*rot_vee(R,rot_log(R,RRef)),...
% %                     -kd*rot_vee(R,rot_log(R,RRef))'*w + kv*(w'*w),...
% %                     0]*gains_dot;            
% %     end
% %     plot(t,controlRate);
% %     if flagEnableTitle
% %         title([preTitle '(' add_labelname{aa} '): Control Rate']);
% %     end
    
    %% Plot Convergence rate
    if flagPlotConvergenceRate
        figure(7)
        if flagStaticParamsSolveConvRate
            % The system was simulated w/out any parameters, but we want to see
            % what is the optimal convergence rate when solved point-wise
            b = zeros(length(t),1);
            for ii = 1:length(t)
                b(ii) = rotRefOptAll_optProblem(x(:,ii),'opt_convergence_rate_only');
            end
        end
        plot(t,b,'LineWidth',LINEWIDTH);
        hold on
        legend_7{aa} = add_labelname{aa};
    end

%% Plot max eigenvalues of contraction matrix
% This plot should be done after convergence rate in case we want to check
% for updated convergence rates for the static parameter case
    if flagPlotContractionMat
        figure(8)
        maxEval = zeros(length(t),1);
        % Create function handle to contraction matrix on TSO(3)xSO(3) as
        % function of (R,w,RRef,kd,kv,kref,M_nn,b)
        M_contract_params = @(b,kd,kv,kref,M_nn) rotRef_contractionMat2(b,kd,kv,kref,M_nn,'sym');
        for ii = 1:length(t)
            if isfield(DataStruct,'GainScheduledParams') && norm(rot_log(RRef(:,:,ii),eye(3))) <= 1e-5
                % If Gain Schedule Controller and RRef = eye(3), then local
                % controller is being used and the contraction matrix on TSO(3)
                % should be used
                temp_M_nn = M_nn(:,:,ii);
                M_contraction = rotBundle_contractionMat(@(R) R*hat3(w(:,ii)),b(ii),kd(ii),kv(ii),[temp_M_nn(1,1),temp_M_nn(1,2),temp_M_nn(2,2)],'sym');
                M_contract = M_contraction(R(:,:,ii));
            else
                % Contraction matrix as function of (R,w,RRef)
                M_contract_state = M_contract_params(b(ii),kd(ii),kv(ii),kref(ii),M_nn(:,:,ii));
                % Numerical contraction matrix
                M_contract = M_contract_state(R(:,:,ii),w(:,ii),RRef(:,:,ii));
            end
            maxEval(ii) = max(eig(M_contract));
        end
        plot(t,maxEval,'LineWidth',LINEWIDTH);
        hold on
        legend_8{aa} = add_labelname{aa};
    end

%% Plot metric eigenvalue
    if flagPlotMetric
        % These can be plotted separately
        figure
        metric_Eig = zeros(length(t),1);
        for ii = 1:length(t)
            metric_Eig(ii) = min(eig(M_nn(:,:,ii)));
        end
        subplot(1,2,1)
        plot(t,metric_Eig,'LineWidth',LINEWIDTH);
        set(gca,'TickLabelInterpreter','latex')
        ylabel('Min Eignvalue','Interpreter','latex')
        if flagEnableTitle
            title([preTitle ': Min Eigenvalue of Metric'])
        end
        subplot(1,2,2)
        m1_all = M_nn(1,1,:); m1_all = m1_all(:);
        m2_all = M_nn(1,2,:); m2_all = m2_all(:);
        m3_all = M_nn(2,2,:); m3_all = m3_all(:);
        m4_all = M_nn(3,3,:); m4_all = m4_all(:);
        m5_all = M_nn(2,3,:); m5_all = m5_all(:);
        m6_all = M_nn(1,3,:); m6_all = m6_all(:);
        plot(t,[m1_all, m2_all, m3_all, m4_all, m5_all, m6_all],'LineWidth',LINEWIDTH);
        set(gca,'TickLabelInterpreter','latex')
        ylabel('Metric Values','Interpreter','latex')
        legend('m_1','m_2','m_3','m_4','m_5','m_6');
        if flagEnableTitle
            title([preTitle '(' add_labelname{aa} '): Metric Parameters'])
        end
        if flagSaveFigures
            customSaveFig(prefixSaveFigName+"_"+add_labelname{aa}+"_Metric");
        end
    end
%% Plot the gains
    if flagPlotGains
        % These can be plotted separately
        figure
        plot(t,[kd;kv;kref],'LineWidth',LINEWIDTH);
        set(gca,'TickLabelInterpreter','latex')
        ylabel('Gains','Interpreter','latex')
        legend('k_d','k_v','k_{ref}','FontSize',14);
        if flagEnableTitle
            title([preTitle '(' add_labelname{aa} '): Gains over Time'])
        end
        if flagSaveFigures
            customSaveFig(prefixSaveFigName+"_"+add_labelname{aa}+"_Gains");
        end
    end
end

%% Add legends and labels
figure(1)
set(gca,'TickLabelInterpreter','latex')
ylabel('$\|(\log_{I_{3}}R)^{\vee}\|$','Interpreter','latex')
xlabel('$\mathrm{Time } (s)$','Interpreter','latex')
legend(legend_1);

figure(2)
set(gca,'TickLabelInterpreter','latex')
ylabel('$\|(\log_{R}R_{ref})^{\vee}\|$','Interpreter','latex')
xlabel('$\mathrm{Time } (s)$','Interpreter','latex')
legend(legend_2);

figure(3);
set(gca,'TickLabelInterpreter','latex')
ylabel('$\|(\log_{I_{3}}R_{ref})^{\vee}\|$','Interpreter','latex')
xlabel('$\mathrm{Time } (s)$','Interpreter','latex')
legend(legend_3);

figure(4);
set(gca,'TickLabelInterpreter','latex')
ylabel('$\|\omega\|$','Interpreter','latex')
xlabel('$\mathrm{Time } (s)$','Interpreter','latex')
legend(legend_4);

figure(5);
% subplot(2,2,[1 3]);
% set(gca,'TickLabelInterpreter','latex');
% ylabel('$\Gamma$','Interpreter','latex');
% xlabel('$\mathrm{Time } (s)$','Interpreter','latex')
% legend([legend_5_1{:}]); %[cell{:}] puts all elements into one array
subplot(2,1,1);
set(gca,'TickLabelInterpreter','latex');
ylabel('$\|\Gamma\|$','Interpreter','latex')
xlabel('$\mathrm{Time } (s)$','Interpreter','latex')
legend(legend_5_2);
subplot(2,1,2);
set(gca,'TickLabelInterpreter','latex');
ylabel('$\sum\|\Gamma\|$','Interpreter','latex')
xlabel('$\mathrm{Time } (s)$','Interpreter','latex')
legend(legend_5_4);

figure(6);
set(gca,'TickLabelInterpreter','latex')
xlabel('$\sum \|u\|$','Interpreter','latex')
ylabel('$\|(\log_{I_{3}}R)^{\vee}\|$','Interpreter','latex')
legend(legend_6);

figure(7);
set(gca,'TickLabelInterpreter','latex')
ylabel('Convergence Rate','Interpreter','latex')
xlabel('$\mathrm{Time } (s)$','Interpreter','latex')
legend(legend_7);

figure(8);
set(gca,'TickLabelInterpreter','latex')
ylabel('$\lambda_{max}$','Interpreter','latex')
xlabel('$\mathrm{Time } (s)$','Interpreter','latex')
legend(legend_8,'Location','Best');

%% Add titles and save (optional)
if flagEnableTitle
    figure(1); title([preTitle ': Distance from R to I']);
    figure(2); title([preTitle ': Distance from R to RRef']);
    figure(3); title([preTitle ': Distance from RRef to I']);
    figure(4); title([preTitle ': Angular Speed']);
    figure(5);
%     subplot(2,2,[1 3]); title([preTitle ': Control Inputs']);
    subplot(2,1,1); title([preTitle ': Control Norm']);
    subplot(2,1,2); title([preTitle ': Cumulative Control Norms']);
    figure(6); title([preTitle ': Distance vs Norm(u)']);
    figure(7); title([preTitle ': Convergence Rate over Time']);
    figure(8); title([preTitle ': Max Contraction Eigenvalue']);
end  
if flagSaveFigures
    figure(1); customSaveFig(prefixSaveFigName+"_NormR_I");
    figure(2); customSaveFig(prefixSaveFigName+"_NormR_RRef");
    figure(3); customSaveFig(prefixSaveFigName+"_NormRRef_I");
    figure(4); customSaveFig(prefixSaveFigName+"_NormW");
    figure(5); customSaveFig(prefixSaveFigName+"_Control");
    figure(6); customSaveFig(prefixSaveFigName+"_Distance_vs_normU");
    figure(7); customSaveFig(prefixSaveFigName+"_ConvergRate");
    figure(8); customSaveFig(prefixSaveFigName+"_Eig");
end

end

function [] = customSaveFig(saveFilename)
savefig(saveFilename);
pause(1);
stdFigSize = [3 1.25];
savefigure(char(saveFilename),'epsc',stdFigSize);
pause(1);
end