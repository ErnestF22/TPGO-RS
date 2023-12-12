% POC_Compare_Convergence_Part2.m
% This script runs the simulations and plot the control effort vs distance
% error
% NOTE: THIS SCRIPT SHOULD BE RUN AFTER PART1 W/OUT CLEARING THE WORKSPACE

%% Init
c = parcluster(); % set our cluster to local

%% Collect the optimal gain results
gainResults = cell(2*length(maxGains),5);
for ii = 1:2*length(maxGains)
    tempJob = findJob(c,'ID',cell2mat(jobsID(ii,1)));
    tempResults = fetchOutputs(tempJob);
    tempResultsMat = cell2mat(tempResults(1:4));
    gainResults(ii,1:4) = tempResults; % Store gains as a row vector
    gainResults(ii,5) = jobsID(ii,2); % Store the system type
    % Delete the jobs and clear resources
    delete(tempJob)
    clear tempJob
end
delete c; clear c;

%% Simluate
figure;
% Plot fake points for legend
plot(0,0,'rx');
hold on;
plot(0,0,'go');
optsOde=odeset('MaxStep',0.01);
control_Scale = 1/1e4;
for ii = 1:size(gainResults,1)
    kd = gainResults{ii,1}; kv = gainResults{ii,2}; kp = gainResults{ii,3};
    beta = gainResults{ii,4};

    if strcmpi(gainResults{ii,5},'augmentedsystem')
        control=@(t,x) TR3xR3_rigidControlPD(x, v1, mass, kd, kv, gradf);
        closedLoop=@(t,x) TR3xR3_rigidModel(x, x1, mass, control(t,x), gradf, kp);
        for jj = 1:sim_perGain
            x0 = randn(3,1); x0 = x0_dist*x0/norm(x0); % Random initial condition
            [t,x]=ode45(closedLoop,[0 TFinal],[x0;v0;.75*x0], optsOde);
            % Sum the control effort
            u_sum = 0;
            for kk = 1:length(t)
                u_sum = u_sum + norm(control(t(kk),x(kk,:)'),2);
            end
            % Plot norm(x(TFinal),2)^2 vs. u_sum
            semilogx(norm(x(end,:),2),u_sum*control_Scale,'rx');
            hold on;
        end
    else
        control=@(t,x) contraction3D_rigidControlPD(x, x1, v1, mass, kd, kv, gradf);
        closedLoop=@(t,x) contraction3D_rigidModel(x, mass, control(t,x));
        for jj = 1:sim_perGain
            x0 = randn(3,1); x0 = x0_dist*x0/norm(x0); % Random initial condition
            [t,x]=ode45(closedLoop,[0 TFinal],[x0;v0], optsOde);
            % Sum the control effort
            u_sum = 0;
            for kk = 1:length(t)
                u_sum = u_sum + norm(control(t(kk),x(kk,:)'),2);
            end
            % Plot norm(x(TFinal),2)^2 vs. u_sum
            semilogx(norm(x(end,:),2),u_sum*control_Scale,'go');
            hold on;
        end
    end
end

legend('Augmented','Non-augmented')
xlabel(['norm(x('  num2str(TFinal) '),2)'])
ylabel(['sum(norm(u,2))*' num2str(control_Scale)])

