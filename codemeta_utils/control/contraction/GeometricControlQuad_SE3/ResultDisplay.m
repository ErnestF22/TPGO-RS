function ResultDisplay(timeSteps,stateVariables,desiredStateVariables,positionError,velocityError)

    idx=round(linspace(1,length(timeSteps),round(length(timeSteps))));
    %size(idx)
%%    
    %%% plot trajectory %%%
    figure;
    %subplot(2,2,1);
    plot3(stateVariables(idx,1),stateVariables(idx,2),stateVariables(idx,3),'-r','LineWidth',2);
    hold on;
    plot3(desiredStateVariables(idx,1),desiredStateVariables(idx,2),desiredStateVariables(idx,3),':k','LineWidth',2);
    hold off;
    grid on;
    title('Comparision between Specified Trajectory and Tracking Result');
    legend('trajectory','desired trajectory');
    %axis equal;
    xlabel('x-axis');
    ylabel('y-axis');
    zlabel('z-axis');
    
    figure
    subplot(3,1,1);
    plot(timeSteps(idx),stateVariables(idx,1),'-r','LineWidth',2);%-b
    hold on;
    plot(timeSteps(idx),desiredStateVariables(idx,1),':k','LineWidth',2);%:r
    hold off;
    grid on;
    title('X-Component');
    legend('x position','desired x position');
    %axis equal;
    xlabel('time, sec');
    ylabel('x');
    
    subplot(3,1,2);
    plot(timeSteps(idx),stateVariables(idx,2),'-r','LineWidth',2);
    hold on;
    plot(timeSteps(idx),desiredStateVariables(idx,2),':k','LineWidth',2);%:m
    hold off;
    grid on;
    title('Y-Component');
    legend('y position','desired y position');
    %axis equal;
    xlabel('time, sec');
    ylabel('y');
    
    subplot(3,1,3);
    plot(timeSteps(idx),stateVariables(idx,3),'-r','LineWidth',2);
    hold on;
    plot(timeSteps(idx),desiredStateVariables(idx,3),':k','LineWidth',2);
    hold off;
    grid on;
    title('Z-Component');
    legend('z position','desired z position');
    %axis equal;
    xlabel('time, sec');
    ylabel('z');
%%    
    %%% plot angular velocity %%%
    figure
    subplot(3,1,1);
    plot(timeSteps(idx),stateVariables(idx,16),'-r','LineWidth',2);%-b
    hold on;
    plot(timeSteps(idx),desiredStateVariables(idx,16),':k','LineWidth',2);%:r
    hold off;
    grid on;
    title('X-Component of {\Omega}');
    legend('{\Omega}_x','desired {\Omega}_x');
    %axis equal;
    xlabel('time, sec');
    ylabel('{\Omega}_x');
    
    subplot(3,1,2);
    plot(timeSteps(idx),stateVariables(idx,17),'-r','LineWidth',2);
    hold on;
    plot(timeSteps(idx),desiredStateVariables(idx,17),':k','LineWidth',2);%:m
    hold off;
    grid on;
    title('Y-Component of {\Omega}');
    legend('{\Omega}_y','desired {\Omega}_y');
    %axis equal;
    xlabel('time, sec');
    ylabel('{\Omega}_y');
    
    subplot(3,1,3);
    plot(timeSteps(idx),stateVariables(idx,18),'-r','LineWidth',2);
    hold on;
    plot(timeSteps(idx),desiredStateVariables(idx,18),':k','LineWidth',2);
    hold off;
    grid on;
    title('Z-Component of {\Omega}');
    legend('{\Omega}_z','desired {\Omega}_z');
    %axis equal;
    xlabel('time, sec');
    ylabel('{\Omega}_z');
%%    
    %%% Plot error %%%
    figure;
    subplot(2,1,1);
    plot(timeSteps(idx),positionError(idx),'-r','LineWidth',2);
    grid on;
    title('error in position');
    legend('position error');
    xlabel('Time');
    ylabel('position error');
    subplot(2,1,2);
    plot(timeSteps(idx),velocityError(idx),'-r','LineWidth',2);
    grid on;
    title('error in velocity');
    legend('velocity error');
    xlabel('time, sec');
    ylabel('velocity error');
    
%%
    %%% Plot orientation
    orientation=zeros(numel(idx),3);
    desiredOrientation=zeros(numel(idx),3);
    
    for iStep=idx
        tempRotationMatrix=reshape(stateVariables(iStep,7:15),3,3);
        tempEulerAngles=rotationMatrix2EulerAngles(tempRotationMatrix);
        orientation(iStep,:)=tempEulerAngles(:,1)';    
        tempRotationMatrix=reshape(desiredStateVariables(iStep,7:15),3,3);
        tempEulerAngles=rotationMatrix2EulerAngles(tempRotationMatrix);
        desiredOrientation(iStep,:)=tempEulerAngles(:,1)';
    end
    
    figure
    subplot(3,1,1);
    plot(timeSteps(idx),orientation(idx,1),'-r','LineWidth',2);%-b
    hold on;
    plot(timeSteps(idx),desiredOrientation(idx,1),':k','LineWidth',2);%:r
    hold off;
    grid on;
    title('rolling angle');
    legend('{\phi}','desired {\phi}');
    %axis equal;
    xlabel('time, sec');
    ylabel('{\phi}');
    
    subplot(3,1,2);
    plot(timeSteps(idx),orientation(idx,2),'-r','LineWidth',2);
    hold on;
    plot(timeSteps(idx),desiredOrientation(idx,2),':k','LineWidth',2);%:m
    hold off;
    grid on;
    title('pitch angle');
    legend('{\theta}','desired {\theta}');
    %axis equal;
    xlabel('time, sec');
    ylabel('{\theta}');
    
    subplot(3,1,3);
    plot(timeSteps(idx),orientation(idx,3),'-r','LineWidth',2);
    hold on;
    plot(timeSteps(idx),desiredOrientation(idx,3),':k','LineWidth',2);
    hold off;
    grid on;
    title('yaw angle');
    legend('{\psi}','desired {\psi}');
    %axis equal;
    xlabel('time, sec');
    ylabel('{\psi}');   
    
end