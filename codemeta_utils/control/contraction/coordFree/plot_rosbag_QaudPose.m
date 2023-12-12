function [] = plot_rosbag_QaudPose(rosbagFile,dataStruct,varargin)
% Load rosbagFile, extract, and plot the '/mavros/local_position/pose' and 
% '/mavros/setpoint_pose/local' topics from a quadrotor flight
% INPUTS:
%   dataStruct := A structure containing the following informations
%       .quadName := quadrotor name as a string ""
%       .kd := the geometric controller gain for position error
%       .kv := the geometric controller gain for velocity error
%       .MC_ROLLRATE_MAX := the max roll rate in degree/s from px4 params
%       .MC_PITCHRATE_MAX := the max pitch rate in degree/s from px4 params
%       .MC_YAWRATE_MAX := the max yaw rate in degree/s from px4 params

%% Set Parameters
plotTime_start = 0.0; % Set the time to start the plot
plotTime_end = inf; % Set the time to end the plot (inf == length of rosbag)

%% Load optional parameters
ivarargin=1;
len_varargin = length(varargin);
while ivarargin<=len_varargin
    switch lower(varargin{ivarargin})
        case 'setstarttime'
            % set the ploting time to start this many seconds relative to
            % the rosbag start time
            ivarargin=ivarargin+1;
            plotTime_start = varargin{ivarargin};
        case 'setendtime'
            % set the ploting time to end this many seconds relative to
            % the rosbag start time, should be greated than plotTime_start
            ivarargin=ivarargin+1;
            plotTime_end = varargin{ivarargin};  
        otherwise
            % ignore the entry and increase index until next entry is
            % a string
            while ivarargin+1 <= len_varargin && ~ischar(varargin{ivarargin+1})
                ivarargin=ivarargin+1;
            end
    end
    ivarargin=ivarargin+1;
end

if exist(rosbagFile, 'file')
    bag=rosbag(rosbagFile);
else
    fprintf('File does not exist\n');
    return;
end

if ~exist('dataStruct','var')
    % If no data information, assume no gains
    dataStruct.quadName = "";
    dataStruct.kd = 0.0;
    dataStruct.kv = 0.0;
    dataStruct.MC_ROLLRATE_MAX = 0.0;
    dataStruct.MC_PITCHRATE_MAX = 0.0;
    dataStruct.MC_YAWRATE_MAX = 0.0;
else
    % Fill in missing data
    if ~isfield(dataStruct,'quadName')
        dataStruct.quadName = "";
    end
    if ~isfield(dataStruct,'kd')
        dataStruct.kd = 0.0;
    end
    if ~isfield(dataStruct,'kv')
        dataStruct.kv = 0.0;
    end
    if ~isfield(dataStruct,'MC_ROLLRATE_MAX')
        dataStruct.MC_ROLLRATE_MAX = 0.0;
    end
    if ~isfield(dataStruct,'MC_PITCHRATE_MAX')
        dataStruct.MC_PITCHRATE_MAX = 0.0;
    end
    if ~isfield(dataStruct,'MC_YAWRATE_MAX')
        dataStruct.MC_YAWRATE_MAX = 0.0;
    end
end

%% Extract Data (message lengths may vary for each topic)
% Extract the state (Pose, Orientation) as matlab timeseries data
localPoseSel = select(bag,'Topic',dataStruct.quadName + '/mavros/local_position/pose');
localPoseData_ts = timeseries(localPoseSel,...
    'Pose.Position.X', 'Pose.Position.Y', 'Pose.Position.Z', ...
    'Pose.Orientation.W', 'Pose.Orientation.X', 'Pose.Orientation.Y', 'Pose.Orientation.Z'); % [w,x,y,z] to match quat2rotm function

% Extract the desired state (Pose, Orientation) as matlab timeseries data
localSetPoseSel = select(bag,'Topic',dataStruct.quadName + '/mavros/setpoint_position/local');
localSetPoseData_ts = timeseries(localSetPoseSel,...
    'Pose.Position.X', 'Pose.Position.Y', 'Pose.Position.Z', ...
    'Pose.Orientation.W', 'Pose.Orientation.X', 'Pose.Orientation.Y', 'Pose.Orientation.Z'); % [w,x,y,z] to match quat2rotm function

% Extract the velocity data
localVelSel = select(bag,'Topic',dataStruct.quadName + '/mavros/local_position/velocity_local');
% localVelSel = select(bag,'Topic',dataStruct.quadName + '/mavros/local_position/velocity');
localVelData_ts = timeseries(localVelSel,...
    'Twist.Linear.X', 'Twist.Linear.Y', 'Twist.Linear.Z', ...
    'Twist.Angular.X', 'Twist.Angular.Y', 'Twist.Angular.Z');

% Adjust time to start at 0
poseTime = localPoseData_ts.Time;
setPoseTime = localSetPoseData_ts.Time;
velTime = localVelData_ts.Time;
% Find which is greater and adjust time accordingly
tempTimeHolder = [poseTime(1), setPoseTime(1), velTime(1)];
minTime = min(tempTimeHolder);
% Adjust time to start at 0
poseTime = poseTime - minTime;
setPoseTime = setPoseTime - minTime;
velTime = velTime - minTime;
% set the plotting index
poseTime_start_idx = find(poseTime >= plotTime_start, 1, 'first');
poseTime_end_idx = find(poseTime <= plotTime_end, 1, 'last');
setPoseTime_start_idx = find(setPoseTime >= plotTime_start, 1, 'first');
setPoseTime_end_idx = find(setPoseTime <= plotTime_end, 1, 'last');
velTime_start_idx = find(velTime >= plotTime_start, 1, 'first');
velTime_end_idx = find(velTime <= plotTime_end, 1, 'last');

% Store pose data, each pose is stored along a row
currentPosition = localPoseData_ts.Data(:,1:3);
desiredPosition = localSetPoseData_ts.Data(:,1:3);
currentAngVel = localVelData_ts.Data(:,4:6);
currentQuat = localPoseData_ts.Data(:,4:end);
desiredQuat = localSetPoseData_ts.Data(:,4:end);
currentR = quat2rotm(currentQuat);
desiredR = quat2rotm(desiredQuat);

%% Plots
% Plot the position
figure
ax = axes;
ax.ColorOrder = [1 0 0; 0 1 0; 0 0 1];
hold on
plot(poseTime(poseTime_start_idx:poseTime_end_idx), currentPosition(poseTime_start_idx:poseTime_end_idx,:), '-', 'LineWidth', 3)
plot(setPoseTime(setPoseTime_start_idx:setPoseTime_end_idx), desiredPosition(setPoseTime_start_idx:setPoseTime_end_idx,:), '--', 'LineWidth', 3)
legend('x','y','z','setX','setY','setZ');
title('Position: '+rosbagFile, 'Interpreter', 'none')
xlabel('Time')
ylabel('Position Error','Interpreter','latex')

% Plot the orientation in quaternions (note q=-q, so plot may be useless),
% will plot as norm(log(R,Rd))
figure
ax = axes;
ax.ColorOrder = [0 0 0; 1 0 0; 0 1 0; 0 0 1];
hold on
plot(poseTime(poseTime_start_idx:poseTime_end_idx), currentQuat(poseTime_start_idx:poseTime_end_idx,:), '-', 'LineWidth', 3)
% Note taking negative of the quaternion may not always give the "correct" graph
plot(setPoseTime(setPoseTime_start_idx:setPoseTime_end_idx) , -desiredQuat(setPoseTime_start_idx:setPoseTime_end_idx,:), '--', 'LineWidth', 3) 
legend('w','x','y','z','setW','setX','setY','setZ');
title('Orientation: '+rosbagFile, 'Interpreter', 'none')
xlabel('Time')
ylabel('Attitude Error (Quat.)','Interpreter','latex')

% Compute position and rotation distance errors, and desired control
normR = zeros(size(currentR,3),1);
normPosition = zeros(size(currentR,3),1);
controlVec = zeros(3,size(currentR,3));
for ii = 1:size(currentR,3)
    % Assume that the desiredR is the most recent setpoint before
    % poseTime(ii)
    desiredR_idx = find(setPoseTime<=poseTime(ii),1,'last');
    if isempty(desiredR_idx)
        % if poseTime is before any setPoseTime, use the first entry
        desiredR_idx = 1;
    end
    % If desired is NaN, set to Identity
    if any(any(isnan(desiredR(:,:,desiredR_idx))))
        desiredR(:,:,desiredR_idx) = eye(3);
    end
    normR(ii) = norm(rot_log(currentR(:,:,ii),desiredR(:,:,desiredR_idx)));   
    normPosition(ii) = norm( (currentPosition(ii,:)-desiredPosition(ii,:))' );
    
    % Assume that the current velocity is closest to the current R
    [~,currentVel_idx] = min(abs(poseTime(ii)-velTime));
    uR = dataStruct.kd*rot_vee(currentR(:,:,ii),rot_log(currentR(:,:,ii),desiredR(:,:,desiredR_idx)));
    uW = -dataStruct.kv*currentAngVel(currentVel_idx,:)';
    controlVec(:,ii) = uR+uW;
end

% Plot position distance error on R^3
figure
plot(poseTime(poseTime_start_idx:poseTime_end_idx), normPosition(poseTime_start_idx:poseTime_end_idx), 'LineWidth', 3);
title('Position Error: '+rosbagFile, 'Interpreter', 'none')
xlabel('Time')
ylabel('$\|$Position Error$\|$','Interpreter','latex')

% Plot rotation distance error on SO(3)
figure
plot(poseTime(poseTime_start_idx:poseTime_end_idx), normR(poseTime_start_idx:poseTime_end_idx), 'LineWidth', 3)
title('Rotation Error: '+rosbagFile, 'Interpreter', 'none')
xlabel('Time')
ylabel('$\|$Attitude Error$\|$','Interpreter','latex')

% Plot control vector along with limits on each axis
figure
h=plot(poseTime(poseTime_start_idx:poseTime_end_idx), controlVec(:,poseTime_start_idx:poseTime_end_idx), 'LineWidth', 3);
c = get(h,'Color'); % get the color
hold on
% plot the limits
yline(dataStruct.MC_ROLLRATE_MAX*pi/180.0, '--', 'Color', c{1}, 'LineWidth', 3);
yline(dataStruct.MC_PITCHRATE_MAX*pi/180.0, '--', 'Color', c{2}, 'LineWidth', 3);
yline(dataStruct.MC_YAWRATE_MAX*pi/180.0, '--', 'Color', c{3}, 'LineWidth', 3);
yline(-dataStruct.MC_ROLLRATE_MAX*pi/180.0, '--', 'Color', c{1}, 'LineWidth', 3);
yline(-dataStruct.MC_PITCHRATE_MAX*pi/180.0, '--', 'Color', c{2}, 'LineWidth', 3);
yline(-dataStruct.MC_YAWRATE_MAX*pi/180.0, '--', 'Color', c{3}, 'LineWidth', 3);
title('Control Input: ' +rosbagFile, 'Interpreter', 'none')
xlabel('Time');
ylabel('Control Input','Interpreter','latex')
legend('Roll', 'Pitch', 'Yaw', 'Roll Limit', 'Pitch Limit', 'Yaw Limit');
end

