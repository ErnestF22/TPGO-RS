%%% GEOMETRIC TRACKING CONTROLLER TEST %%%
%%% Literature references:
%%% (1) Geometric Tracking Control of a Quadrotor UAV on SE(3);
%%% (2) Control of Complex Maneuvers for a Quadrotor UAV using Geometric
%%% Methods on SE(3);

clear;
clc;
close all; 

m=[4.34 0.5];
J=[0.082 0.0845 0.1377;0.00557 0.00557 0.0105];

kp=[16 16 16;10 10 9];
kv=[5.6 5.6 5.6;5 5 10];
kr=[8.81 200];
kw=[2.54 1];

paramPlan=1;

%%
%%% System Parameters Initialization %%%

quadrotor.params.mass=m(paramPlan);                 %the total mass
quadrotor.params.J=diag(J(paramPlan));              %inertia tensor in kg(m^2)
quadrotor.params.g=9.81;                            %gravitational acceleration

quadrotor.params.e1=[1;0;0] ;                       %interial frame basis
quadrotor.params.e2=[0;1;0] ;                       %interial frame basis
quadrotor.params.e3=[0;0;1] ;                       %interial frame basis

quadrotor.params.d=0.315;                           %distance of center of mass from fram center in m        
quadrotor.params.c=8.004*10e-4;                     %fixed constant in m

quadrotor.params.Kp=diag(kp(paramPlan+1,:));        %gains
quadrotor.params.Kv=diag(kv(paramPlan+1,:));
quadrotor.params.Kr=kr(paramPlan+1);
quadrotor.params.Kw=kw(paramPlan+1);  

%parameter of trajectory
quadrotor.trajectory.flag=true;                     %flag indicates whether or not the function generateTrajectory() is to be called
quadrotor.trajectory.desiredPosition=nan;           %desired position
quadrotor.trajectory.desiredAngularVelocity=[0;0;0];%desired angular velocity

quadrotor.trajectory.velocity=nan;
quadrotor.trajectory.acceleration=nan;
quadrotor.trajectory.bodyDir1=nan;   

%%% State Variable Initialization %%%
initialPosition=[0;0;0];                            %intial spatial position   
initialLinearVelocity=zeros(3,1);                   %intial linear velocity
initialRotationMatrix=eye(3);                       %initial rotation matrix from the body-fixed frame to the inertial frame
initialAngularVelocity=zeros(3,1);                  %intial angular velocity in the body-fixed frame

%concatenate the entire initial condition into a 18-by-1 vector
initialStateVariables=[initialPosition;initialLinearVelocity;
    reshape(initialRotationMatrix,9,1);initialAngularVelocity];

%%%% Simulation Starts %%%
odeopts=odeset('RelTol',1e-8,'AbsTol',1e-9);       %set up relative error tolerance and absolute error tolerance
[timeSteps,stateVariables]=ode15s(@quadrotorModel,[0 20],initialStateVariables,odeopts,quadrotor) ;

idx=round(linspace(1,length(timeSteps),round(1*length(timeSteps))));

desiredStateVariables=zeros(size(stateVariables));
positionError=zeros(numel(idx));
velocityError=zeros(numel(idx));
Force=zeros(numel(idx));
Moment=zeros(numel(idx),3);

for iStep=idx
   [~,desiredVariableTranspose,tempForce,tempMoment]=quadrotorModel(timeSteps(iStep),stateVariables(iStep,:)',quadrotor);
   desiredStateVariables(iStep,:) = desiredVariableTranspose';
   Force(iStep,1)=tempForce;
   Moment(iStep,:)=tempMoment';
   positionError(iStep)=norm(stateVariables(iStep,1:3)-desiredStateVariables(iStep,1:3));%compute position error
   velocityError(iStep)=norm(stateVariables(iStep,4:6)-desiredStateVariables(iStep,4:6));%compute velocity error
end

ResultDisplay(timeSteps,stateVariables,desiredStateVariables,positionError,velocityError);%plot results








