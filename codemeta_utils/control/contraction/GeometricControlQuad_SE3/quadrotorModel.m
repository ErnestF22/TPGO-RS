function [deltaStateVariables,desiredStateVariables,force,moment]=quadrotorModel(timeSteps,currentStateVariables,quadrotor)

    %%% Get Desired State Variables %%%
    if(quadrotor.trajectory.flag)
        
        [desiredTrajectory]=generateTrajectory(timeSteps);
        desiredPosition=desiredTrajectory.position;
        desiredLinearVelocity=desiredTrajectory.velocity;
        desiredAcceleration=desiredTrajectory.acceleration;
        bodyDir1=desiredTrajectory.bodyDir1;
        desiredAngularVelocity=quadrotor.trajectory.desiredAngularVelocity;
        desiredDomegadt=zeros(3,1);
        
    else
        
        desiredPosition=quadrotor.trajectory.desiredPosition;
        desiredLinearVelocity=[0;0;0];
        desiredAcceleration=[0;0;0];
        bodyDir1=quadrotor.trajectory.basis1Direction;
        %desiredRotation=eye(3);
        desiredAngularVelocity=quadrotor.trajectory.desiredAngularVelocity;
        desiredDomegadt=zeros(3,1);

    end
    
    mass=quadrotor.params.mass;
    inertiaTensor=quadrotor.params.J;
    g=quadrotor.params.g;
    %e1=quadrotor.params.e1;
    %e2=quadrotor.params.e2;
    e3=quadrotor.params.e3;

    currentPosition=currentStateVariables(1:3);
    currentLinearVelocity=currentStateVariables(4:6);
    currentRotationMatrix=reshape(currentStateVariables(7:15),3,3);
    currentAngularVelocity=currentStateVariables(16:18);
    bodyDir3=currentRotationMatrix(:,3);
   
%%
    %%% Position Control %%%
    
    positionError=currentPosition-desiredPosition;
    velocityError=currentLinearVelocity-desiredLinearVelocity;
    
    % gains for position error and velocity error
    Kp=quadrotor.params.Kp;
    Kv=quadrotor.params.Kv;
    temp=(-Kp*positionError-Kv*velocityError+mass*desiredAcceleration-mass*g*e3);
    desiredBasis3Direction=-temp/norm(temp);
    force=dot(-temp,bodyDir3);
    
%%
    %%% Orientation Control %%%

    bodyDir2=cross(desiredBasis3Direction,bodyDir1);
    bodyDir2=bodyDir2/norm(bodyDir2);
    %projection_b1d=-cross(desiredBasis3Direction,basis2Direction);
    desiredRotation=[cross(bodyDir2,desiredBasis3Direction) bodyDir2 desiredBasis3Direction];    
    
    %gains for orientation error and angular velocity error
    Kr = quadrotor.params.Kr;
    Kw = quadrotor.params.Kw;
    
    %compute error in angular velocity and orientation
    orientationError=vee(desiredRotation'*currentRotationMatrix-currentRotationMatrix'*desiredRotation)/2;
    
    angularVelocityError=currentAngularVelocity-currentRotationMatrix'*desiredRotation*desiredAngularVelocity;
    
    moment=-Kr*orientationError-Kw*angularVelocityError+cross(currentAngularVelocity,inertiaTensor*currentAngularVelocity)...
        -inertiaTensor*(hat(currentAngularVelocity)*currentRotationMatrix'*desiredRotation*desiredAngularVelocity...
        -currentRotationMatrix'*desiredRotation*desiredDomegadt);

    %%% flight dynamics model %%%
    dXdt=currentLinearVelocity;                                     %x_dot=v
    dVdt=g*e3-(force/mass)*currentRotationMatrix*e3;                %m*v_dot=m*g*e3-f*R*e3
    dRdt=currentRotationMatrix*hat(currentAngularVelocity);         %R_dot=R*omega_hat
    dOmegadt=inertiaTensor\(moment-cross(currentAngularVelocity,inertiaTensor*currentAngularVelocity));%J*omega_dot+omega x (J*omega)=M
    
    desiredStateVariables = [desiredPosition; desiredLinearVelocity;%compute desired state variables
        reshape(desiredRotation,9,1); desiredAngularVelocity];
    deltaStateVariables=[dXdt;dVdt;reshape(dRdt,9,1);dOmegadt;];    %compute state variables' changes in a time step dt


end