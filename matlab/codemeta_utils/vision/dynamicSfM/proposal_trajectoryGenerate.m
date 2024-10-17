function proposal_trajectoryGenerate
T0=[0;-1;1];
T1=[1;0;0];
T2=[-1;1;0];
T3=T2;
R0=eye(3);
R1=rot([0;0;1]*20/180*pi);
R2=rot(-[0;0;1]*20/180*pi);
GReference_R=cat(3,R0,R1,R2,R0);
GReference_T=[T0 T1 T2 T3]/2;
GReference_t=[0 5 10 15]/2;
nameSequence='proposal_trajectory';

quadrotor_controlPDRT_generateDataset(nameSequence,GReference_R,GReference_T,GReference_t)
dynamicSfM_trajectoryInertial('proposal_trajectory',false)
dynamicSfM_trajectoryInertialVisual('proposal_trajectoryInertial')
dynamicSfm_testNoise
