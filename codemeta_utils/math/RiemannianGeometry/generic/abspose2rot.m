%function R=abspose2rot(yaw,pitch,roll)
%Pass from Euler angles to rotation matrix

%%AUTORIGHTS%%

function R=abspose2rot(yaw,pitch,roll)
if(exist('yaw','var')==0)
    yaw=0;
end
if(exist('pitch','var')==0)
    pitch=0;
end
if(exist('roll','var')==0)
    roll=0;
end

R=rot(yaw*[0;0;1])*rot(pitch*[1;0;0])*rot(roll*[0;1;0]);

