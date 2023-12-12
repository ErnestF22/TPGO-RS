%Pack the state of a rigid body transformation and its velocities in a vector
%function x=rigidDyn_statePack(R,T,w,v)
%The vectors x are divided as follows:
%   x(1:9)      Rotation R ([3 x 3] matrix)
%   x(10:12)    Angular velocity w ([3 x 1] vector)
%   x(13:15)    Translation T ([3 x 1] vector)
%   x(16:18)    Linear velocity v ([3 x 1] vector)
%Note: the vector is ordered in this way to be cross-compatible with the
%rotDyn_* functions.
function x=rigidDyn_statePackRT(R,w,T,v)
NPoses=size(R,3);
x=zeros(18,NPoses);
x(1:12,:)=rotDyn_statePack(R,w);
x(13:18,:)=[T;v];
