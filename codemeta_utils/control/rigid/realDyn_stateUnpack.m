%Unpack state for translation dynamics
%function [T,v]=realDyn_stateUnpack(x)
%The vector x is divided as follows:
%   x(end-5:end-3)  Translation T
%   x(end-2:end)    Translation velocity v
%Note: the definition of the indeces starting from the end is used to make
%this function cross-compatible with the rigidDyn_* functions.
function [T,v]=realDyn_stateUnpack(x)
T=x(end-5:end-3,:);
v=x(end-2:end,:);
