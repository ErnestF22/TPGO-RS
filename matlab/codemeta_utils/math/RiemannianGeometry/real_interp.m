%Interpolation of values in R^d
%function real_interp(t,x,tQuery,varargin)
%This function is essentially a wrapper for interp1. It interprets the
%values in x as a sequence, where the last dimension correspond to the time
%index (that is, the length of t and the length along the last dimension of
%x must coincide).
function xQuery=real_interp(t,x,tQuery,varargin)
sz=size(x);
if sz(2)==1
    sz(2)=[];
end
NSamples=sz(end);
if length(t)~=NSamples
    error('Number of time istants must be equal to the number of samples (length along the last dimension of x)');
end
x=reshape(x,[],NSamples);
xQuery=interp1(t,x',tQuery,varargin{:})';
if length(sz)>2
    xQuery=reshape(xQuery,[sz(1:end-1) length(tQuery)]);
end


