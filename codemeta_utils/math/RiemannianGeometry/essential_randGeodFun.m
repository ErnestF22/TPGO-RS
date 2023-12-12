function [Qt,vt,Q0,v0,v0Vec]=essential_randGeodFun(Q0,varargin)

if ~exist('Q0','var') || isempty(Q0)
    Q0=essential_randn();
end

v0=essential_randTangentNormVector(Q0);
[Qt,vt]=essential_geodFun(Q0,v0,varargin{:});
if nargout>4
    v0Vec=essential_vee(Q0,v0);
end