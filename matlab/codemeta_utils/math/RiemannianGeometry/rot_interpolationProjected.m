%Rotation interpolation by using interpolation in the ambient space and projections
%function rot_interpolationProjected(t,R,tQuery,varargin)
function RQuery=rot_interpolationProjected(t,R,tQuery,varargin)
QQuery=real_interp(t,reshape(R,9,[]),tQuery,varargin{:});
RQuery=rot_proj(reshape(QQuery,3,3,[]));
