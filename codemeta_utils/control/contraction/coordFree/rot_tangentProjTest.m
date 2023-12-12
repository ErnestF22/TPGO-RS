function [r, varargout ] = rot_tangentProjTest( R, X, Y, t )
% Check if the covariant derivative of Y wrt X creates another vector field
% in the tangent space

if ~exist('t','var') || isempty(t)
    t=[];
end

D_X_Y = rot_covar(X, Y);
[r,t]=funEval(@(t) norm(D_X_Y(R(t)) - rot_tangentProj(R(t), ...
    D_X_Y(R(t))), 'fro'),t);

if nargout > 1
  varargout{1} = t;
end

end

