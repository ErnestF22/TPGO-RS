%Computes the value sectional curvature for two orthonormals vectors
%function W=rot_sectionalCurvature(R,X,Y)
%Inputs
%   R       the rotation at which the tangent vectors are based
%   X,Y     two ORTHONORMAL tangent vectors
%
%Reference
%   Do Carmo, "Riemannian Geometry", p.103

%Note: this formula is only valid on Lie groups with bi-invariant metric
function K=rot_sectionalCurvature(R,X,Y)
bracket=@(A,B) A*B-B*A;
X=R'*X;
Y=R'*Y;
Z=bracket(X,Y);
K=rot_metric(eye(3),Z,Z)/4;
