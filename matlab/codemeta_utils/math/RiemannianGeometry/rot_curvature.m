%Computes the value of R(X,Y)Z, where R is the curvature tensor
%function w=rot_curvature(R,X,Y,Z,W)
%Inputs
%   R       the rotation at which the tangent vectors are based
%   X,Y,Z,W four tangent vectors
%
%If only three tangent vectors are supplied, return R(X,Y)Z
%If the fourth tangent vector (W) is supplied, return <R(X,Y)Z,W>   
%Reference
%   Do Carmo, "Riemannian Geometry", p.103

%Note: this formula is only valid on Lie groups with bi-invariant metric
function w=rot_curvature(R,X,Y,Z,W)
if (isa(X, 'function_handle'))
    X = X(R);
end
if (isa(Y, 'function_handle'))
    Y = Y(R);
end
if (isa(Z, 'function_handle'))
    Z = Z(R);
end
bracket=@(A,B) A*B-B*A;
X=R'*X;
Y=R'*Y;
Z=R'*Z;
w=-R*bracket(bracket(X,Y),Z)/4;
if exist('W','var')
    if (isa(X, 'function_handle'))
        W = W(R);
    end
    w=rot_metric(R,w,W);
end
