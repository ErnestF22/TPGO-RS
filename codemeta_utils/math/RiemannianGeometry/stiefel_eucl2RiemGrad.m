%Projects H onto the tangent space of the Stiefel manifold at Y
function H=stiefel_eucl2RiemGrad(Y,H)
H=H-Y*(H'*Y);