function rconn = riemannian_connection_stiefel(x,u,problem)
%RIEMANNIAN_CONNECTION_STIEFEL
% implementing \nabla_u V = Proj_x D \overline{V}(x)[u]
% V = function i.e., f
% D = eucl. grad.
% x \in Stiefel
% u \in Tangent space to stiefel at x
xStack=matStack(x);
g=matUnstack((problem.L+problem.L')*xStack+problem.P,problem.sz(1));


rconn = stp_manopt(x, g);

end

