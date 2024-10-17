function r = CurvatureR(u,v,w)

% Riemannian curvature for sphere

r = (w'*u)*v-(w'*v)*u;


end

