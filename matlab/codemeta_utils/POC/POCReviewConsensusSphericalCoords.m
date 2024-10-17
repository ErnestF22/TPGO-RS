function POCReviewConsensusSphericalCoords
syms psi_i psi_j phi_i phi_j
sphereCoords=@(psi,phi) [cos(psi)*cos(phi); sin(psi)*cos(phi); sin(phi)];
DsphereCoords=@(psi,phi) [-sin(psi)*cos(phi) -cos(psi)*sin(phi); cos(psi)*cos(phi) -sin(psi)*sin(phi); 0 cos(phi)];
hat=@(v) [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];

gamma_i=sphereCoords(psi_i,phi_i);
gamma_j=sphereCoords(psi_j,phi_j);


% %test DsphereCoords
% [psi,~,~,dpsi]=real_randGeodFun(0);
% [phi,~,~,dphi]=real_randGeodFun(0,'speed',rand);
% gamma=@(t) sphereCoords(psi(t),phi(t));
% dgamma=@(t) DsphereCoords(psi(t),phi(t))*[dpsi;dphi];
% funCheckDer(gamma,dgamma)

dgamma_i=hat(gamma_i)^2*gamma_j;
dpsi_i=(-sin(psi_j-psi_i)*cos(phi_j))/cos(phi_i);
dphi_i=sin(phi_i)*cos(phi_j)*cos(psi_i-psi_j)-cos(phi_i)*sin(phi_j);

dgamma_iB=DsphereCoords(psi_i,phi_i)*[dpsi_i; dphi_i];

simplify(dgamma_i-dgamma_iB)
