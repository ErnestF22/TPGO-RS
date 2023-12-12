randn('state',0)
R0=eye(3);
v1=rot_randTangentNormVector(R0);
U1=rot_exp(R0,0.1*v1);
U2=rot_exp(R0,-0.1*v1);

fun=@(R) rot_dist(R(:,:,1),R(:,:,2))^2;
der=@(R) cat(3,-rot_log(R(:,:,1),R(:,:,2)),-rot_log(R(:,:,2),R(:,:,1)));

Ropt=lie_minimize(rot_funs(), fun, der, cat(3,U1,U2),'showcost');
disp(Ropt)

% A=randn(3);
% fun=@(R) trace(A'*(R1-R2));
% der=@(R) 4*cat(3,rot_tangentProj(R(:,:,1),A),-rot_tangentProj(R(:,:,2),A));
% 
% Ropt=lie_minimize(rot_funs(), fun, der, cat(3,U1,U2),'showcost');
% disp(Ropt)

