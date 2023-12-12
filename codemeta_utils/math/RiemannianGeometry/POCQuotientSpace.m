function POCQuotientSpace
R0=eye(5);
[Q,dQ]=rot_randGeodFun(eye(3));

R=@(t) R0*blkdiag(eye(2),Q(t));
dR=@(t) R0*blkdiag(zeros(2),dQ(t));

%funCheckDer(R,dR)

%disp(rot_vee(R(0),dR(0)))

bVecVert=[eye(3);zeros(7,3)];
bVecHoriz=[zeros(3,7); eye(7)];

bVert=rot_hat(R0,bVecVert);
bHoriz=rot_hat(R0,bVecHoriz);


X=rot_hat(R0,bVecHoriz*randn(7,1));
Y=rot_hat(R0,bVecHoriz*randn(7,1));
V=rot_hat(R0,bVecVert*randn(3,1));
W=rot_hat(R0,bVecVert*randn(3,1));
T=rot_randTangentNormVector(R0);

%rot_exp(R0,V)

brXY=rot_bracket(R0,X,Y);
brXV=rot_bracket(R0,V,X);
brVW=rot_bracket(R0,V,W);
brVT=rot_bracket(R0,V,T);

%disp('<[X,V],W>=')
%disp(rot_metric(R0,brXV,W))
disp('X [X,V]')
disp([rot_vee(R0,X) rot_vee(R0,brXV)])

disp('V [V,W]')
disp([rot_vee(R0,V) rot_vee(R0,brVW)])

% disp('[V,T]')
% disp(rot_vee(R0,brVT))

% disp('[X,Y]')
% disp(rot_vee(R0,brXY))

%disp(rot_vee(R0,reshape(rot_bracket(R0,bHoriz,bHoriz),5,5,[])))

