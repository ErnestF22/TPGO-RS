function POCProjectiveRecovery
%These tests eventually ended up in euclideanUpgradeWithNormConstraints.m
resetRands()

S=randn(3,4);
A=squeeze(sphere_randn([1;0;0],[],30));
NA=size(A,2);

SA=S*homogeneous(A,4);

invS=inv([S; 0 0 0 1]);
invS=invS(1:3,:);

% %result should be all ones
% disp(cnorm(invS*homogeneous(SA,4)))

Q=invS'*invS;
% Qext=[invS eye(3)]'*[invS eye(3)];

B=[];
for iA=1:NA
    B=[B;vec(homogeneous(SA(:,iA),4)*homogeneous(SA(:,iA),4)')'];
end

% disp((B*vec(Qext(1:4,1:4)))')

% B=[];
% for iA=1:NA
%     SAi=SA(:,iA);
%     Bi=[vec(SAi*SAi'); 2*SAi; 1];
%     B=[B;Bi'];
% end
% 
% QB=[vec(Qext(1:3,1:3));Qext(5:7,4); Q(4,4)];
% disp((B*QB)')

B=[];
for iA=1:NA
    a=SA(:,iA);
    Bi=[a(1)^2, 2*a(1)*a(2), 2*a(1)*a(3), 2*a(1)*1, a(2)^2, 2*a(2)*a(3), 2*a(2)*1, a(3)^2, 2*a(3)*1];
    B=[B;Bi];
end
b=1;
u=b*ones(NA,1);

% tic
% qSym=sym('q',[10 1]);
% QSym=[qSym(1) qSym(2) qSym(3) qSym(4); qSym(2) qSym(5) qSym(6) qSym(7); qSym(3) qSym(6) qSym(8) qSym(9); qSym(4) qSym(7) qSym(9) qSym(10)];
% %Split the matrix into two components
% QSym1=QSym;
% QSym1(end,end)=0;
% QSym2=sym(zeros(size(Q)));
% QSym2(end,end)=1;
% a=sym('a');
% Qa=a*QSym1+(1-a)*QSym2;
% %[pol,powers]=coeffs((prod(eig(Qa))),a);
% %pol=prod(eig(Qa));
% pol=det(Qa);
% toc

qEst=[B\u;b];
QEst1=[qEst(1) qEst(2) qEst(3) qEst(4); qEst(2) qEst(5) qEst(6) qEst(7); qEst(3) qEst(6) qEst(8) qEst(9); qEst(4) qEst(7) qEst(9) 0];
QEst2=diag([0 0 0 qEst(10)]);
p=determinantLinearCombinationCoeff(QEst1,b);
a=-p(2)/p(1);
QEst=a*QEst1+(1-a)*QEst2;
det(QEst)

[U,S,V]=svd(QEst);
invSEst=diag(sqrt(diag(S(1:3,1:3))))*V(:,1:3)';

AEst=invSEst*homogeneous(SA,4);

cnorm(AEst)


% pol=subs(pol,qSym,qEst);
% toc
% [polCoeffs,powers]=coeffs(pol);
% polCoeffs=double(polCoeffs);
% toc
% polEst=double(subs(pol,qEst));
% roots(polCoeffs)

%determinantLinearCombination_test


keyboard
% %CVX solution
% cvx_begin
% variable QEst(4,4) symmetric
% r=B*vec(QEst)-ones(NA,1);
% minimize(r'*r);
% QEst == hermitian_semidefinite(4)
% cvx_end
% 
% disp(QEst)
% svd(QEst)'
% rank(Q)
% rank(QEst)
% [Q./QEst]

function determinantLinearCombinationCoeff_testRoots()
A=randn(5);
E=diag([0 0 0 0 1]);
b=rand;
p=determinantLinearCombinationCoeff(A,b);
x=-p(2)/p(1);
disp([det(a*A+(1-a)*E*b) determinantLinearCombination(a,A,b)])

function p=determinantLinearCombinationCoeff(A,b)
% compute the coefficients of the polynomial in x given by the 
% determinat of x*A+(1-x)*E*b, where E=diag([0....0 1]), from the
% determinant of A, x and b are scalars
N=size(A,1);
c1=det(A);
c2=b*det(A(1:N-1,1:N-1));
p1=c1-c2;
p2=c2;
%d=p1*a^N+p2*a^(N-1);
p=[p1 p2 zeros(1,N-1)];

function determinantLinearCombination_test
A=randn(5);
E=diag([0 0 0 0 1]);
b=rand;
a=rand;
disp([det(a*A+(1-a)*E*b) determinantLinearCombination(a,A,b)])

function d=determinantLinearCombination(a,A,b)
% compute the determinat of a*A+(1-a)*E*b, where E=diag([0....0 1]), from the
% determinant of A, a and b are scalars

% N=size(A,1);
% %d=a^N*det(A)+(1-a)*a^(N-1)*b*det(A(1:N-1,1:N-1));
% c1=det(A);
% c2=b*det(A(1:N-1,1:N-1));
% p1=c1-c2;
% p2=c2;
% %d=p1*a^N+p2*a^(N-1);
% p=[p1 p2 zeros(1,N-1)];
p=determinantLinearCombinationCoeff(A,b);
d=polyval(p,a);

function determinantLinearCombination2_test
A=randn(5);
b=rand;
E=diag([0,0,0,0,1]);
disp([det(A+b*E) determinantLinearCombination2(A,b)])

function d=determinantLinearCombination2(A,b)
% compute the determinat of A+E*b, where E=diag([0....0 1]), from the
% determinant of A
N=size(A,1);
d=det(A)+b*det(A(1:N-1,1:N-1));

function buildDeterminantPolynomial()
qSym=sym('q',[10 1]);
QSym=[qSym(1) qSym(2) qSym(3) qSym(4); qSym(2) qSym(5) qSym(6) qSym(7); qSym(3) qSym(6) qSym(8) qSym(9); qSym(4) qSym(7) qSym(9) qSym(10)];
%Split the matrix into two components
QSym1=QSym;
QSym1(end,end)=0;
QSym2=sym(zeros(size(Q)));
QSym2(end,end)=1;
a=sym('a');
