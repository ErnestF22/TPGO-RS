function POCProjectiveRecovery()
Ti=randn(3,1);
Tj=randn(3,2);
Ri=rot_randn(eye(3));
Tij=Tj-Ti;
[tij,lij]=cnormalize(Ri'*Tij);
Lij=diag(lij);

disp('Check consistency between all generated data (28)')
disp(norm(Ri*tij*Lij-Tij,'fro'))

% Go to higher dimension
dimStair=5;
Rtilde=rot_randn(eye(dimStair));
%Rtilde=eye(4);
Tij_tilde_pre=Rtilde*[Tij;zeros(dimStair-3,2)];
Ri_tilde_pre=Rtilde*[Ri;zeros(dimStair-3,3)];

disp('Align to get to the case where the last row is zero')
Qalign=align3d(Tij_tilde_pre);
Tij_tilde=Qalign*Tij_tilde_pre;
Ri_tilde=Qalign*Ri_tilde_pre;


disp('Check consistency between all generated data with high dimension (28)')
disp(norm(Ri_tilde*tij*Lij-Tij_tilde,'fro'))

disp('Test computation of Qx')
disp(['It should align the last' num2str(dimStair-2) ' rows of Tijtilde to zero'])
Qx=align2d_EoF(Tij_tilde);

disp(Qx*Tij_tilde)


% Generate a random Rb and the associated Qb
Rb=rot_randn(eye(dimStair-2));
Qb=blkdiag(eye(2),Rb);

disp('Check that Qx''*Qb*Qx leaves Tijtilde invariant')
disp(norm(Qx'*Qb*Qx*Tij_tilde-Tij_tilde))
disp(norm(Ri_tilde*tij*Lij-Qx'*Qb*Qx*Tij_tilde))

disp('And that the new Ritilde generated from the ambiguity')
disp('still gives the same measurements.')
Ri_tilde2=Qx'*Qb*Qx*Ri_tilde;
disp(norm(Ri_tilde*tij*Lij-Ri_tilde2*tij*Lij))

disp('But Ritilde2 does not have the last row to zero')
disp(Ri_tilde2)

% test for creating Jacobian to test for ambiguities
s=randn(3,1);
mu=randn(2,1);
b=randn(dimStair-3,3);
Ri_tilde2_ort=orthComplement(Ri_tilde2);

tangentRi=Ri_tilde2*hat3(s)+Ri_tilde2_ort*b;

disp('Test that tangentRi is in the tangent space at Ri_tilde2')
disp(norm(tangentRi'*Ri_tilde2-tangentRi'*Ri_tilde2,'fro'))

%equation (38) on Overleaf
tangentAmbiguity=vec((Ri_tilde2*hat3(s)+Ri_tilde2_ort*b)*tij*Lij+Ri_tilde2*tij*diag(mu));
tangentAmbiguity2= ...
    vec(Ri_tilde2*hat3(s)*tij*Lij+Ri_tilde2_ort*b*tij*Lij+Ri_tilde2*tij*diag(mu));
M1=zeros(dimStair*2,3);
I3=eye(3);
for iM1=1:3
    M1(:,iM1)=vec(Ri_tilde2*hat3(I3(:,iM1))*tij*Lij);
end
M2=kron((tij*Lij)',Ri_tilde2_ort);
M3=blkdiag(Ri_tilde2*tij(:,1),Ri_tilde2*tij(:,2));
tangentAmbiguity3=M1*s+M2*vec(b)+M3*mu;
M=[M1 M2 M3];
v=[s;vec(b);mu];
tangentAmbiguity4=M*v;

disp('Check derivations of (38)')
disp(norm(tangentAmbiguity-tangentAmbiguity4,'fro'))

N=null(M);
disp('Nullspace of M')
disp(N)

%% checking (51)

Qx=align2d(Tij_tilde);
QxRitilde2Bot=Qx(3:4,:)*Ri_tilde2;
[U,~,~]=svd(QxRitilde2Bot,'econ');
c=U(:,2);

QLastRight=Qx(3:4,4)';

RbEst=procrustesRb(c,QLastRight');
d = 3;
lhs_51 = Qx' * blkdiag(eye(d), RbEst') * Qx * Ri_tilde2;
disp("lhs_51")
disp(lhs_51)

keyboard()

end % file function

function Qx=align2d_EoF(v)
Q=fliplr(orthComplement(v));
Qx=flipud(orthCompleteBasis(Q)');
end

function Qalign=align3d(v)
[U,S,V]=svd(v);
Qalign=fliplr(orthCompleteBasis(U(:,4)))';

end 

% function [RitildeEst1,RitildeEst2,Qx, RbEst]=recoverRitilde(Ritilde2,Tijtilde)
% Qx=align2d(Tijtilde);
% QxRitilde2Bot=Qx(3:4,:)*Ritilde2;
% [U,~,~]=svd(QxRitilde2Bot,'econ');
% c=U(:,2);
% 
% QLastRight=Qx(3:4,4)';
% 
% RbEst=procrustesRb(c,QLastRight');
% RitildeEst1=Qx'*blkdiag(eye(2),-RbEst')*Qx*Ritilde2;
% RitildeEst2=Qx'*blkdiag(eye(2),RbEst')*Qx*Ritilde2;
% end

function RbEst=procrustesRb(c,q)
[U,~,V]=svd(c*q');
diagmat = eye(size(U,1));
diagmat(end, end) = det(U*V');
RbEst=U*diagmat*V';
end
