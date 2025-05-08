function POCProjectiveRecovery2_backup()
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

%Embed Ri, Tij in p-dimensional space
Tij=[Tij;zeros(dimStair-3,2)];
Ri=[Ri;zeros(dimStair-3,3)];

Rtilde=rot_randn(eye(dimStair));
%Rtilde=eye(dimStair);
Tij_tilde_pre=Rtilde*Tij;
Ri_tilde_pre=Rtilde*Ri;

disp('Align to get to the case where the last row is zero')
Qalign=align3d(Tij_tilde_pre);
%Qalign=fliplr(orthCompleteBasis(null(Tij_tilde_pre')))';
Tij_tilde=Qalign*Tij_tilde_pre;
Ri_tilde=Qalign*Ri_tilde_pre;


disp('Check consistency between all generated data with high dimension (28)')
disp(norm(Ri_tilde*tij*Lij-Tij_tilde,'fro'))

disp('Test computation of Qx')
disp(['It should align the last ' num2str(dimStair-2) ' rows of Tijtilde to zero'])
Qx=align2d_EoF(Tij_tilde);

disp(Qx*Tij_tilde)


% Generate a random Rb and the associated Qb
Rb=rot_randn(eye(dimStair-2));
Qb=blkdiag(eye(2),Rb);

disp('Check that Qx''*Qb*Qx leaves Tijtilde invariant')
disp(norm(Qx'*Qb*Qx*Tij_tilde-Tij_tilde))
disp(norm(Ri_tilde*tij*Lij-Qx'*Qb*Qx*Tij_tilde))

%Eq. (50) defining \doubletilde{R}_i on Overleaf
disp('And that the new Ritilde generated from the ambiguity')
disp('still gives the same measurements.')
Ri_tilde2=Qx'*Qb*Qx*Ri_tilde;
disp(norm(Ri_tilde2*tij*Lij-Tij_tilde))

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
tangentAmbiguity2=vec(Ri_tilde2*hat3(s)*tij*Lij+Ri_tilde2_ort*b*tij*Lij+Ri_tilde2*tij*diag(mu));
M1=zeros(dimStair*2,3);
I3=eye(3);
for iM1=1:3
    M1(:,iM1)=vec(Ri_tilde2*hat3(I3(:,iM1))*tij*Lij);
end
M2=kron((tij*Lij)',Ri_tilde2_ort);
M3=blkdiag(Ri_tilde2*tij(:,1),Ri_tilde2*tij(:,2));
%tangentAmbiguity3=M1*s+M2*vec(b)+M3*mu;
M=[M1 M2 M3];
v=[s;vec(b);mu];
tangentAmbiguity4=M*v;

disp('Check derivations of (38)')
disp(norm(tangentAmbiguity-tangentAmbiguity4,'fro'))

N=null(M);
disp('Last two rows of nullspace of M')
disp('If zero, it means that there is no ambiguity involving the lambdas')
disp(N(end-1:end,:))


% To get recovery
disp('Any matrix of the form Qx''*blkdiag(eye(2),Rb)*Qx can change Ri_tilde2')
disp('without changing Tij_tilde')
Rbrand=rot_randn(eye(3));
disp(norm(Qx'*blkdiag(eye(2),Rbrand)*Qx*Ri_tilde2*tij*Lij-Ri_tilde2*tij*Lij,'fro'))

%equation (51) on Overleaf
disp('Now we want to find Rb such that Qx''*blkdiag(eye(2),Rb)*Qx*Ri_tilde2 has last p-3 rows to zero')

disp('Generate splits Qbot, Rcal, Rcal_top, Rcal_bot, Qbot_left, Qbot_right')
pm3=dimStair-3;
Qbot=Qx(:,end-(dimStair-3)+1:end)';
assert(size(Qbot,1)==dimStair-3)
Rcal=Qx*Ri_tilde2;

Rcal_top=Rcal(1:2,:);
Rcal_bot=Rcal(3:end,:);
pm2=dimStair-2;
assert(size(Rcal_bot,1)==pm2)

Qbot_left=Qbot(:,1:2);
Qbot_right=Qbot(:,3:end);

disp('Verify that Qbot_left*Rcal_top+Qbot_right*Rb*Rcal_bot is equal to')
disp('the last p-3 rows of Qx''*blkdiag(eye(2),Rb)*Qx*Ri_tilde2 for random Rb')
disp('(LHS of (53) on Overleaf)')
A1=Qbot_left*Rcal_top+Qbot_right*Rbrand'*Rcal_bot;
A2=[zeros(2,pm2),eye(2)]*Qx'*blkdiag(eye(2),Rbrand')*Qx*Ri_tilde2;
disp(norm(A1-A2,'fro'))
disp('Which is the same as Qbot_right*Rb*Rcal_bot, because Qbot_left is zero')
disp(norm(Qbot_right*Rbrand'*Rcal_bot-A2,'fro'))

disp('Since RHS of (53) is zero, check that Rcal_bot has p-2 sval equal to zero')
sRcal_bot=flipud(svd(Rcal_bot));
disp(sRcal_bot(1:pm3))


[U,~,~]=svd(Rcal_bot);
Rcal_bot_N=U(:,end-1:end);
Rb_est=procrustes_R(Qbot_right',Rcal_bot_N);

disp('Check Qbot_right*Rb*Rcal_bot==0 when Rb is found by Procrustes')
disp(norm(Qbot_right*Rb_est'*Rcal_bot,'fro'))

disp('Last p-2 rows of Ri_est=Qx''*blkdiag(eye(2),Rb_est'')*Qx*Ri_tilde2 should be equal to zero')
Ri_est=Qx'*blkdiag(eye(2),Rb_est')*Qx*Ri_tilde2;
disp(Ri_est)

disp('Ri_est should fit the measurements')
disp(Ri_est*tij*Lij-Tij_tilde)

disp('Check that output of RbRecovery matches calculations above')
Ri_est2=RbRecovery(Ri_tilde2,Tij_tilde);
disp(norm(Ri_est-Ri_est2,'fro'))


%keyboard()


function Ri_est=RbRecovery(Ri_tilde2,Tij_tilde)
%Given Ri_tilde2 such that Ri_tilde2*tij*Lij-Tij_tilde, where Tij_tilde has
%last p-3 rows to zero, return Ri_est that satisfies the same equation, but
%has the last p-3 rows to zero
if norm(Tij_tilde(4:end,:),'fro')/numel(Tij_tilde(4:end,:))>1e-5
    error('Tij_tilde expected to have p-3 lines equal to zero')
end
Qx=align2d_EoF(Tij_tilde);
Qbot=Qx(:,4:end)';
Rcal_bot=Qx(3:end,:)*Ri_tilde2;
Qbot_right=Qbot(:,3:end);
[URCal_bot,~,~]=svd(Rcal_bot);
Rcal_bot_N=URCal_bot(:,end-1:end);
Rb_est=procrustes_R(Qbot_right',Rcal_bot_N);
Ri_est=Qx'*blkdiag(eye(2),Rb_est')*Qx*Ri_tilde2;

%disp('Get rotation that aligns last rows of Ri_tilde2 to zero')
%Qy=fliplr(orthCompleteBasis(N2))';
%disp(Qy*Ri_tilde2)

%disp('Test to undo the invariance using a known Rb')
%Qx'*blkdiag(eye(2),Rb')*Qx*Qy*Ri_tilde2

function R=procrustes_R(X,Y)
nb_dim=size(X,1);
[U,~,V]=svd(Y*X');
R=U*diag([ones(nb_dim-1,1);det(U*V')])*V';


function Qx=align2d_EoF(v)
Q=fliplr(orthComplement(v));
Qx=flipud(orthCompleteBasis(Q)');

function Qalign=align3d(v)
[U,S,V]=svd(v);
Qalign=fliplr(orthCompleteBasis(U(:,4:end)))';