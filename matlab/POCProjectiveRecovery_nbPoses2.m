function POCProjectiveRecovery_nbPoses2()
nbPoses=10;
d = 3;
Ti=randn(d,1,nbPoses);
Tj=randn(d,2,nbPoses);
Ri=rot_randn(eye(d),[],nbPoses);
Tij=Tj-Ti;
[tij,lij]=cnormalize(multiprod(multitransp(Ri),Tij));
Lij=multidiag(lij);

disp('Check consistency between all generated data (28)')
disp(normvec(multiprodN(Ri,tij,Lij)-Tij))

% Go to higher dimension
dimStair=9;

%Embed Ri, Tij in p-dimensional space
low_deg = 2;
Tij=[Tij;zeros(dimStair-d,low_deg,nbPoses)];
Ri=[Ri;zeros(dimStair-d,d,nbPoses)];

Rtilde=rot_randn(eye(dimStair));
%Rtilde=eye(dimStair);
Tij_tilde_pre=multiprod(Rtilde,Tij);
Ri_tilde_pre=multiprod(Rtilde,Ri);

disp('Align to get to the case where the last row is zero')
Qalign=align3d(Tij_tilde_pre);
%Qalign=fliplr(orthCompleteBasis(null(Tij_tilde_pre')))';
Tij_tilde=multiprod(Qalign,Tij_tilde_pre);
Ri_tilde=multiprod(Qalign,Ri_tilde_pre);


disp('Check consistency between all generated data with high dimension (28)')
disp(normvec(multiprodN(Ri_tilde,tij,Lij)-Tij_tilde))

disp('Test computation of Qx')
disp(['It should align the last ' num2str(dimStair-2) ' rows of Tijtilde to zero'])
Qx=align2d_nbPoses(Tij_tilde);

disp(multiprod(Qx,Tij_tilde))

% Generate a random Rb and the associated Qb
Rb=rot_randn(eye(dimStair-2),[],nbPoses);
Qb=multiblkdiag(repmat(eye(2),[1 1 nbPoses]),Rb);

disp('Check that Qx''*Qb*Qx leaves Tijtilde invariant')
disp(normvec(multiprodN(multitransp(Qx),Qb,Qx,Tij_tilde)...
    -Tij_tilde))
disp(normvec(multiprodN(Ri_tilde,tij,Lij)...
    -multiprodN(multitransp(Qx),Qb,Qx,Tij_tilde)))

%Eq. (50) defining \doubletilde{R}_i on Overleaf
disp('And that the new Ritilde generated from the ambiguity')
disp('still gives the same measurements.')
Ri_tilde2=multiprodN(multitransp(Qx),Qb,Qx,Ri_tilde);
disp(normvec(multiprodN(Ri_tilde2,tij,Lij)-Tij_tilde))

disp('But Ritilde2 does not have the last row to zero')
disp(Ri_tilde2)

if nbPoses==1
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
    % tangentAmbiguity2=vec(Ri_tilde2*hat3(s)*tij*Lij+Ri_tilde2_ort*b*tij*Lij+Ri_tilde2*tij*diag(mu));
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
end

% To get recovery
disp('Any matrix of the form Qx''*blkdiag(eye(2),Rb)*Qx can change Ri_tilde2')
disp('without changing Tij_tilde')
Rbrand=rot_randn(eye(3),[],nbPoses);
% disp(normvec(...
%     multiprodN(multitransp(Qx),multiblkdiag(repmat(eye(2),[1 1 nbPoses]),Rbrand),Qx,Ri_tilde2,tij,Lij)...
%     -multiprodN(Ri_tilde2,tij,Lij)...
%     ))

if nbPoses==1
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
end

Ri_est2=RbRecovery(Ri_tilde2,Tij_tilde);
disp('Last p-2 rows of Ri_est2 should be equal to zero')
disp(Ri_est2)

disp('Ri_est2 should fit the measurements')
disp(normvec(multiprodN(Ri_est2,tij,Lij)-Tij_tilde))

% disp('det(Ri_est(1:d, 1:d))')
% disp(det(Ri_est(1:d, 1:d)))

disp('det(Ri_est2(1:d, 1:d))')
disp(det(Ri_est2(1:d, 1:d)))

end %file function

