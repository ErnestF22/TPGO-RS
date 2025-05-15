function Ri_est=RbRecovery(Ri_tilde2,Tij_tilde)
%Given Ri_tilde2 such that Ri_tilde2*tij*Lij-Tij_tilde, where Tij_tilde has
%last p-3 rows to zero, return Ri_est that satisfies the same equation, but
%has the last p-3 rows to zero
%If the inputs have multiple slices, apply the algorithm independently for
%each slice
nbPoses=size(Ri_tilde2,3);
d = size(Ri_tilde2, 2);
p = size(Ri_tilde2, 1);
if nbPoses>1
    Ri_est=zeros(size(Ri_tilde2));
    for iPose=1:nbPoses
        Ri_est(:,:,iPose)=RbRecovery(Ri_tilde2(:,:,iPose),Tij_tilde(:,:,iPose));
    end
else
    % base case, for single pose
    if norm(Tij_tilde(4:end,:),'fro')/numel(Tij_tilde(4:end,:))>1e-5
        error('Tij_tilde expected to have p-3 lines equal to zero')
    end
    Qx=align2d_nbPoses(Tij_tilde);
    Qbot=Qx(:,d+1:end)';
    Rcal_bot=Qx(3:end,:)*Ri_tilde2;
    Qbot_right=Qbot(:,d:end);
    [URCal_bot,~,~]=svd(Rcal_bot);
    Rcal_bot_N=URCal_bot(:,2:end);
    Rb_est=procrustes_R(Qbot_right',Rcal_bot_N);
    Ri_est=Qx'*blkdiag(eye(2),Rb_est')*Qx*Ri_tilde2;
end
%disp('Get rotation that aligns last rows of Ri_tilde2 to zero')
%Qy=fliplr(orthCompleteBasis(N2))';
%disp(Qy*Ri_tilde2)

%disp('Test to undo the invariance using a known Rb')
%Qx'*blkdiag(eye(2),Rb')*Qx*Qy*Ri_tilde2

end %file function


function R=procrustes_R(X,Y)
nb_dim=size(X,1);
[U,~,V]=svd(Y*X');
R=U*diag([ones(nb_dim-1,1);det(U*V')])*V';
end %function

