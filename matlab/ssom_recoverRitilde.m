function [RitildeEst1,RitildeEst2,Qx, RbEst]=ssom_recoverRitilde(Ritilde2,Tijtilde)
Qx=align2d(Tijtilde);
QxRitilde2Bot=Qx(3:end,:)*Ritilde2;
[U,~,~]=svd(QxRitilde2Bot,'econ');
c=U(:,2);

QLastRight=Qx(3:end,end);

RbEst=ssom_procrustesRb(c,QLastRight);
RitildeEst1=Qx'*blkdiag(eye(2),-RbEst')*Qx*Ritilde2;
RitildeEst2=Qx'*blkdiag(eye(2),RbEst')*Qx*Ritilde2;
end





function RbEst=ssom_procrustesRb(c,q)
[U,~,V]=svd(c*q');
tmp = eye(size(U));
tmp(end,end) = det(U*V');
RbEst=U*tmp*V';
end