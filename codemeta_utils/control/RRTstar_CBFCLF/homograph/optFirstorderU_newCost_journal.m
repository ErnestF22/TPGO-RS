%this function finds K1 for distance measurments:
% y = L-x and u = Ky+k_added
function Kres = optFirstorderU_newCost_journal(M,lambda,flag_loop)
%Wb:steering factor for cbf
%Wl:steering factor for clf
%Cb:cbf coefficoient
%Cl:clf coefficoient
%z:exit direction
%Ah:barrier matrix
%bh:barrier bias
%Au:velocity matrix
%bu:velocity bias
%Ax:convex cell matrix
%bx:convex cell bias
%Y: vectorized of landmarks
%d:dimestion


L = M(1).L;
I =kron(ones(size(L,2),1),eye(2));
L = reshape(L,[],1);
nL = length(L);
XE = kron(ones(1,nL/2),[28;25]);
XE = reshape(XE,[],1);
Ax1 = M(1).Ax;
Ah1 = M(1).Ah;
nSets = size(M,2);
nB = size(Ah1,1);
nC = size(Ax1,1);


cvx_begin
variables   Sly(nSets) sb(nSets) K(2*nSets,nL) k_added(2,nSets)...
    lB(nC*nSets,nB) lL(nC,nSets)


% c1(:,i) = K1*(reshape((y-xc1*ones(1,size(nL,2))),[],1))+k1_added;
% c2(:,i) = K2*(reshape((y-xc1*ones(1,size(nL,2))),[],1))+k2_added;
% c3(:,i) = K1*(reshape((y-xc2*ones(1,size(nL,2))),[],1))+k1_added;
% c4(:,i) = K2*(reshape((y-xc2*ones(1,size(nL,2))),[],1))+k2_added;
% (norm((c1-c2),2)+norm((c3-c4),2))
obj = 0;
for j=1:nSets-1
    y = M(j).L;
    z = M(j).z;
    P = eye(2)-lambda*z*z';
    xcorner1 = M(j).xc(:,1);
    xcorner2 = M(j).xc(:,2);
    c1 = K(2*(j-1)+1:2*j,:)*(reshape((y-xcorner1*ones(1,size(nL,2))),[],1))+k_added(:,j);%K1(y-x1)
    c2 = K(2*(j)+1:2*(j+1),:)*(reshape((y-xcorner1*ones(1,size(nL,2))),[],1))+k_added(:,j+1);%K2(y-x1)
    c3 = K(2*(j-1)+1:2*j,:)*(reshape((y-xcorner2*ones(1,size(nL,2))),[],1))+k_added(:,j);%K1(y-x2)
    c4 = K(2*(j)+1:2*(j+1),:)*(reshape((y-xcorner2*ones(1,size(nL,2))),[],1))+k_added(:,j+1);%K2(y-x2)
    obj = obj+(norm((c1-c2),2)+norm((c3-c4),2))+norm(P*c1-(xcorner2-xcorner1))+norm(P*c3-(xcorner1-xcorner2));
    %
end

minimize (obj)

subject to
for i=1:nSets
    Ax = M(i).Ax;
    bx = M(i).bx;
    Ah = M(i).Ah;
    bh = M(i).bh;
    xe = M(i).xe;
    z = M(i).z;
    y = M(i).L;
    L = reshape(y,[],1);
    Sb(:,i) = kron(ones(nB,1),sb(i));
    Cb = M(i).Cb;
    Cl = M(i).Cl;
    
    %constraints for safety:
    bx'*lB((i-1)*nC+1:i*nC,:) <= (Sb(:,i)+Ah*K(2*i-1:2*i,:)*L+Cb*bh)'+(Ah*k_added(:,i))'
    Ax'*lB((i-1)*nC+1:i*nC,:) == (Ah*K(2*i-1:2*i,:)*I-Ah*Cb)'
    
    %constraints for stability:
    if ~flag_loop
        if i==6
            K(2*i-1:2*i,:)*(L-XE)+k_added(:,i)==[0;0];
        else
            bx'*lL(:,i)<= (Sly(i)-z'*K(2*i-1:2*i,:)*L+Cl*z'*xe)'-z'*k_added(:,i)
            Ax'*lL(:,i) == (-z'*K(2*i-1:2*i,:)*I+z'*Cl)'
        end
    else
        bx'*lL(:,i)<= (Sly(i)-z'*K(2*i-1:2*i,:)*L+Cl*z'*xe)'-z'*k_added(:,i)
        Ax'*lL(:,i) == (-z'*K(2*i-1:2*i,:)*I+z'*Cl)'
    end
    
    
end

lB>=0
lL>=0
Sb<=0
Sly<=0
Sb>=-0.1
Sly>=-0.1

k_added<=2

cvx_end

for j=1:nSets
    Kres(j).K = K(2*j-1:2*j,:);
    Kres(j).added = k_added(:,j);
end
end



