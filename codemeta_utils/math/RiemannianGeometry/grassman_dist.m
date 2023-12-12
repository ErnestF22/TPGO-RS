function d=grassman_dist(Y1,Y2)
H12=grassman_log(Y1,Y2);
N=size(Y2,3);
d=zeros(N,1);
for iN=1:N
    [U,S,V]=svd(H12(:,:,iN));
    d(iN)=norm(S);
end

% N2=size(Y2,3);
% d=zeros(1,N2);
% for iN2=1:N2
%     A=Y2(:,:,iN2)'*Y1;
%     B=Y2(:,:,iN2)-Y1*Y1'*Y2(:,:,iN2);
%     [~,~,~,C,S]=gsvd(A,B);
%     theta=atan2(diag(S),diag(C));
%     d(iN2)=norm(theta);
% end
 
