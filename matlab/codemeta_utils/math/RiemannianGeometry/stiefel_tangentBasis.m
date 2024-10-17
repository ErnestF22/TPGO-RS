%generates a basis for the tangent space of Stiefel(n,p) at Y
function X=stiefel_tangentBasis(Y)
[n,p]=size(Y);
%X=zeros(n,p,(n-p)*p+p*(p-1)/2);
np=(n-p)*p+p*(p-1)/2;
X=zeros(n,p*np);

YY=orthCompleteBasis(Y);
cnt=1;

%skew symmetric part (first p by p block)
for ii=1:p 
    for jj=ii+1:p 
        X(ii,jj+(cnt-1)*p)=1;
        X(jj,ii+(cnt-1)*p)=-1;
        cnt=cnt+1;
    end
end

%other part (lower n-p by p block)
for ii=p+1:n 
    for jj=1:p 
        X(ii,jj+(cnt-1)*p)=1;
        cnt=cnt+1;
    end
end

X=YY*X;
X=reshape(X,[n,p,np]);

% %skew symmetric part (first p by p block)
% for(ii=1:p)
%     for(jj=ii+1:p)
%         X(ii,jj,cnt)=1;
%         X(jj,ii,cnt)=-1;
%         X(:,:,cnt)=YY*X(:,:,cnt);
%         cnt=cnt+1;
%     end
% end
% 
% %other part (lower n-p by p block)
% for(ii=p+1:n)
%     for(jj=1:p)
%         X(ii,jj,cnt)=1;
%         X(:,:,cnt)=YY*X(:,:,cnt);
%         cnt=cnt+1;
%     end
% end
