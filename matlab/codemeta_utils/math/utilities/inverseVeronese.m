function [x,powers] = inverseVeronese(y,n,K,powers)

[Mn,N] = size(y);

if(nargin<4)
    powers = exponent(n,K);
end

[l,m]=find(powers==n);

%naive inversion, can be improved
signs=[-1 1];
switch n
    case 2
        x=abs(y(l,:).^(1/n));
        sx=zeros(K,1);
        for i=1:N
            sy=sign(y(:,i));
            for isigns=1:2;
                sx1=signs(isigns);
                sx(isigns,1)=sx1;
                sx(isigns,2:K)=sy(2:K)/sx(1);
                err(isigns)=sum((y-veronese(sx(isigns,:)'.*x(:,i),n)).^2);
            end
            if(err(1)<err(2))
                x(:,i)=sx(1,:)'.*x(:,i);
            else
                x(:,i)=sx(2,:)'.*x(:,i);
            end
        end
                    
                
    case 3
        x=abs(y(l,:).^(1/n)).*sign(y(l,:));
    otherwise
        error('Inverse of veronese not implemented for this value of n')
end


% numerical optimization
% fcost=inline('sum((veronese(x,n)-y).^2)','x','n','y');
% 
% for(idx=1:N)
%     x(:,idx)=fminunc(@(x) fcost(x,n,y(:,idx)), x(:,idx));
% end
