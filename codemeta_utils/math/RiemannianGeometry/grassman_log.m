function [H]=grassman_log(Y1,Y2,varargin)
method='gsvd';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'method'
            ivarargin=ivarargin+1;
            method=varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

%recursive call if Y2 contains multiple points
N2=size(Y2,3);
if(N2>1)
    H=zeros([size(Y1) N2]);
    for(iN2=1:N2)
        H(:,:,iN2)=grassman_log(Y1,Y2(:,:,iN2),varargin{:});
    end
else
    [n,p]=size(Y1);

    switch method
        case 'iterative'
            YY1=orthCompleteBasis(Y1);
            YY2=orthCompleteBasis(Y2);

            YY2=grassman_log_findOrthGrad(YY1,YY2,p,'linearInit');
            H=rot_log(YY1,YY2);
            H=H(:,1:p);
        case 'gsvd'
            H=logGsvd(Y1,Y2,n,p);
    end
end

function H=logGsvd(Y1,Y2,n,p)
A=Y1'*Y2;
B=Y2-Y1*Y1'*Y2;

[U,V,~,C,S] = gsvd(A,B);
theta=atan2(diag(S),diag(C));
H=V*[diag(theta); zeros(n-p,p)]*U';
