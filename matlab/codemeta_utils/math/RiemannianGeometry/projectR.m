%Project a vector of matrix on SO(3)
%function R=projectR(M,varargin)
%
%Optional arguments
%   'cholesky'              Use Cholesky factorization instead of SVD
%   'positiveDeterminant'   Flip sign if the determinant is negative
%

%%AUTORIGHTS%%

function R=projectR(M,varargin)
flagCholesky=false;
flagDeterminant=false;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch varargin{ivarargin}
        case 'cholesky'
            flagCholesky=true;
        case 'positiveDeterminant'
            flagDeterminant=true;
        otherwise
            error ('Unknown option for projectR')
    end
    ivarargin=ivarargin+1;
end

NM=size(M,3);
R=zeros(size(M));
for iM=1:NM
    if flagDeterminant && det(M(:,:,iM))<0
        M(:,:,iM)=-M(:,:,iM);
    end
    if(~flagCholesky)
        [U,~,V]=svd(M(:,:,iM));
        if(det(M(:,:,iM))<0)
            %flip sign of last column of V
            %(using the last row of U would be the same)
            V(:,3)=-V(:,3);
        end
        R(:,:,iM)=U*V';
    else
        A=chol(M(:,:,iM)*M(:,:,iM)');
        R(:,:,iM)=A'\M(:,:,iM);
        if(det(R(:,:,iM))<0)
            R(:,:,iM)=-R(:,:,iM);
        end
    end
end