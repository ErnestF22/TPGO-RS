%function [U,D]=udut(M)
% Given M, a symmetric positive-definit m x m matrix, U and D, modified
% Cholesky factors of M, are computed, such that U is a unit upper
% triangular matrix, D is a diagonal matrix, and M=UDU'

function [U,D]=udut(M)
m=size(M,1);

for(j=m:-1:1)
    for(i=j:-1:1)
        sigma=M(i,j);
        for(k=j+1:m)
            sigma=sigma-U(i,k)*D(k,k)*U(j,k);
        end
        if i==j
            D(j,j)=sigma;
            U(j,j)=1;
        else
            U(i,j)=sigma/D(j,j);
        end
    end
end