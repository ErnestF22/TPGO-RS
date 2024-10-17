function [v]=combIndexGenIdx(n,k,Idx)
%% This function generates the Combinatorial Numbering system series from ...
% Idx, such that the first value is 1:k and the last value is n-k+1:n
% and Idx ranges from 1 to nchoosek(n,k)
if Idx<1 || Idx>nchoosek(n,k)
    error('Idx out of bounds nchoosek(n,k)!!!');
end
N=Idx-1;
C=0:k-1;
for Kidx=k:-1:1
    if N>0
        C(Kidx)=C(Kidx)+1;
        while nchoosek(C(Kidx),Kidx)<=N
            C(Kidx)=C(Kidx)+1;
        end
        C(Kidx)=C(Kidx)-1;
        N=N-nchoosek(C(Kidx),Kidx);
    end
end
v=C+1;
end