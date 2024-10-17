%function a=vecsym(A)
%VECSYM generates the vector for A which linearizes the quadratic form
%x'*A*x, meaning
%   x'*A*x=vecsym(A)'*veronese(x,2)
function [a,idxPowers]=vecsym(A,idxPowers)
K=size(A,1);

if ~exist('idxPowers','var') || isempty(idxPowers)
    powers = exponent(2,K);
    m=size(powers,1);
    idxPowers=struct('idx',cell(m,1),'idxLen',cell(m,1));
    for ia=1:m
        idxPowers(ia).idx=find(powers(ia,:)>0);
        idxPowers(ia).idxLen=length(idxPowers(ia).idx);
    end
else
    m=length(idxPowers);
end

a=zeros(m,1);

for ia=1:m
    idx=idxPowers(ia).idx;
    switch idxPowers(ia).idxLen
        case 2
            a(ia)=A(idx(1),idx(2))+A(idx(2),idx(1));
        case 1
            a(ia)=A(idx,idx);
    end
end
