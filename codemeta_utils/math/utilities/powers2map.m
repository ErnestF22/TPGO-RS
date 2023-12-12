%function map=powers2map(K)
%Creates a matrix MAP where the extries (i,j) indicates what is the entry
%in the vector veronese(x,2) corresponding to x_i*x_j
function [map,m]=powers2map(K)
powers=exponent(2,K);
m=size(powers,1);

map=zeros(K);

for(it=1:m);
    idx=find(powers(it,:)>0);
    switch(length(idx))
        case 1
            map(idx,idx)=it;
        case 2
            map(idx(1),idx(2))=it;
            map(idx(2),idx(1))=it;
    end
end

