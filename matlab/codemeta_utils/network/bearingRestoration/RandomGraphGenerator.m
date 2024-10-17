function [E,x]=RandomGraphGenerator(n,m,varargin)
Range=100;
D=2;
%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'Range'
            ivarargin=ivarargin+1;
            Range=varargin{ivarargin};
        case 'dimension'
            ivarargin=ivarargin+1;
            D=varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

if m>nchoosek(n,2)
    m=nchoosek(n,2);
elseif m<n-1
    m=n-1;
end

x=rand([D,n])*Range;

E=[];
Set=1;
for idx=2:n
    node=randi(idx-1);
    E=[E;sort([idx node])];
end

comb=nchoosek(n,2)-(n-1);
V=nchoosek(1:n,2);
V=setdiff(V,E,'rows');

for idx=n:m
    randnum=randi(comb);
    v=V(randnum,:);
    E=[E;v];
    V=setdiff(V,E,'rows');
    comb=comb-1;
end
end