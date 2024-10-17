function [E,w]=adj2edges(A,varargin)
methodDirection='directed';
flagWeighted=nargout>1;
%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case {'undirected','oriented','directed'}
            methodDirection=lower(varargin{ivarargin});
        case 'weighted'
            flagWeighted=true;
        otherwise
            error('MATLAB:ArgumentInvalid',['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

N=size(A,1);

%make sure A is simmetric if undirected or oriented
switch methodDirection
    case 'directed'
        %nothing to do
    case {'undirected','oriented'}
        A=(A+A')/2;
    otherwise
        error('methodDirection invalid')
end
        
[I,J]=find(A~=0);
E=[I J];
idxValid=E(:,1)<E(:,2);
E=E(idxValid,:);

switch methodDirection
    case 'directed'
        E=[E;fliplr(E)];
    case {'undirected','oriented'}
        %nothing to do
    otherwise
        error('methodDirection invalid')
end

if flagWeighted
    w=A(sub2ind(size(A),E(:,1),E(:,2)));
end
