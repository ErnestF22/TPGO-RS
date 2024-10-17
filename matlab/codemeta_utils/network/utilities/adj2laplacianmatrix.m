%Computes the Laplacian matrix of a graph given the adjacency matrix
function L=adj2laplacianmatrix(A,varargin)
flagNormalized=false;
%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'normalized'
            flagNormalized=true;
%         case 'methodabsoluteposes'
%             ivarargin=ivarargin+1;
%             methodAbsolutePoses=lower(varargin{ivarargin});
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

d=sum(A,2);
D=diag(d);

if ~flagNormalized
    L=D-A;
else
    normD=diag(1./sqrt(d));
    L=eye(size(A))-normD*A*normD;
end
