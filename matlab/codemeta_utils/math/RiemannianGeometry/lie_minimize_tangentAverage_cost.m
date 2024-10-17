function fCost=lie_minimize_tangentAverage_cost(lf,y,yi,varargin)
LNorm=2;
Ny=size(y,3);

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'norm'
            ivarargin=ivarargin+1;
            LNorm=varargin{ivarargin};
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

fCost=zeros(Ny,1);
for iy=1:Ny
    fCost(iy)=sum(lf.dist(y(:,:,iy),yi).^LNorm)/LNorm;
end
