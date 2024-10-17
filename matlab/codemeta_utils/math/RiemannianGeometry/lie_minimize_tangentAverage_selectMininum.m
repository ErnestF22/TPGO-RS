function y=lie_minimize_tangentAverage_selectMininum(lf,yi,varargin)
LNorm=2;
NMin=1;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'norm'
            ivarargin=ivarargin+1;
            LNorm=varargin{ivarargin};
        case 'nmin'
            ivarargin=ivarargin+1;
            NMin=varargin{ivarargin};
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

fCost=lie_minimize_tangentAverage_cost(lf,yi,yi,'norm',LNorm);
[~,idxSort]=sort(fCost);
y=yi(:,:,idxSort(1:NMin));
