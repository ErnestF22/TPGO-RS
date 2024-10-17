function yInit=lie_minimize_tangentAverage_getInitialPoint(lf,yi,varargin)
LNorm=2;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'method'
            ivarargin=ivarargin+1;
            method=varargin{ivarargin};
        case 'norm'
            ivarargin=ivarargin+1;
            LNorm=varargin{ivarargin};
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

if strcmpi(method,'twominaverage') && size(yi,3)==1
    method='min';
end

switch lower(method)
    case 'min'
        yInit=lie_minimize_tangentAverage_selectMininum(lf,yi,'norm',LNorm);
    case 'twominaverage'
        yInitPre=lie_minimize_tangentAverage_selectMininum(lf,yi,'norm',LNorm,'nmin',2);
        yInit=lf.exp(yInitPre(:,:,1),0.5*lf.log(yInitPre(:,:,1),yInitPre(:,:,2)));
end
