function P=G2P(G,varargin)
methodAbsolutePoses='reference';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'references'
            methodAbsolutePoses='reference';
        case 'poses'
            methodAbsolutePoses='pose';
        case 'methodabsoluteposes'
            ivarargin=ivarargin+1;
            methodAbsolutePoses=lower(varargin{ivarargin});
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

switch methodAbsolutePoses
    case 'reference'
        G=invg(G);
    case 'pose'
        %nothing to do
    otherwise
        error('methodAbsolutePoses not valid')
end

P=RTK2P(G2R(G),G2T(G),eye(3));