function C=xcorrCoeff(A,B,varargin)
if ~exist('B','var') || isempty(B)
    B=A;
end

if size(A,1)~=size(B,1)
    error('A and B have incompatible number of rows')
end

outputType='matrix';

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case {'vector','matrix'}
            outputType=lower(varargin{ivarargin});
        case 'outputtype'
            ivarargin=ivarargin+1;
            outputType=lower(varargin{ivarargin});
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

ACentered=removeMean(A);
BCentered=removeMean(B);

ANormSq=sum(ACentered.^2);
BNormSq=sum(BCentered.^2);

switch outputType
    case 'matrix'
        C=(ACentered'*BCentered)./sqrt(ANormSq'*BNormSq);
    case 'vector'
        if size(A,2)~=size(B,2)
            error('A and B have incompatible number of columns')
        end
        
        C=sum(ACentered.*BCentered)./sqrt(ANormSq.*BNormSq);
end

function ACentered=removeMean(A)
ACentered=A-ones(size(A,1),1)*mean(A);
