%Evaluate Lagrangian as a function of z
function c=admmMedoids_auxiliaryCost(muijk,lambdaijk,rho,varargin)

zEval=muijk;
flagDimSelect=false;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'zeval'
            ivarargin=ivarargin+1;
            zEval=varargin{ivarargin};
        case 'dim'
            ivarargin=ivarargin+1;
            dim=varargin{ivarargin};
            flagDimSelect=true;
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

if flagDimSelect
    muijk=muijk(dim,:);
    lambdaijk=lambdaijk(dim,:);
    zEval=zEval(dim,:);
end

c=sum(lambdaijk,2)'*zEval;
for iNeighbor=1:2
    e=zEval-muijk(:,iNeighbor);
    c=c+rho*sum(abs(e));
end
