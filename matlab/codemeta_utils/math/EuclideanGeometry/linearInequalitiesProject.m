%project a vector onto the sets defined by each row of C*x>d
function Y=linearInequalitiesProject(x,C,d,varargin)
flagNormalized=false;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'normalized'
            flagNormalized=true;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

if ~flagNormalized
    [C,d]=linearInequalitiesNormalize(C,d);
end

nbInequalities=size(C,1);
Y=repmat(x,1,nbInequalities);
for iInequality=1:nbInequalities
    ci=C(iInequality,:);
    di=d(iInequality);
    pi=ci*x;
    
    if pi<di
        Y(:,iInequality)=x+(di-pi)*ci';
    end
    %note: if the inequality is satisfied, Y(:,iInequality) alread contains x
end
