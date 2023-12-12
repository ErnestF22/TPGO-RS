%Computes squared euclidean distance for matching
function [D,DIndicator]=matchingEuclideanDistMatrix(X,varargin)
flagMembershipPrior=false;
flagNBest=false;
%parse optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'membershipprior'
            ivarargin=ivarargin+1;
            membershipPrior=varargin{ivarargin};
            flagMembershipPrior=true;
        case 'nbest'
            ivarargin=ivarargin+1;
            NBest=varargin{ivarargin};
            flagNBest=true;
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

if flagNBest
    if ~flagMembershipPrior
        D=euclideanDistMatrix(X);
        [i,j,v,m,n]=buildSparseElements(D,NBest,'removeDiagonal');
        D=sparse(i,j,v,m,n);
        DIndicator=sparse(i,j,true(size(v)),m,n);
    else
        m=length(membershipPrior);
        iD=[];
        jD=[];
        vD=[];
        classes=unique(membershipPrior);
        for iClass=classes
            for jClass=classes
                iIdx=find(membershipPrior==iClass);
                jIdx=find(membershipPrior==jClass);
                Dij=euclideanDistMatrix(X(:,iIdx),X(:,jIdx));
                if iClass==jClass
                    opts={'removeDiagonal'};
                else
                    opts={};
                end
                [i,j,v]=buildSparseElements(Dij,NBest,opts{:});
                iiIdx=iIdx(i);
                jjIdx=jIdx(j);
                iD=[iD iiIdx];
                jD=[jD jjIdx];
                vD=[vD;v];
            end
        end
        D=sparse(iD,jD,vD,m,m);
        DIndicator=sparse(iD,jD,true(size(vD)),m,m);
    end
else
    D=euclideanDistMatrix(X);
    DIndicator=full(size(D));
end

function [i,j,v,m,n]=buildSparseElements(D,NBest,varargin)
flagRemoveDiagonal=false;
methodNBest='sort';

%parse optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'removediagonal'
            flagRemoveDiagonal=true;
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

%size of the matrix
[m,n]=size(D);

NSelect=min(NBest,n);
if flagRemoveDiagonal && NSelect~=n
    NSelect=NSelect+1;
end


%find NSelect elements
switch methodNBest
    case 'sort'
        [DSort,idxSort]=sort(D,1);
        i=vec(idxSort(1:NSelect,:));
        v=vec(DSort(1:NSelect,:));
    case 'min'
        [DSort,idxSort]=mink(D,NSelect);
        i=idxSort(:);
        v=DSort(:);
    otherwise
        error('Method for selecting top N values not recognized')
end
j=vec(ones(NSelect,1)*(1:n));

if flagRemoveDiagonal
    %remove elements on the diagonal
    flagValid=i~=j;
    i=i(flagValid);
    j=j(flagValid);
    v=v(flagValid);
end

