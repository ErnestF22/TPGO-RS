%Generates the walk tensor from edges
function EWalk=edges2walks(E,k,varargin)
flagSymmetrize=false;
if ~exist('k','var') || isempty(k)
    k=1;
end

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'symmetrize'
            flagSymmetrize=true;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

EWalk=E;
for ik=2:k
    NWalks=size(EWalk,1);
    EWalkNew=[];
    for iIdx=1:NWalks
        thisWalk=EWalk(iIdx,:);
        node=thisWalk(end);
        neighbors=E(E(:,1)==node,2)';
        if flagSymmetrize
            neighbors=union(neighbors, E(E(:,2)==node,1)');
        end
        EWalkNew=[EWalkNew; [repmat(thisWalk,length(neighbors),1) neighbors']];
    end
    EWalk=EWalkNew;
end

