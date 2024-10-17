function C=walks2indicatorVectors(E,EWalks,varargin)
E=E(:,1:2);

flagDirected=false;
%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'directed'
            flagDirected=true;
        case 'undirected'
            flagDirected=false;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

NWalks=size(EWalks,1);
NEdges=size(E,1);
C=zeros(NEdges,NWalks);
for iWalk=1:NWalks
    C(:,iWalk)=findEdges(E,EWalks(iWalk,:),flagDirected);
end

function c=findEdges(E,thisWalk,flagDirected)
c=false(size(E,1),1);
for k=1:length(thisWalk)-1
    c=c | (E(:,1)==thisWalk(k) & E(:,2)==thisWalk(k+1));
    if ~flagDirected
        c=c | (E(:,2)==thisWalk(k) & E(:,1)==thisWalk(k+1));
    end
end
