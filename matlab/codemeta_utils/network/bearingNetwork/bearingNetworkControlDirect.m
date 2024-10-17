%Compute control law for single-integrator model
%function dx=bearingNetworkControlDirect(EBearings,yBearings,ygBearings,funsBearings,varargin) 
%Inputs
%   EBearings       edge list for bearing vectors measurements
%   yBearings       bearing vectors corresponding to EBearings
%   ygBearings      desired values for yBearings
%   funsBearings    reshaping function structure for bearing terms
%Optional Inputs
%   'Ranges',ERanges,yRanges,ygRanges,nygRanges,funsRanges
%                   Add control from range measurements
%       ERanges         edge list for bearing vectors and range measurements
%       yRanges         bearing vectors corresponding to ERanges
%       ygRanges        desired values for yRanges
%       nyRanges        range measurements corresponding to ERanges
%       nygRanges       desired values for ygRanges
%       funsRanges      reshaping function structure for range terms
%
%   'Alpha',alpha   [1 x 2] vector with weights for the two types of terms
%
%   't',t           Current wall time (used only for time-sensitive
%                   operations)
%   'project',idxEdge
%                   Project the computed control for the two nodes
%                   determined by EBearings(idxEdge,:) on the orthogonal
%                   complement of the bearing direction for that edge
%   'leader',idxLeader,dxLeader
%                   Set dx to dxLeader for the nodes idxLeader
%   'preconditioner', T
%                   Multiply gradient by preconditioner T
%   'preconditionerDelay', tT
%                   Apply preconditioner only for t>tT
%Output
%   dx  control velocities
%
function dx=bearingNetworkControlDirect(EBearings,yBearings,ygBearings,funsBearings,varargin)
flagUseRanges=false;
flagProject=false;
flagSetLeader=false;
flagPreconditioner=false;
tPreconditioner=0;
t=0;

methodBearing='gradient';
alpha=1;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'ranges'
            flagUseRanges=true;
            ivarargin=ivarargin+1;
            ERanges=varargin{ivarargin};
            ivarargin=ivarargin+1;
            yRanges=varargin{ivarargin};
            ivarargin=ivarargin+1;
            ygRanges=varargin{ivarargin};
            ivarargin=ivarargin+1;
            nyRanges=varargin{ivarargin};
            ivarargin=ivarargin+1;
            nygRanges=varargin{ivarargin};
            ivarargin=ivarargin+1;
            funsRanges=varargin{ivarargin};
        case 'alpha'
            ivarargin=ivarargin+1;
            alpha=varargin{ivarargin};
        case 'methodbearing'
            ivarargin=ivarargin+1;
            methodBearing=varargin{ivarargin};
        case 'project'
            flagProject=true;
            ivarargin=ivarargin+1;
            idxEBearingProject=varargin{ivarargin};
        case 'leader'
            flagSetLeader=true;
            ivarargin=ivarargin+1;
            idxXLeader=varargin{ivarargin};
            ivarargin=ivarargin+1;
            dxLeader=varargin{ivarargin};
        case 'preconditioner'
            flagPreconditioner=true;
            ivarargin=ivarargin+1;
            T=varargin{ivarargin};
        case 'preconditionerdelay'
            flagPreconditioner=true;
            ivarargin=ivarargin+1;
            tPreconditioner=varargin{ivarargin};
        case 't'
            ivarargin=ivarargin+1;
            t=varargin{ivarargin};
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

if flagUseRanges && length(alpha)==1
    alpha=[1 alpha];
end

%compute pure bearing contribution
switch lower(methodBearing)
    case 'gradient'
        dx=-bearingNetworkCostGradient(EBearings,yBearings,ygBearings,funsBearings);
    case 'geodesic'
        dx=bearingNetworkField(EBearings,yBearings,ygBearings,'geodesic',funsBearings);
    case 'projection'
        dx=bearingNetworkField(EBearings,yBearings,ygBearings,'projection');
end
dx=alpha(1)*dx;

%compute bearing+range contribution if necessary
if flagUseRanges
    NNodes=max(EBearings(:));
    gRanges=bearingNetworkCostRangesGradient(ERanges,yRanges,ygRanges,nyRanges,nygRanges,funsRanges,'NNodes',NNodes);
    dx=dx-alpha(2)*gRanges;
end

%project on orthogonal complement of one of the bearings if necessary
if flagProject
    Py=orthComplementProjector(yBearings(:,idxEBearingProject));
    iNode=EBearings(idxEBearingProject,1);
    jNode=EBearings(idxEBearingProject,2);
    dx(:,[iNode jNode])=Py*dx(:,[iNode jNode]);
end

%set dx to the one provided for the leaders if necessary
if flagSetLeader
    if isa(dxLeader, 'function_handle')
        dx(:,idxXLeader)=dxLeader(t);
    else
        dx(:,idxXLeader)=dxLeader;
    end
end

%precondition dx, if the inverse preconditioner T is provided
if flagPreconditioner && t>tPreconditioner
    if size(T,1)==size(yBearings,1)
        %a preconditioner for each node
        for iNode=1:size(dx,2)
            dx(:,iNode)=T(:,:,iNode)*dx(:,iNode);
        end
    else
        %general preconditioner
        dx=reshape(T*dx(:),size(dx));
    end
end
