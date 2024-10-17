%function testNetworkDisplay(S,varargin)
%If S is a 4x4xNCameras matrix, display cameras corresponding to each one
%of the poses S(:,:,iCamera). If S is a t_node structure, display cameras
%corresponding to each one of the poses t_node().gitruth (or another member
%if the option 'Member' is specified). Also, the matrix A is obtained from
%the fields t_node().aij.
%
%Optional arguments
%   'Points',X              Add the 3-D points X to the plot
%   'PointStyle',s          Style to use to draw the points (see PLOT)
%   'flagDisplayPoints',b   Turn on (b=true) or off (b=false) the
%                           visualization of the 3-D points (works only in
%                           conjunction with 'Points', default on)
%   'AdjMatrix',A           Pass an adjacency matrix and turn on
%                           visualization of edges
%   'flagDisplayEdges',b    Turn on/off visualization of edges (default on)
%   'DisplayEdges'          Same as 'flagDisplayEdges',true
%   'Member',member         If a t_node structure is passed, use
%                           t_node().(member) to get the camera reference
%                           frames. If a cell array is passed, all the
%                           specified members are split into ([member 'R')
%                           and ([member T]), but by default only the first
%                           one is used. The default is 'gitruth'.Has no
%                           effect if S is an array
%   'SplitMembers', memberR, memberT
%                           Build poses from the specified separate members
%                           for rotations and translations. This options is
%                           most useful when combined with the option
%                           'member'. The defaults are [member{1} 'R'] and
%                           [member{1} 'T']. Has no effect if S is an array
%   'EdgesStyle',style      Plot style for drawing the edges
%   'OptionsDrawCamera',opt Cell array to pass as options to
%                           draw3dcameraFromPose
%   'FlipZ'                 Adds 'FlipZ' to the options for
%                           draw3dcameraFromPose
%   'Estimated'             Same as 'OptionsDrawCamera',
%                           {'Color1','red','Color2','red'},
%                           'EdgeStyle','r','DisplayEdges'
%
%Note: everything regarding members is not available if S is an array.
%
%See also draw3dcameraFromPose, draw3dcameraFromAxesAndCenter

%%AUTORIGHTS%%

function testNetworkDisplay(S,varargin)
%3D points
X=[];                       %coordinates
flagDisplayPoints=true;     %flag
pointsStyle='xk';           %marker style
%Edges
E=[];                       %edge list
flagDisplayEdges=true;      %flag
edgesStyle='b';             %line style
%Geotags
flagDisplayGeoTags=true;    %flag
geoTagStyle={'s','MarkerSize',10,'MarkerFaceColor',[0 0 1]};

%Cameras
flagDisplayCameras=true;    %flag
green=[15934	35723	14392]/65535/0.6;           %camera color
optionsDrawCamera={'Color1',green,'Color2',green};  %options to pass to drawCamera
flagDisplayCameraNumber=false;

member='gitruth';
memberR='gitruthR';
memberT='gitruthT';
flagSetSplitMembers=false;

methodAbsolutePoses='reference';    %how to interpret poses

%Input checking
if ~isstruct(S) && (~(size(S,1)==4 || size(S,1)==3) || size(S,2)~=4)
    error('S must be a t_node structure array or a 4x4xNCameras matrix')
end

if isstruct(S)
    %S is a t_node structure
    t_node=S;
else
    %S is a 4x4xN matrix of poses
    G=S;
    if size(G,1)==3
        G=[G; zeros(1,3,size(G,3)) ones(1,1,size(G,3))];
    end
    methodAbsolutePoses='pose';
end
scale=1;

%get points from structure if necessary
if isfield(t_node,'X')
    X=t_node.X;
end

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'scale'
            ivarargin=ivarargin+1;
            scale=varargin{ivarargin};
        case 'points'
            ivarargin=ivarargin+1;
            X=varargin{ivarargin};
        case 'flagdisplaypoints'
            ivarargin=ivarargin+1;
            flagDisplayPoints=varargin{ivarargin};
        case 'pointsstyle'
            ivarargin=ivarargin+1;
            pointsStyle=varargin{ivarargin};
        case 'adjmatrix'
            ivarargin=ivarargin+1;
            E=adj2edges(varargin{ivarargin});
            flagDisplayEdges=true;
        case 'flagdisplayedges'
            ivarargin=ivarargin+1;
            flagDisplayEdges=varargin{ivarargin};
        case 'displayedges'
            flagDisplayEdges=true;
        case 'edgesstyle'
            ivarargin=ivarargin+1;
            edgesStyle=varargin{ivarargin};
        case 'member'
            ivarargin=ivarargin+1;
            member=varargin{ivarargin};
        case 'flagdisplaycameras'
            ivarargin=ivarargin+1;
            flagDisplayCameras=varargin{ivarargin};
        case 'flagdisplaycameranumber'
            ivarargin=ivarargin+1;
            flagDisplayCameraNumber=varargin{ivarargin};
        case 'edgesonly'
            flagDisplayCameras=false;
            flagDisplayEdges=true;
        case 'optionsdrawcamera'
            ivarargin=ivarargin+1;
            if ~iscell(varargin{ivarargin})
                newArgs=cell(1);
                newArgs{1}=varargin{ivarargin};
            else
                newArgs=varargin{ivarargin};
            end
            optionsDrawCamera=[optionsDrawCamera newArgs];
        case 'flipz'
            optionsDrawCamera=[optionsDrawCamera 'flipz'];
        case 'estimated'
            red=[65535	8567	0]/65535;
            optionsDrawCamera=[optionsDrawCamera 'Color1' red 'Color2' red];
            edgesStyle='r';
            flagDisplayEdges=true;
            flagDisplayGeoTags=false;
        case 'references'
            methodAbsolutePoses='reference';
        case 'poses'
            methodAbsolutePoses='pose';
        case 'methodabsoluteposes'
            ivarargin=ivarargin+1;
            methodAbsolutePoses=varargin{ivarargin};
        case 'splitmembers'
            ivarargin=ivarargin+1;
            memberR=varargin{ivarargin};
            ivarargin=ivarargin+1;
            memberT=varargin{ivarargin};
            flagSetSplitMembers=true;
           
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end


%convert values to cell arrays if necessary
if ~iscell(edgesStyle)
    edgesStyle={edgesStyle};
end

if ~iscell(pointsStyle)
    pointsStyle={pointsStyle};
end

if ~iscell(geoTagStyle)
    geoTagStyle={geoTagStyle};
end

if ~iscell(member)
    member={member};
end

%set default value of name of members for splitted rotation and
%translation part
if ~flagSetSplitMembers
    memberR=[member{1} 'R'];
    memberT=[member{1} 'T'];
end

if isstruct(S)
    %get structure type
    structType=testNetworkDetectStructType(t_node);

    %note, we already set t_node=S
    
    %if necessary, pass from poses to references
    switch structType
        case 'single'
            N=t_node.NNodes;
            switch methodAbsolutePoses
                case 'reference'
                    %nothing to do
                case 'pose'
                    for iMember=1:length(member)
                        t_node.(member{iMember})=invg(t_node.(member{iMember}));
                    end
                otherwise
                    error('Absolute pose method specification not valid')
            end
            
        case 'array'
            N=length(t_node);
            G=zeros(4,4,N);

            for iCamera=1:N
                for iMember=1:length(member)
                    switch methodAbsolutePoses
                        case 'reference'
                            %nothing to do
                        case 'pose'
                            t_node(iCamera).(member{iMember})=invg(t_node(iCamera).(member{iMember}));
                        otherwise
                            error('Absolute pose method specification not valid')
                    end
                end
            end
    end
    if isempty(E)
        E=testNetworkGetEdges(S);
    end
    
    %split members, so that we can access them
    for iMember=1:length(member)
        t_node=splitgi(t_node,member{iMember},[member{iMember} 'R'], [member{iMember} 'T']);
    end
    
    %assemble references from splitted members
    switch structType
        case 'single'
            G=RT2G(t_node.(memberR), t_node.(memberT));
        case 'array'
            for iNode=1:N
                G(:,:,iNode)=RT2G(t_node(iNode).(memberR), t_node(iNode).(memberT));
            end
    end
else
    N=size(G,3);
    switch methodAbsolutePoses
        case 'reference'
            %do nothing
        case 'pose'
            G=invg(G);
        otherwise
            error('Absolute pose method specification not valid')
    end
end

%save "hold" status of the plot
flagHold=ishold();

%scale data
G(1:3,4,:)=scale*G(1:3,4,:);
if ~isempty(X)
    X=scale*X;
end

%display cameras
if flagDisplayCameras
    for iCamera=1:N
        draw3dcameraFromAxesAndCenter(G(1:3,1:3,iCamera),G(1:3,4,iCamera),optionsDrawCamera{:});
        hold on
        if flagDisplayCameraNumber
            text(G(1,4,iCamera),G(2,4,iCamera),G(3,4,iCamera),num2str(iCamera));
        end
    end
end
view(3)

%display 3-D points
if ~isempty(X) && flagDisplayPoints
    plot3(X(1,:),X(2,:),X(3,:),pointsStyle{:})
    ax=axis();
    ax=[min(ax(1),min(X(1,:))) max(ax(2),max(X(1,:)))...
        min(ax(3),min(X(2,:))) max(ax(4),max(X(2,:)))...
        min(ax(5),min(X(3,:))) max(ax(6),max(X(3,:)))];
    axis(ax);
end

%display edges
if ~isempty(E) && flagDisplayEdges
    %get camera centers in world coordinates
    T=G2T(G);
    plotArrows(T(:,E(:,1)),T(:,E(:,2)),edgesStyle{:})
    
%     %get (non-duplicated) edges
%     A(sub2ind([N,N],1:N,1:N))=0;        %remove self-loops
%     [I,J]=find((abs(A)+abs(A'))~=0);    %symmetrize
%     flagUnique=I<J;                     %take only one direction
%     E=[I(flagUnique) J(flagUnique)];
%     %draw
%     for iEdge=1:size(E,1)
%         TT=[T(:,E(iEdge,1)) T(:,E(iEdge,2))];
%         plot3(TT(1,:),TT(2,:),TT(3,:),edgesStyle{:});
%         hold on
%     end
end

if flagDisplayGeoTags && isstruct(S)
    [geoTagIdx,geoTag]=testNetworkGetGeoTags(t_node);
    if ~isempty(geoTagIdx)
        plot3(geoTag(1,:), geoTag(2,:), zeros(size(geoTagIdx)), geoTagStyle{:});
    end
end

%restore "hold" status of the plot
if ~flagHold
    hold off
end

%set view
xlabel('x')
ylabel('y')
zlabel('z')

axis equal
view(3)
