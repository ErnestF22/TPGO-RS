%2-D or 3-D graph drawing function
%
%function graphShow(E,x,R,varargin)
%Optional arguments
%   'Normals',normals       vectors to show at each node
%   'Segmentation',members  color edges according to given segmentation
%
function graphShow(E,x,R,varargin)
resetRands()
flagDrawNormals=false;
flagDrawLocalAxes=false;
flagSegmentation=false;
segmentsColorMap=[
    1 0 0;
    0 0.8 0;
    0 0 1;
    1 1 0;
    0 1 1;
    1 0 1;
    0 0 0;
    rand(5,3);
    ];

styleNormals={'m'};
markerOpts={'o','MarkerSize',20};
lineOpts={'b-','LineWidth',1.5};
textOpts={'HorizontalAlignment','Center','FontWeight','Bold',...
    'FontSize',15,'BackgroundColor',[1 1 1]};

if exist('R','var') && ~isempty(R)
    flagDrawLocalAxes=true;
end

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'normals'
            ivarargin=ivarargin+1;
            n=varargin{ivarargin};
            flagDrawNormals=true;
        case 'segmentation'
            ivarargin=ivarargin+1;
            members=varargin{ivarargin};
            flagSegmentation=true;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

[dimData,NNodes]=size(x);

% if flagSegmentation && isempty(members)
%     members=ones(1,NNodes);
% end

switch dimData
    case 2
        if ~flagSegmentation
            plot([x(1,E(:,1)); x(1,E(:,2))],[x(2,E(:,1)); x(2,E(:,2))],lineOpts{:})
            hold on
        else
            if ~isempty(members)
                for iSegment=1:size(members,1)
                    group=find(members(iSegment,:));
                    ESegment=E(ismember(E(:,1),group) & ismember(E(:,2),group),:);
                    plot([x(1,ESegment(:,1)); x(1,ESegment(:,2))],[x(2,ESegment(:,1)); x(2,ESegment(:,2))],lineOpts{:},'Color',segmentsColorMap(iSegment,:))
                    hold on
                end
            end
        end
        plot(x(1,:),x(2,:),markerOpts{:})
        for iNode=1:NNodes
            text(x(1,iNode),x(2,iNode),num2str(iNode),textOpts{:})
            if flagDrawLocalAxes
                plotframe(x(:,iNode),R(:,:,iNode),'scale',0.3)
            end
        end
        if flagDrawNormals
            quiver(x(1,:),x(2,:),n(1,:),n(2,:),styleNormals{:})
        end
    case 3
        if ~flagSegmentation
            plot3([x(1,E(:,1)); x(1,E(:,2))],[x(2,E(:,1)); x(2,E(:,2))],[x(3,E(:,1)); x(3,E(:,2))],lineOpts{:})
            hold on
        else
            if ~isempty(members)
                for iSegment=1:size(members,1)
                    group=find(members(iSegment,:));
                    ESegment=E(ismember(E(:,1),group) & ismember(E(:,2),group),:);
                    plot3([x(1,ESegment(:,1)); x(1,ESegment(:,2))],[x(2,ESegment(:,1)); x(2,ESegment(:,2))],[x(3,ESegment(:,1)); x(3,ESegment(:,2))],lineOpts{:},'Color',segmentsColorMap(iSegment,:))
                    hold on
                end
            end
        end
        plot3(x(1,:),x(2,:),x(3,:),markerOpts{:});
        for iNode=1:NNodes
            text(x(1,iNode),x(2,iNode),x(3,iNode),num2str(iNode),textOpts{:})
            if flagDrawLocalAxes
                plotframe(x(:,iNode),R(:,:,iNode),'scale',0.3)
            end
        end
        if flagDrawNormals
            quiver3(x(1,:),x(2,:),x(3,:),n(1,:),n(2,:),n(3,:),styleNormals{:})
        end
end
hold off
