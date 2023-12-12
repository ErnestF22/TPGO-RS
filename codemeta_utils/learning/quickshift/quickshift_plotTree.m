%Plot a quickshift tree
function quickshift_plotTree(X,treeVectorMembership,varargin)
flagArrows=true;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    if isstring(varargin{ivarargin})
        switch lower(varargin{ivarargin})
            case 'noarrows'
                flagArrows=false;
                varargin(ivarargin)=[];
                ivarargin=ivarargin-1;
            otherwise
        end
    end
    ivarargin=ivarargin+1;
end


XTo=X(:,treeVectorMembership);
%XDiff=XTo-X;
switch size(X,1)
    case 2
        if flagArrows
            quiver(X(1,:),X(2,:),XTo(1,:)-X(1,:),XTo(2,:)-X(2,:),0,'k-',varargin{:})
        else
            plot([X(1,:);XTo(1,:)],[X(2,:); XTo(2,:)],'k-',varargin{:})
        end
    case 3
        if flagArrows
            quiver3(X(1,:),X(2,:),X(3,:),XTo(1,:)-X(1,:),XTo(2,:)-X(2,:),XTo(3,:)-X(3,:),0,'k-',varargin{:})
        else
            plot3([X(1,:);XTo(1,:)],[X(2,:); XTo(2,:)],[X(3,:); XTo(3,:)],'k-',varargin{:})
        end
end
