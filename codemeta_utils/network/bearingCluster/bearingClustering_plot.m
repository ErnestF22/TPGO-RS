function bearingClustering_plot(x,E,membership,varargin)
NEdges=size(E,1);
flagVertexNumber=false;
%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'vertexnumber'
            flagVertexNumber=true;
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
    
end

if ~exist('membership','var') || isempty(membership)
    membership=ones(1,NEdges);
end

c=unique(membership);
K=length(c);
cmap=jet(K);

plotPoints(x)
hold on

for iEdge=1:NEdges    
    plotPoints(x(:,E(iEdge,:)),'-','Color',cmap(membership(iEdge)==c,:))
end

if flagVertexNumber
    dim=size(x,1);
    NNodes=size(x,2);
    switch dim
        case 2
            for ix=1:NNodes
                text(x(1,ix),x(2,ix),num2str(ix))
            end
        case 3
            for ix=1:NNodes
                text(x(1,ix),x(2,ix),x(3,ix),num2str(ix))
            end
    end
end

hold off
grid on
axis equal

% function bearingClustering_plot(x,E,membership)
% NEdges=size(E,1);
% 
% if ~exist('membership','var') || isempty(membership)
%     membership=ones(1,NEdges);
% end
% 
% c=unique(membership);
% K=length(c);
% cmap=jet(K);
% 
% plotPoints(x)
% hold on
% for iEdge=1:NEdges    
%     plotPoints(x(:,E(iEdge,:)),'-','Color',cmap(membership(iEdge)==c,:))
% end
% hold off
% grid on
% axis equal
