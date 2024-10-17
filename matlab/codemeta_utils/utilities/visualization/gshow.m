%function gshow(E)
%Show the graph given by the edge list E
function [C]=gshow(E,varargin)
cyclestyles=char('b','r','g','c','m');
flagplot=false;

D=size(E,3);
if(D>1)
    for(d=1:D)
        subplot(ceil(D/6),min(D,6),d)
        [C]=gshow(E(:,:,d),varargin{:},'style',cyclestyles(mod(d-1,size(cyclestyles,1))+1,:));
    end
else
    N=[];
    flagcoordonly=false;    %compute only coordinates w/o showing anything
    prob=[];                %numerical values to show near the nodes
    C=[];                   %coordinates for the nodes
    style='b';              %line style
    cycles=[];              %cycles
    cols=0;                 %number of columns for the grid disposition

    %optional parameters
    ivarargin=1;
    while(ivarargin<=length(varargin))
        switch(varargin{ivarargin})
            case 'coordonly'
                flagcoordonly=true;
            case 'adj'  %E is, in fact, an adjacency matrix
                [I,J]=find(E>0);
                N=size(E,1);
                E=[I J];
            case 'prob'
                ivarargin=ivarargin+1;
                prob=varargin{ivarargin};
            case 'coords'
                ivarargin=ivarargin+1;
                C=varargin{ivarargin};
            case 'style'
                ivarargin=ivarargin+1;
                style=varargin{ivarargin};
            case 'cycles'
                ivarargin=ivarargin+1;
                cycles=varargin{ivarargin};
            case 'grid'
                ivarargin=ivarargin+1;
                cols=varargin{ivarargin};
            case 'simpleplot'
                flagplot=true;
            otherwise
                disp(varargin{ivarargin})
                error('Argument not valid!')
        end
        ivarargin=ivarargin+1;
    end


    if(isempty(C))
        if(~isempty(E))
            N=max(E(:));
        end
        if(cols<=0)
            xcoord=cos(2*pi*(1:N)'/N+1);
            ycoord=sin(2*pi*(1:N)'/N+1);
        else
            [x,y]=meshgrid(1:cols, ceil(N/cols):-1:1);
            xcoord=x(1:N)';
            ycoord=y(1:N)';
        end
        C=[xcoord ycoord];
    else
        N=size(C,1);
        xcoord=C(:,1);
        ycoord=C(:,2);
    end

    if(~flagcoordonly)
        if(isempty(cycles))
            if(flagplot)
                plot([xcoord(E(:,1))';xcoord(E(:,2))'],[ycoord(E(:,1))';ycoord(E(:,2))'],style)
                hold on
                plot(xcoord,ycoord,[style 'o'])
                hold off
            else
                quiver(xcoord(E(:,1)),ycoord(E(:,1)),xcoord(E(:,2))-xcoord(E(:,1)),ycoord(E(:,2))-ycoord(E(:,1)),0,style)
                for n=1:N
                    if(isempty(prob))
                        text(xcoord(n)+0.01,ycoord(n)+0.01,num2str(n));
                    else
                        text(xcoord(n)+0.01,ycoord(n)+0.01,num2str(prob(n)));
                    end
                end
            end            
            axis([min(xcoord)-0.5 max(xcoord)+0.5 min(ycoord)-0.5 max(ycoord)+0.5])
%            axis equal
        else
            L=size(cycles,2);
            f=ishold;
            for l=1:L
                gshow_cycle(E,C,cycles(:,l),deblank(cyclestyles(mod(l-1,size(cyclestyles,1))+1,:)))
                hold on
            end
            if(~f)
                hold off
            end

        end
    end
    axis('equal')
end
