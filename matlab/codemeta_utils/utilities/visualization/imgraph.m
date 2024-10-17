function imgraph(t_node,varargin)
N=length(t_node);

imgW=0.1;
imgH=0.1;
imgR=0.8;
maxLineWidth=10;

normMode='all';

displaySet=[1:N];
%flagsequence=false;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(varargin{ivarargin})
        case 'norm'
            ivarargin=ivarargin+1;
            normMode=varargin{ivarargin};
        case 'display'
            ivarargin=ivarargin+1;
            displaySet=varargin{ivarargin};
%         case 'sequence'
%             flagsequence=true;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

axes('position',[0  0  1  1],'XTick',[],'YTick',[])
axis([0 1 0 1])

for(inode=1:N)
    x(inode)=(imgR*cos(2*pi/N*(inode-1))+1)/2;
    y(inode)=(imgR*sin(2*pi/N*(inode-1))+1)/2;
end

if(isfield(t_node(1),'sim'))
    sim=reshape([t_node.sim],N,N);
    switch(normMode)
        case 'all'
            maxSim=max(sim(:))*ones(N,N);
        case 'single'
            maxSim=max(sim,[],2)*ones(1,N);
        case 'sqrt'
            Nfeats=sqrt([t_node.Nfeats]);
            maxSim=Nfeats'*Nfeats;
        case 'nfeats'
            Nfeats=sqrt([t_node.Nfeats]);
            maxSim=Nfeats'*ones(1,N);
        otherwise
            error('Normalization method not valid!')
    end
    for(inode=displaySet)
        for(jnode=1:N)
            if(inode~=jnode && t_node(inode).sim(jnode)>0)
                relSim=t_node(inode).sim(jnode)/maxSim(inode,jnode);
                lineColor=(1-relSim)*[1 1 0]+relSim*[0 0 1];
                line([x(inode) x(jnode)],[y(inode) y(jnode)],'LineWidth',maxLineWidth*relSim,'Color',lineColor);
            end
        end
    end
end
for(inode=1:N)
    axes('position',[x(inode)-0.5*imgW y(inode)-0.5*imgH imgW imgH])
    imshow(t_node(inode).img,[])
end
