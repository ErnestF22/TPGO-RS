function tnodeshow(t_node,varargin)
    flagcompensate=false;
    flagdisplaytruth=false;
    flagdisplaystructure=false;
    
    %optional parameters
    ivarargin=1;
    while(ivarargin<=length(varargin))
        if(ischar(varargin{ivarargin}))
            switch(varargin{ivarargin})
                case 'compensate'
                    flagcompensate=true;
                    ivarargin=ivarargin+1;
                    t_node_ref=varargin{ivarargin};
                case 'displaytruth'
                    flagdisplaytruth=true;
                case 'structure'
                    ivarargin=ivarargin+1;
                    X=varargin{ivarargin};
                    flagdisplaystructure=true;
                otherwise
                    if(~flagsilentparse)
                        error(['Argument ' varargin{ivarargin} ' not valid!'])
                    else
                        warning(['Argument ' varargin{ivarargin} ' not valid!'])
                    end                
            end
        else
            if(~flagsilentparse)
                error(['Argument ' varargin{ivarargin} ' not valid!'])
            else
                warning(['Argument ' varargin{ivarargin} ' not valid!'])
            end
        end
        ivarargin=ivarargin+1;
    end

    N=length(t_node);

    A=zeros(N);
    for(inode=1:N)
        A(inode,:)=t_node(inode).aij;
    end
    %find edge list
    [E(:,1),E(:,2)]=find(A~=0);
    
    t_node_comp=t_node;
    if(flagcompensate)
            t_node_comp=compensate(t_node_comp,t_node_ref,'rotationsRigid');
%            t_node_comp=compensate(t_node_comp,t_node_ref,'rotations');
        t_node_comp=compensate(t_node_comp,t_node_ref,'translations');
        t_node_comp=compensate(t_node_comp,t_node_ref,'scale');
    end
    holdstate=ishold();
    
    for(inode=1:N)
        draw3dcamera(t_node_comp(inode).gi(1:3,1:3),t_node_comp(inode).gi(1:3,4))
        hold on
        if(inode==1)
            ti=t_node_comp(inode).gi(1:3,4);
            text(ti(1),ti(2),ti(3)+1.1,'Camera 1','HorizontalAlignment','Right');
        end
    end
    for(iedge=1:size(E,1))
        l(1,iedge,:)=t_node_comp(E(iedge,1)).gi(1:3,4);
        l(2,iedge,:)=t_node_comp(E(iedge,2)).gi(1:3,4);
        ltruth(1,iedge,:)=t_node_comp(E(iedge,1)).gitruth(1:3,4);
        ltruth(2,iedge,:)=t_node_comp(E(iedge,2)).gitruth(1:3,4);
    end
    plot3(l(:,:,1),l(:,:,2),l(:,:,3),'b','LineWidth',1)
    if(flagdisplaytruth)
        plot3(ltruth(:,:,1),ltruth(:,:,2),ltruth(:,:,3),'r','LineWidth',1)
    end
    if(flagdisplaystructure)
        plot3(X(1,:),X(2,:),X(3,:),'k.')
    end
    hold off

    axis equal
    view(3)
    set(gcf,'Renderer','zbuffer')
    if(holdstate)
        hold on
    else
        hold off
    end
end

function [t_node]=compensate(t_node,t_node_ref,type)
    N=length(t_node);
    switch type
        case 'rotations'
            R1=t_node(1).gi(1:3,1:3);
            R1ref=t_node_ref(1).gi(1:3,1:3);
            for(inode=1:N)
                t_node(inode).gi(1:3,1:3)=R1ref*R1'*t_node(inode).gi(1:3,1:3);
            end
        case 'rotationsRigid'
            R1=t_node(1).gi(1:3,1:3);
            R1ref=t_node_ref(1).gi(1:3,1:3);
            for(inode=1:N)
                t_node(inode).gi=[R1ref*R1' zeros(3,1);zeros(1,4)]*t_node(inode).gi;
            end
        case 'translations'
            t_node=removeTranslationMean(t_node);
            muref=computeTranslationMean(t_node_ref);
            for(inode=1:N)
                t_node(inode).gi(1:3,4)=t_node(inode).gi(1:3,4)+muref;
            end
        case 'scale'
            mu=computeTranslationMean(t_node);
            muref=computeTranslationMean(t_node_ref);
            t_node=removeTranslationMean(t_node,mu);
            t_node_ref=removeTranslationMean(t_node_ref,mu);
            
            scale1=computeTranslationDistanceMean(t_node);
            scaleref=computeTranslationDistanceMean(t_node_ref);
            
            for(inode=1:N)
                t_node(inode).gi(1:3,4)=t_node(inode).gi(1:3,4)/scale1*scaleref+mu;
            end
             
    end
end

function t_node=removeTranslationMean(t_node,mu)
    N=length(t_node);
    if(nargin<2)
        mu=computeTranslationMean(t_node);
    end
    for(inode=1:N)
        t_node(inode).gi(1:3,4)=t_node(inode).gi(1:3,4)-mu;
    end
end

function mu=computeTranslationMean(t_node)
    N=length(t_node);
    mu=0;
    for(inode=1:N)
        mu=mu+t_node(inode).gi(1:3,4);
    end
    mu=mu/N;
end

function scale=computeTranslationDistanceMean(t_node)
    N=length(t_node);
    scale=0;
    for(inode=1:N)
        scale=scale+norm(t_node(inode).gi(1:3,4));
    end
    scale=scale/N;
end
