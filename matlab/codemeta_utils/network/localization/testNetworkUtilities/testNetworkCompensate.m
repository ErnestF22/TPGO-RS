%function t_node=testNetworkCompensate(t_node,varargin)
%Rotate and scale estimates t_node().gi to compare them with
%t_node().gitruth.
%Optional arguments
%   'LeastSquares'      Find lambda, R and T that minimize
%                       norm(Ttruth-lambda*R*Testimated-T*ones(1,N),'fro'), where Ttruth
%                       and Testimated are the camera centers corresponding to
%                       the gi and gitruth fields of t_node. (This is the
%                       default)
%   'Camera', nCamera   Rotate and translate such that t_node(NCamera).gi
%                       coincides with t_node(NCamera).gitruth
%   'Mode', mode        Default: align translations
%       'Translations'  Align the translations
%       'RotationsOnly' Align the rotations only, without affecting the
%                       translations

%%AUTORIGHTS%%

function t_node=testNetworkCompensate(t_node,varargin)
structType=testNetworkDetectStructType(t_node);
N=testNetworkGetNumberOfNodes(t_node);

method='LeastSquares';
mode='Translations';
methodAbsolutePoses='reference';
flagRansac=false;
nCamera=1:N;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'leastsquares'
            method='LeastSquares';
        case 'camera'
            ivarargin=ivarargin+1;
            nCamera=varargin{ivarargin};
            if length(nCamera)==1
                method='Camera';
            end
        case 'mode'
            ivarargin=ivarargin+1;
            mode=varargin{ivarargin};
        case 'references'
            methodAbsolutePoses='reference';
        case 'poses'
            methodAbsolutePoses='pose';
        case 'methodabsoluteposes'
            ivarargin=ivarargin+1;
            methodAbsolutePoses=varargin{ivarargin};
        case 'ransac'
            flagRansac=true;
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end


switch lower(mode)
    case 'translations'
        switch lower(method)
            case 'leastsquares'
                Gtruth=cat(3,t_node.gitruth); %%!!
                Gestimated=cat(3,t_node.gi); %%!!

                if strcmp(methodAbsolutePoses,'pose')
                    for iNode=1:N
                        Gtruth(:,:,iNode)=invg(Gtruth(:,:,iNode));
                        Gestimated(:,:,iNode)=invg(Gestimated(:,:,iNode));
                    end
                end
                
                Ttruth=squeeze(Gtruth(1:3,4,nCamera)); %%!!
                Testimated=squeeze(Gestimated(1:3,4,nCamera)); %%!!
                
                if ~flagRansac
                    [~, ~, transf] = ... %%!!
                        procrustes(Ttruth', Testimated','reflection',false);
                else
                    NTrials=100;
                    NPointsSelected=4;
                    ransacThresh=0.01;
                    cMax=0;
                    for iTrial=1:NTrials
                        idx=randperm(N,NPointsSelected);
                        [~, ~, transf1] = procrustes(Ttruth(:,idx)', Testimated(:,idx)','reflection',false);
                        Tnew = (transf1.b * Testimated' * transf1.T + ones(N,1)*transf1.c(1,:))';
                        c1=sum(sum((Tnew-Ttruth).^2)<ransacThresh);
                        if c1>cMax
                            cMax=c1;
                            transf=transf1;
                        end
                    end
                end
                        
                    

                gGlobal=[transf.T' transf.c(1,:)'; zeros(1,3) 1]; %%!!

                Gnew=Gestimated; %%!!
                for iNode=1:N %%!!
                    Gnew(1:3,4,iNode)=Gnew(1:3,4,iNode)*transf.b; %%!!
                    Gnew(:,:,iNode)=gGlobal*Gnew(:,:,iNode); %%!!
                end
                
                if strcmp(methodAbsolutePoses,'pose')
                    Gnew=invg(Gnew);
                end
                
                switch structType
                    case 'single'
                        t_node.gi=Gnew;
                    case 'array'
                        for iNode=1:N
                            t_node(iNode).gi=Gnew(:,:,iNode);
                        end
                end
            case 'camera'
                %get reference camera poses
                switch structType
                    case 'single'
                        gn=t_node.gi(:,:,nCamera);
                        gntruth=t_node.gitruth(:,:,nCamera);
                    case 'array'
                        gn=t_node(nCamera).gi;
                        gntruth=t_node(nCamera).gitruth;
                end
                switch methodAbsolutePoses
                    case 'reference'
                        gGlobal=gntruth*invg(gn);
                        t_node=testNetworkTransform(t_node,gGlobal,'direction','left');
                    case 'pose'
                        gGlobal=invg(gn)*gntruth;
                        t_node=testNetworkTransform(t_node,gGlobal,'direction','right');
                    otherwise
                        error('Absolute pose method specification not valid')
                end
        end
    case 'rotationsonly'
        switch lower(method)
            case 'leastsquares'
                error('To be implemented')
            case 'camera'
                switch structType
                    case 'single'
                        Rn=t_node.gi(1:3,1:3,nCamera);
                        Rntruth=t_node.gitruth(1:3,1:3,nCamera);
                    case 'array'
                        gn=t_node(nCamera).gi(1:3,1:3);
                        gntruth=t_node(nCamera).gitruth(1:3,1:3);
                end
                
                switch methodAbsolutePoses
                    case 'reference'
                        gGlobal=[Rtruth*Rn' zeros(3,1); zeros(1,3) 1];
                        t_node=testNetworkTransform(t_node,gGlobal,'direction','left');
                    case 'pose'
                        gGlobal=[Rn'*Rtruth zeros(3,1); zeros(1,3) 1];
                        t_node=testNetworkTransform(t_node,gGlobal,'direction','right');
                    otherwise
                        error('Absolute pose method specification not valid')
                end
        end
    otherwise
        error(['Mode ' mode ' not valid']) 
end

