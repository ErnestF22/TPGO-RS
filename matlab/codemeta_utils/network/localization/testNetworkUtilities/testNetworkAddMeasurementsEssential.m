function t_node=testNetworkAddMeasurementsEssential(t_node,varargin)
method='truth';
sigmaNoiseE=0.001;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'method'
            ivarargin=ivarargin+1;
            method=varargin{ivarargin};
        case 'sigmanoise'
            ivarargin=ivarargin+1;
            sigmaNoiseE=varargin{ivarargin};
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

structType=testNetworkDetectStructType(t_node);

E=testNetworkGetEdges(t_node);
N=testNetworkGetNumberOfNodes(t_node);

if ~isfield(t_node,'Eij')
    [t_node.Eij]=deal(repmat(zeros(3),[1,1,N]));
end
switch lower(method)
    case 'truth'
        [t_node.Eij]=deal(t_node.Eijtruth);
        [t_node.QEij]=deal(t_node.QEijtruth);%essential_fromG(t_node.gitruth(:,:,E(:,1)),t_node.gitruth(:,:,E(:,2)));
    case 'noisy'
        s=sigmaNoiseE;
        switch structType
            case 'single'
                if isempty(s)
                    s=t_node.dispersionMatE;
                end
                t_node.QEij=essential_randn(t_node.QEijtruth,s);
                t_node.Eij=essential_getE(t_node.QEij);
            case 'array'
                for iNode=1:N
                    if isempty(sR)
                        sR=t_node(iNode).dispersionMatR;
                    end
                    if isempty(sT)
                        sT=t_node(iNode).dispersionMatT;
                    end
                    t_node(iNode).gij=noiserigid(t_node(iNode).gijtruth,sR,sT);
                end
        end
end


