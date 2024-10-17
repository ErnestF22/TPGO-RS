%function t_node=testNetworkSymmetrizeEdges(t_node)
%Adds reverse edges using t_node.E. Adds the fields 'direction' and
%'idxRevE'. Also updates the fields 'NNodes', 'NEdges' and (if present) 'gij'.
function t_node=testNetworkSymmetrizeEdges(t_node)
structType=testNetworkDetectStructType(t_node);

switch structType
    case 'single'
        E=t_node.E;

        t_node.NNodes=max(E(:));
        t_node.NEdges=2*size(E,1);

        t_node.direction='directed';
        t_node.E=cat(1,E,fliplr(E));
        if isfield(t_node,'gij')
            t_node.gij=cat(3,t_node.gij,invg(t_node.gij));
        end
        t_node.A=sparse(t_node.E(:,1),t_node.E(:,2),ones(t_node.NEdges,1));
        %dispersion matrices
        if isfield(t_node,'dispersionMat')
            t_node.dispersionMat=repmat(t_node.dispersionMat,[1 1 2]);
            t_node.dispersionMatR=repmat(t_node.dispersionMatR,[1 1 2]);
            t_node.dispersionMatT=repmat(t_node.dispersionMatT,[1 1 2]);
        end
        %data for pose estimation
        if isfield(t_node,'xij')
            t_node.xij=[t_node.xij;t_node.xij];
            t_node.Xij=[t_node.Xij;t_node.Xij];
            t_node.flagReversePoseEstimate=[t_node.flagReversePoseEstimate; ~t_node.flagReversePoseEstimate];
        end
    case 'array'
        error('MATLAB:NotImplemented','Symmetrization of edges not implemented for array structure')
end
