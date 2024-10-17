function dispMat=generateDispersionMatrices(NEdges,method,var,dim)
if ~exist('dim','var')
    dim=3;
end

dispMat=zeros(dim,dim,NEdges);
dispersion=fliplr(1./var);
for iEdge=1:NEdges
    if var(1)==0
        dispMat(:,:,iEdge)=Inf(3);
    else
        switch lower(method)
            case 'anisotropic'
                D=diag(dispersion(1)+rand(dim,1)*(dispersion(2)-dispersion(1)));
                U=rot_randn(eye(dim));
                dispMat(:,:,iEdge)=U*D*U';
            case 'anisotropicdiagonal'
                dispMat(:,:,iEdge)=diag(dispersion(1)+rand(dim,1)*(dispersion(2)-dispersion(1)));
            case 'isotropic'
                dispMat(:,:,iEdge)=eye(dim)*(dispersion(1)+rand*(dispersion(2)-dispersion(1)));
            case 'identity'
                dispMat(:,:,iEdge)=eye(dim);
            case 'zero'
                dispMat(:,:,iEdge)=zeros(dim);
        end
    end
end
