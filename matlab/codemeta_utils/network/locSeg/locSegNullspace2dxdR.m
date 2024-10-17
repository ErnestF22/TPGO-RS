function [dx,dR]=locSegNullspace2dxdR(bN,R,dimData,dimCoordinates)
NNodes=size(bN,1)/dimCoordinates;
idxTransl=reshape((1:dimData)'*ones(1,NNodes)+ones(dimData,1)*((0:NNodes-1)*dimCoordinates),[],1);
dx=reshape(bN(idxTransl,:),dimData,NNodes,[]);
idxRot=reshape((dimData+1:dimCoordinates)'*ones(1,NNodes)+ones(dimCoordinates-dimData,1)*((0:NNodes-1)*dimCoordinates),[],1);
dR=reshape(bN(idxRot,:),dimCoordinates-dimData,NNodes,[]);
if dimData==3
    for iMode=1:size(dR,3)
        for iNode=1:NNodes
            dR(:,iNode,iMode)=R(:,:,iNode)*dR(:,iNode,iMode);
        end
    end
end