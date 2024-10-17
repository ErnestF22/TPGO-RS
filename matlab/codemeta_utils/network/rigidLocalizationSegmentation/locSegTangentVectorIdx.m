function [idxTranslation,idxRotation,idxCoordinates]=locSegTangentVectorIdx(dimData,NNodes)
[dimCoordinates,dimRotations]=locSegDimData2DimCoordinates(dimData);
uNodes=ones(NNodes,1);
idxBase=(0:NNodes-1)*dimCoordinates;
idxRotation=(dimData+1:dimCoordinates)'*uNodes'+ones(dimRotations,1)*idxBase;
idxTranslation=(1:dimData)'*uNodes'+ones(dimData,1)*idxBase;
idxCoordinates=[idxTranslation;idxRotation];
