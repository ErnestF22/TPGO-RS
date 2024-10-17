function Ri=sfm_rawRotationsDevectorize(RiVec)
NRotations=size(RiVec,1)/3;
Ri=reshape(RiVec',[3 3 NRotations]);
