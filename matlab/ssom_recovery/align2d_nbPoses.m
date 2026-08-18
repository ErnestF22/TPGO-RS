function Qx=align2d_nbPoses(v)
nbPoses=size(v,3);
Qx=repmat(zeros(size(v,1)),[1 1 nbPoses]);
for iPose=1:nbPoses
    Q=fliplr(orthComplement(v(:,:,iPose)));
    Qx(:,:,iPose)=flipud(orthCompleteBasis(Q)');
end
end %file function