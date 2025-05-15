function Qalign=align3d(v)
vFlat=reshape(v,size(v,1),[]);
[U,~,~]=svd(vFlat);
Qalign=fliplr(orthCompleteBasis(U(:,4:end)))';
end %file function