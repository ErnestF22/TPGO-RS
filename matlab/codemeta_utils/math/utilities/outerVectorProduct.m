%Compute the outer product of tensors
%function T=outerVectorProduct(v1,v2,...,vN)
%Compute the matrix T of dimension [length(v1) x ... x length(vN)] where
%each entry is computed as T(i1,...,iN)=v1(i1)*...*vN(iN)
function T=outerVectorProduct(varargin)
vArray=varargin;
Nv=length(vArray);
dim=cellfun(@length,vArray);
T=zeros(dim);
NT=numel(T);
h=getSeqComb(dim);
for iT=1:NT
    idx=h();
    t=1;
    for iv=1:Nv
        t=t*vArray{iv}(idx(iv));
    end
    idxCell=num2cell(idx);
    T(idxCell{:})=t;
end
