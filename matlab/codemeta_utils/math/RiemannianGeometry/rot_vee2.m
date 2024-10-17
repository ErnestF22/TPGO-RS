function v=rot_vee2(R,A)
T=rot_tangentBasis(R);
d=size(T,3);
N=size(A,3);
v=zeros(d,N);
for iN=1:N
    for iT=1:d
        v(iT,iN)=rot_metric(R,A(:,:,iN),T(:,:,iT));
    end
end
