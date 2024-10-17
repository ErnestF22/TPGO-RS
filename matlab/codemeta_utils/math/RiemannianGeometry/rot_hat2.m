function A=rot_hat2(R,v)
NR=size(R,3);
N=size(v,2);

if NR~=1 && NR~=N
    error('MATLAB:argumentCheck','Number of rotations and number of tangent vectors should be consistent')
end

T=rot_tangentBasis(R(:,:,1));
d=size(T,3);
A=zeros([size(R,1) size(R,2) N]);
for iN=1:N
    if NR>1
        T=rot_tangentBasis(R(:,:,iN));
    end
        
    for iT=1:d
        A(:,:,iN)=A(:,:,iN)+v(iT,iN)*T(:,:,iT);
    end
end
