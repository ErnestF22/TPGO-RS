function d=rot_distn(R1,R2)
A12=rot_log(R1,R2);
N=size(R2,3);
d=zeros(N,1);
for iN=1:N
    d(iN)=sqrt(rot_metric(R1,A12(:,:,iN),A12(:,:,iN)));
end
