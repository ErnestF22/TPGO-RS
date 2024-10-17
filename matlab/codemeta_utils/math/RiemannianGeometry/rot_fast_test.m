function rot_fast_test
N=10000;
v=randn(3,N);

disp('Naive Implementation')
tic
R1=zeros(3,3,N);
for iNode=1:N
    R1(:,:,iNode)=rot(v(:,iNode));
end
toc

disp('Vectorized implementation')
tic
[vnorm,nv]=cnormalize(v);
R2=rot_fast(vnorm,nv);
toc

disp('Max of abs difference')
disp(max(abs(R1(:)-R2(:))))
