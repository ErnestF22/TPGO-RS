function fastMult3x3MatMat_test
N=500000;
A=randn(3,3,N);
B=randn(3,3,N);

disp('Naive implementation')
tic
C1=zeros(3,3,N);
for iN=1:N
    C1(:,:,iN)=A(:,:,iN)*B(:,:,iN);
end
toc

disp('Vectorized implementation')
tic
C2=fastMult3x3MatMat(A,B);
toc

disp('C implementation')
tic
C3=mexFastMult3x3MatMat(A(:),B(:));
toc

disp('Max of abs difference')
disp(max(abs(C1(:)-C2(:))))
disp(max(abs(C1(:)-C3(:))))
