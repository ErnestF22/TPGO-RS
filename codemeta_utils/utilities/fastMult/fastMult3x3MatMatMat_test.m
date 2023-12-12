function fastMult3x3MatMatMat_test
N=100000;
A=randn(3,3,N);
B=randn(3,3,N);
C=randn(3,3,N);

disp('Naive implementation')
tic
D1=zeros(3,3,N);
for iN=1:N
    D1(:,:,iN)=A(:,:,iN)*B(:,:,iN)*C(:,:,iN);
end
toc

disp('Vectorized implementation')
tic
D2=fastMult3x3MatMatMat(A,B,C);
toc

disp('C implementation, two steps')
tic
D3=mexFastMult3x3MatMat(mexFastMult3x3MatMat(A,B),C);
toc

disp('C implementation')
tic
D4=mexFastMult3x3MatMatMat(A,B,C);
toc

disp('Max of abs difference')
disp(max(abs(D1(:)-D2(:))))
disp(max(abs(D1(:)-D3(:))))
disp(max(abs(D1(:)-D4(:))))
