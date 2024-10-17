function fastMult3x3MatIdxTranspMatMatIdx_test
N=100;
NE=10000;
A=randn(3,3,N);
B=randn(3,3,NE);
C=randn(3,3,N);

eA=randint(1,NE,[1 N]);
eC=randint(1,NE,[1 N]);

disp('Naive implementation')
tic
D1=zeros(3,3,N);
for iNE=1:NE
    D1(:,:,iNE)=A(:,:,eA(iNE))'*B(:,:,iNE)*C(:,:,eC(iNE));
end
toc

disp('Matrix-expansion vectorized implementation')
tic
AA=permute(A(:,:,eA),[2 1 3]);
CC=C(:,:,eC);
D2=fastMult3x3MatMatMat(AA,B,CC);
toc

disp('Vectorized implementation')
tic
D3=fastMult3x3MatIdxTranspMatMatIdx(A,eA,B,C,eC);
toc

disp('C implementation')
tic
D4=mexFastMult3x3MatIdxTranspMatMatIdx(A,eA',B,C,eC');
toc

disp('Max of abs difference')
disp(max(abs(D1(:)-D2(:))))
disp(max(abs(D1(:)-D3(:))))
disp(max(abs(D1(:)-D4(:))))
