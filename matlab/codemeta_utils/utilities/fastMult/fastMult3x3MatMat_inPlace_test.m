function fastMult3x3MatMat_inPlace_test
N=100000;
A=randn(3,3,N);
B=randn(3,3,N);

disp('C implementation')
tic
for iN=1:N
    C1=mexfastMult3x3MatMat(A(:,:,iN),B(:,:,iN));
end
toc

disp('C inplace implementation')
tic
C2=zeros(3,3);
for iN=1:N
    mexfastMult3x3MatMat_inPlace(A(:,:,iN),B(:,:,iN),C2);
end
toc

disp('Max of abs difference')
disp(max(abs(C1(:)-C2(:))))
