function fastMult3x3MatHat_test
N=100000;
A=randn(3,3,N);
v=randn(3,N);

disp('Naive implementation')
tic
B1=zeros(3,3,N);
v1=permute(v,[1 3 2]);
z=zeros(1,1,N);
hatv=[z v1(3,1,:) -v1(2,1,:); -v1(3,1,:) z v1(1,1,:); v1(2,1,:) -v1(1,1,:) z];
for iN=1:N
    B1(:,:,iN)=A(:,:,iN)*hatv(:,:,iN);
end
toc

disp('Vectorized implementation')
tic
B2=fastMult3x3MatHat(A,v);
toc

disp('C implementation')
tic
B3=mexFastMult3x3MatHat(A,v);
toc

disp('Max of abs difference')
disp(max(abs(B1(:)-B2(:))))
disp(max(abs(B1(:)-B3(:))))
