function POCOptimalWLaplacian
load POCBlockJacobi_inputData %this gives a sample H for a 2-d problem

L=makeL(H);
sL=eig(L);
sLThresholded=sL(sL>1e-4);
wOpt=(max(sLThresholded)+min(sLThresholded))/2;
X=makeX(H,wOpt);
disp(sort(eig(X)))


function L=makeL(H)
D=blkZeroOffDiag(H,2);
Dhalf=sqrtPSDMatrix(D);
L=(Dhalf\H)/Dhalf;

function X=makeX(H,w)
D=w*blkZeroOffDiag(H,2);
B=D-H;
Dhalf=sqrtPSDMatrix(D);

X=(Dhalf\B)/Dhalf;
