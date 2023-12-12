%Compute determinant of a matrix using Leibniz formula
%function d=detLeibniz(A)
%This function is much slower and less precise than det(), but it has a
%minimum memory footprint
function d=detLeibniz(A)
N=size(A,1);
d=0;
h=getSeqPerm(1:N,'method','sjc');
b=1:N:N*N;
for iN=1:factorial(N)
    [sigma,sigmaSign]=h();
    d=d+sigmaSign*prod(A(b+sigma-1));
end
