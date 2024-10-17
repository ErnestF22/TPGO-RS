%function EA=rot3_expDiffSeries(A)
%Computes (in closed form) a series used by rot3_expDiff
function EA=rot3_expDiffSeries(A)
A=A';
[UA,SA,VA]=svd(A);
% A2=UA'*A*UA;
% A2=A2(1:2,1:2);

a=sign(UA(:,1)'*VA(:,2))*SA(1,1);
if a==0
    invA2expA2minusIdentity=eye(2);
else
    invA2expA2minusIdentity=[sin(a) cos(a)-1; -cos(a)+1 sin(a)]/a;
end

EA=UA*blkdiag(invA2expA2minusIdentity,1)*UA';
