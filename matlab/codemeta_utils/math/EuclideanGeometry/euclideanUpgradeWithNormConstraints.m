function [invS,ARectified]=euclideanUpgradeWithNormConstraints(A,desiredNorm)
flagAHomogeneousCoords=false;

%if A is in homogeneous coordinates, de-homogenize it
if size(A,1)==4
    flagAHomogeneousCoords=true;
    A=A(1:3,:)./(ones(3,1)*A(4,:));
end

%First we need to find Q=invS'*invS such that [a;1]'*Q*[a;1]=desiredNorm^2
%We vectorize the equation, and build the corresponding linear system
NA=size(A,2);
B=zeros(NA,9);
for iA=1:NA
    a=A(:,iA);
    B(iA,:)=[a(1)^2, 2*a(1)*a(2), 2*a(1)*a(3), 2*a(1)*1, a(2)^2, 2*a(2)*a(3), 2*a(2)*1, a(3)^2, 2*a(3)*1];
end
b=desiredNorm^2;
u=b*ones(NA,1);

%find particular solution for Q with Q(end,end)=0
qEst=B\u;
QEst1=[qEst(1) qEst(2) qEst(3) qEst(4); qEst(2) qEst(5) qEst(6) qEst(7); qEst(3) qEst(6) qEst(8) qEst(9); qEst(4) qEst(7) qEst(9) 0];

%find combination of the solution above with a solution in the nullspace so
%that the linear combination is rank deficient (zero determinant)
QEst2=diag([0 0 0 b]);
p=determinantLinearCombinationCoeff(QEst1,b);
a=-p(2)/p(1);
QEst=a*QEst1+(1-a)*QEst2;

%factorize the solution to get (a solution for) invS
[~,S,V]=svd(QEst);
invS=diag(sqrt(diag(S(1:3,1:3))))*V(:,1:3)';
invS=[invS; 0 0 0 1];

%apply transformation to original data, if requested
if nargout>1
    ARectified=invS*[A;ones(1,NA)];
    if ~flagAHomogeneousCoords
        ARectified=ARectified(1:3,:);
    end
end

end

function p=determinantLinearCombinationCoeff(A,b)
% compute the coefficients of the polynomial in x given by the 
% determinat of x*A+(1-x)*E*b, where E=diag([0....0 1]), from the
% determinant of A, x and b are scalars
N=size(A,1);
c1=det(A);
c2=b*det(A(1:N-1,1:N-1));
p1=c1-c2;
p2=c2;
%d=p1*a^N+p2*a^(N-1);
p=[p1 p2 zeros(1,N-1)];
end
