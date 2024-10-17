function POCReviewPatternFormation
%
%Check of CDC2018 submission, eq. (7), (8).
%%
resetRands()
n=4;
gp=1/2;
kp=1/3;
kv=1/4;

A=rand(n,n)>0.5;

%%
for iA=1:n
    A(iA,iA)=0;
    %the following ensures that we have a spanning tree
    if iA~=n
        A(iA,iA+1)=1;
    end
end

%%

L=adj2laplacianmatrix(A);
B=[zeros(n,n) eye(n,n); -gp*eye(n,n)-kp*L -kv*L];

%%
[ti,Dmui]=eig(L);
mui=-diag(Dmui);
%%

lambdan=(kv*mui+sqrt(kv^2*mui.^2+4*(kp*mui-gp)))/2;
lambdap=(kv*mui-sqrt(kv^2*mui.^2+4*(kp*mui-gp)))/2;

[qsi,Dlambdai]=eig(B);
lambdai=diag(Dlambdai);

[sort([lambdan;lambdap]) sort(lambdai)]

%check that evects of B have the structure [t;lambda t]
disp(qsi(n+1:2*n,:)-qsi(1:n,:)*Dlambdai)

%%

lambdaiRoots=[];
for id=1:n
    lambdaiRoots=[lambdaiRoots; roots([1 -kv*mui(id) -(-gp+kp*mui(id))])];
end

disp([sort(lambdaiRoots) sort(lambdai)])

