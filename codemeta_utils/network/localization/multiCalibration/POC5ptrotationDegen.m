function POC5ptrotationDegen
%%
load rotation5pt_testData
d=4;
normals=n;
lines=l;
NPoints=size(normals,2);
if NPoints<9-d
    error('This algorithm requires at least %d normal-line correspondences',9-d)
end
MLin=zeros(NPoints,6);
for iPoint=1:NPoints
    MLin(iPoint,:)=kron(lines(2:3,iPoint),normals(:,iPoint));
end

[U,S,V]=svd(MLin);
nullMLin=V(:,end-d+4:end);
nullityMLin=size(nullMLin,2);

A=cat(3,reshape([eye(3); zeros(6,3)],3,3,3),...
    [zeros(3,1,nullityMLin) reshape(nullMLin,3,2,[])]);

for i=1:4
    disp(diag(normals'*A(:,:,i)*lines))
end

disp(subspace(R(:),reshape(A,9,[])))

aTruth=reshape(A,9,[])\R(:);
aTruthEmb=veronese(aTruth,2);
disp(norm(R-combineMatrices(A,aTruth),'fro'))

[M,m]=buildSystemEmbeddedRotationQuadraticConstraints(A);

%M gives 12 equations, but only 9 are linearly independent
rankM=sum(svd(M)>1e-6);
[UM,SM,VM]=svd(M);
aEmb1=VM(:,1:rankM)*(SM(1:rankM,1:rankM)\UM(:,1:rankM)')*m;
aEmb2=VM(:,rankM+1:end);

kerM=[aEmb1 aEmb2];

%lambda=lambda/lambda(1);
lambdaTruth=kerM\aTruthEmb;

[x,lambda,K,lambdabar]=solveVeronese2LinearCombination(kerM,d);
%scale solution to account the fact that lambda(1) should be one
a=x/sqrt(abs(lambda(1)));
aEmb=kerM*lambda/lambda(1);

disp('[lambda/lambda(1) lambdaTruth]')
disp([lambda/lambda(1) lambdaTruth])
disp('[aTruthEmb aEmb]')
disp([aTruthEmb aEmb])

REst=combineMatrices(A,a);
REst=REst*det(REst);
disp([R REst R-REst])

function B=combineMatrices(A,a)
[d1,d2,n]=size(A);
B=reshape(reshape(A,[],n)*a,d1,d2);

function [M,m]=buildSystemEmbeddedRotationQuadraticConstraints(A)
d=size(A,3);
D=nchoosek(2+d-1,2);

M1=zeros(6,D);
M2=zeros(6,D);
m12=zeros(6,1);
cntlk=1;
for k=1:3
    for l=k:3
        cntij=1;
        for i=1:d
            for j=i:d
                if i==j
                    M1(cntlk,cntij)=A(:,k,i)'*A(:,l,j);
                    M2(cntlk,cntij)=A(k,:,i)*A(l,:,j)';
                else
                    M1(cntlk,cntij)=A(:,k,i)'*A(:,l,j)+A(:,k,j)'*A(:,l,i);
                    M2(cntlk,cntij)=A(k,:,i)*A(l,:,j)'+A(k,:,j)*A(l,:,i)';
                end
                cntij=cntij+1;
            end
        end
        if k==l
            m12(cntlk)=1;
        else
            m12(cntlk)=0;
        end        
        cntlk=cntlk+1;
    end
end
M=[M1;M2];
m=[m12;m12];
