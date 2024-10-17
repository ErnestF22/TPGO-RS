%Find a rotation R using at least 5 equations of the form n'*R*l=0
%function REst=rotation5pt(n,l)
function REst=rotation5pt(normals,lines)
d=4;
NPoints=size(normals,2);
if NPoints<9-d
    error('This algorithm requires at least %d normal-line correspondences',9-d)
end
MLin=zeros(NPoints,9);
for iPoint=1:NPoints
    MLin(iPoint,:)=kron(lines(:,iPoint),normals(:,iPoint));
end
[U,S,V]=svd(MLin);
A=reshape(V(:,end-3:end),3,3,d);

% D=nchoosek(2+d-1,2);
% 
% M1=zeros(6,D);
% M2=zeros(6,D);
% m12=zeros(6,1);
% cntlk=1;
% for k=1:3
%     for l=k:3
%         cntij=1;
%         for i=1:d
%             for j=i:d
%                 if i==j
%                     M1(cntlk,cntij)=A(:,k,i)'*A(:,l,j);
%                     M2(cntlk,cntij)=A(k,:,i)*A(l,:,j)';
%                 else
%                     M1(cntlk,cntij)=A(:,k,i)'*A(:,l,j)+A(:,k,j)'*A(:,l,i);
%                     M2(cntlk,cntij)=A(k,:,i)*A(l,:,j)'+A(k,:,j)*A(l,:,i)';
%                 end
%                 cntij=cntij+1;
%             end
%         end
%         if k==l
%             m12(cntlk)=1;
%         else
%             m12(cntlk)=0;
%         end        
%         cntlk=cntlk+1;
%     end
% end
% M=[M1;M2];
% m=[m12;m12];

[M,m]=buildSystemEmbeddedRotationQuadraticConstraints(A);
a=M\m;

%M gives 12 equations, but only 9 are linearly independent
% rankM=9;
% [UM,SM,VM]=svd(M);
% aEmb1=VM(:,1:rankM)*(SM(1:rankM,1:rankM)\UM(:,1:rankM)')*m;
% aEmb2=VM(:,rankM+1:end);
% 
% kerM=[aEmb1 aEmb2];
% 
% [x,lambda]=solveVeronese2LinearCombination(kerM,d);
% %scale solution to account the fact that lambda(1) should be one
% a=x/sqrt(abs(lambda(1)));


REst=combineMatrices(A,a(1:d));
REst=rot_proj(REst);

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

