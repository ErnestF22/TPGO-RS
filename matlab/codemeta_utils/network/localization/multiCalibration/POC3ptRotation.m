function POC3ptRotation
R=rot_randn();
a=[2;3;4;1];
d=length(a);

A=zeros(3,3,d);
A(:,:,d)=R;
for i=1:d-1
    A(:,:,i)=randn(3);
    A(:,:,d)=A(:,:,d)-a(i)*A(:,:,i);
end

RRec=zeros(3);
for i=1:d
    RRec=RRec+A(:,:,i)*a(i);
end
disp(RRec'*RRec)


aEmb=veronese(a,2);
D=length(aEmb);

M1=zeros(6,D);
M2=zeros(6,D);
m=zeros(6,1);
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
            m(cntlk)=1;
        else
            m(cntlk)=0;
        end        
        cntlk=cntlk+1;
    end
end
M=[M1;M2];
m=[m;m];
disp(M*aEmb-m)
disp(size(M,2)-rank(M))
%disp([aEmb [M1;M2]\[m;m]])
