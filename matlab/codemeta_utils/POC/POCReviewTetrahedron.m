function POCReviewTetrahedron
%%
syms p1 s1 p2 s2 p3 s3 p4 s4 f df
F=sym('F',[4 4]);
p=[p1;p2;p3;p4];
s=[s1;s2;s3;s4];
z=[p1;s1;p2;s2;p3;s3;p4;s4];

tij=@(i,j) acos(cos(p(i))*cos(p(j))*cos(s(i)-s(j))+sin(p(i))*sin(p(j)));
dpij=@(i,j) 1/cos(p(i))*F(i,j)*sin(s(i)-s(j))*cos(p(j));
dsij=@(i,j) F(i,j)*(sin(p(i))*cos(s(i)-s(j))-cos(p(i))*sin(p(j)));

dp=sym(zeros(4,1));
ds=sym(zeros(4,1));

for iNode=1:4
    for jNode=1:4
        if iNode~=jNode
            dp(iNode)=dp(iNode)+dpij(iNode,jNode);
            ds(iNode)=ds(iNode)+dsij(iNode,jNode);
        end
    end
end
dz=[dp(1) ds(1) dp(2) ds(2) dp(3) ds(3) dp(4) ds(4)];
A1=sym(zeros(8));
A2=sym(zeros(8));
for iVar=1:8
    for jVar=1:8
        A1(iVar,jVar)=diff(dz(iVar),z(jVar));
        for iNode=1:4
            for jNode=1:4
                A2(iVar,jVar)=A2(iVar,jVar)+diff(dz(iVar),F(iNode,jNode))*diff(tij(iNode,jNode),z(jVar));
            end
        end
    end
end
for iNode=1:4
    for jNode=1:4
        A1=subs(A1,F(iNode,jNode),f);
    end
end
A1=simplify(A1/f);
%A2=simplify(A2);
save([mfilename '_data'])
%%
