%function D=fastMult3x3MatMatMat(A,B,C,D)
%Similar to fastMult3x3MatMat, but computes 
%   D(:,:,j)=A(:,:,eA(j))'*B(:,:,j)*C(:,:,eC(j)).

%It makes the variable for intermediate computations persistent, so that
%(ideally) memory is allocated only at the first call. But it not always
%work
function D=fastMult3x3MatIdxTranspMatMatIdx(A,eA,B,C,eC,D)
persistent E

eA=(eA-1)*9;
eC=(eC-1)*9;

E=[...
    A(eA+1).*B(1:9:end)+A(eA+2).*B(2:9:end)+A(eA+3).*B(3:9:end);
    A(eA+4).*B(1:9:end)+A(eA+5).*B(2:9:end)+A(eA+6).*B(3:9:end);
    A(eA+7).*B(1:9:end)+A(eA+8).*B(2:9:end)+A(eA+9).*B(3:9:end);
    A(eA+1).*B(4:9:end)+A(eA+2).*B(5:9:end)+A(eA+3).*B(6:9:end);
    A(eA+4).*B(4:9:end)+A(eA+5).*B(5:9:end)+A(eA+6).*B(6:9:end);
    A(eA+7).*B(4:9:end)+A(eA+8).*B(5:9:end)+A(eA+9).*B(6:9:end);
    A(eA+1).*B(7:9:end)+A(eA+2).*B(8:9:end)+A(eA+3).*B(9:9:end);
    A(eA+4).*B(7:9:end)+A(eA+5).*B(8:9:end)+A(eA+6).*B(9:9:end);
    A(eA+7).*B(7:9:end)+A(eA+8).*B(8:9:end)+A(eA+9).*B(9:9:end);
   ];
D=[...
    E(1:9:end).*C(eC+1)+E(4:9:end).*C(eC+2)+E(7:9:end).*C(eC+3);
    E(2:9:end).*C(eC+1)+E(5:9:end).*C(eC+2)+E(8:9:end).*C(eC+3);
    E(3:9:end).*C(eC+1)+E(6:9:end).*C(eC+2)+E(9:9:end).*C(eC+3);
    E(1:9:end).*C(eC+4)+E(4:9:end).*C(eC+5)+E(7:9:end).*C(eC+6);
    E(2:9:end).*C(eC+4)+E(5:9:end).*C(eC+5)+E(8:9:end).*C(eC+6);
    E(3:9:end).*C(eC+4)+E(6:9:end).*C(eC+5)+E(9:9:end).*C(eC+6);
    E(1:9:end).*C(eC+7)+E(4:9:end).*C(eC+8)+E(7:9:end).*C(eC+9);
    E(2:9:end).*C(eC+7)+E(5:9:end).*C(eC+8)+E(8:9:end).*C(eC+9);
    E(3:9:end).*C(eC+7)+E(6:9:end).*C(eC+8)+E(9:9:end).*C(eC+9);
   ];
