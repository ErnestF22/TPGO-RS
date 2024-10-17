%function D=fastMult3x3MatMatMat(A,B,C,D)
%Similar to fastMult3x3MatMat, but computes D=A*B*C.
%
%It makes the variable for intermediate computations persistent, so that
%memory is allocated only at the first call
function D=fastMult3x3MatMatMat(A,B,C,D)
persistent E

E=[...
    A(1:9:end).*B(1:9:end)+A(4:9:end).*B(2:9:end)+A(7:9:end).*B(3:9:end);
    A(2:9:end).*B(1:9:end)+A(5:9:end).*B(2:9:end)+A(8:9:end).*B(3:9:end);
    A(3:9:end).*B(1:9:end)+A(6:9:end).*B(2:9:end)+A(9:9:end).*B(3:9:end);
    A(1:9:end).*B(4:9:end)+A(4:9:end).*B(5:9:end)+A(7:9:end).*B(6:9:end);
    A(2:9:end).*B(4:9:end)+A(5:9:end).*B(5:9:end)+A(8:9:end).*B(6:9:end);
    A(3:9:end).*B(4:9:end)+A(6:9:end).*B(5:9:end)+A(9:9:end).*B(6:9:end);
    A(1:9:end).*B(7:9:end)+A(4:9:end).*B(8:9:end)+A(7:9:end).*B(9:9:end);
    A(2:9:end).*B(7:9:end)+A(5:9:end).*B(8:9:end)+A(8:9:end).*B(9:9:end);
    A(3:9:end).*B(7:9:end)+A(6:9:end).*B(8:9:end)+A(9:9:end).*B(9:9:end);
   ];
D=[...
    E(1:9:end).*C(1:9:end)+E(4:9:end).*C(2:9:end)+E(7:9:end).*C(3:9:end);
    E(2:9:end).*C(1:9:end)+E(5:9:end).*C(2:9:end)+E(8:9:end).*C(3:9:end);
    E(3:9:end).*C(1:9:end)+E(6:9:end).*C(2:9:end)+E(9:9:end).*C(3:9:end);
    E(1:9:end).*C(4:9:end)+E(4:9:end).*C(5:9:end)+E(7:9:end).*C(6:9:end);
    E(2:9:end).*C(4:9:end)+E(5:9:end).*C(5:9:end)+E(8:9:end).*C(6:9:end);
    E(3:9:end).*C(4:9:end)+E(6:9:end).*C(5:9:end)+E(9:9:end).*C(6:9:end);
    E(1:9:end).*C(7:9:end)+E(4:9:end).*C(8:9:end)+E(7:9:end).*C(9:9:end);
    E(2:9:end).*C(7:9:end)+E(5:9:end).*C(8:9:end)+E(8:9:end).*C(9:9:end);
    E(3:9:end).*C(7:9:end)+E(6:9:end).*C(8:9:end)+E(9:9:end).*C(9:9:end);
   ];
