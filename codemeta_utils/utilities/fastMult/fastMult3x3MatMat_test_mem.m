function fastMult3x3MatMat_test_mem
expNum=4;
N=1e5;
A=randn(9,N);
B=randn(9,N);

switch expNum
    case 1
        C=fastMult3x3MatMat(A,B);
        disp(mean(C(:)))
        C=fastMult3x3MatMat(A,B,C);
        disp(mean(C(:)))
    case 2
        C=zeros(9,N);
        C=fastMult3x3MatMat(A,B,C);
        disp(mean(C(:)))
        C=fastMult3x3MatMat(A,B,C);
        disp(mean(C(:)))
    case 3
        C=randn(9,N);
        D=randn(9,N);
        [C,D]=f1(A,B,C,D);
        disp(mean(C(:)))
        [C,D]=f1(A,B,C,D);
        disp(mean(C(:)))
    case 4
        C=randn(9,N);
        D=randn(9,N);
        D=f2(A,B,C,D);
        disp(mean(D(:)))
        D=f2(A,B,C,D);
        disp(mean(D(:)))
        D=f2(A,B,C,D);
        disp(mean(D(:)))
end

function [C,D]=f1(A,B,C,D)
D=fastMult3x3MatMat(C,B,D);
C=fastMult3x3MatMat(D,A,C);

function D=f2(A,B,C,D)
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
