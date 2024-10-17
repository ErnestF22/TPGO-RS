function bidiagonalize3x3_test
A=randn(3);
[D,U,V]=bidiagonalize3x3(A);
disp([D U*(A*V)])
