function POCSecondOrderPositionControlLyapunovDerivative
A1=@(bx) [0 bx;bx 0];
A2=@(a1,a2) [2*a2 a1; a1 0];
A3=@(a1,a2,mx) [-2*mx a2; a2 2*a1];
B=@(bx,kx,cx,a1,a2,mx) A1(bx)+kx*A2(a1,a2)+cx*A3(a1,a2,mx);

syms bx kx cx a1 a2 mx
Bsym=B(bx,kx,cx,a1,a2,mx);
disp(Bsym)
%disp(eig(Bsym))




% figure(1)
% funPlot(@(t) eig(B(3.2,-1.25+t,-0.4,1,1)),linspace(-1,1))
% grid on
% 
% eig(B(1,3.2,-1.25,-0.4,1,1))

mx=5;
bx=1;
cx=1;
a1=1;
kx=1;
a2=2*cx*mx/kx;

while bx-kx*a1-cx*a1>0
    kx=2*kx;
    cx=2*cx;
end

disp(B(bx,kx,cx,a1,a2,mx))
