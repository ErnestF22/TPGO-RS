function POCPointMassRouthHurwitz
kv=sym('kv');
kx=sym('kx');
cx=sym('cx');
m=sym('m');
mux=sym('mux');

Wx1=[kx*cx/m kv*cx/(2*m);kv*cx/(2*m) kv-mux*cx];
disp('Wx1')
evaluateRh(Wx1)
Wx2=[kx*mux cx; cx m];
disp('Wx2')
evaluateRh(Wx2)

    function evaluateRh(Wx1)
    disp('det()')
    disp(det(Wx1))
    disp('Condition on cx, rh(cx)<0')
    rh=collect(-4*m^2*det(Wx1)/cx,cx);
    disp(rh<0)
    cxMax= 4*kv*kx*m/(kv^2 + 4*kx*m*mux);
    disp('Value for rh(cxMax)')
    disp(subs(rh,cx,cxMax))
    end
end