% Calculate gradient of CLF
clear all;

% Define CLF
syms x1 x2 x3 x1s x2s x3s real
Zc = [x1;x2];
Zs = [x1s;x2s];
V_phi = ((Zc-Zs)/norm(Zc-Zs)-[cos(x3-x3s);sin(x3-x3s)])'*((Zc-Zs)/norm(Zc-Zs)-[cos(x3-x3s);sin(x3-x3s)]);
V_p = (Zc-Zs)'*(Zc-Zs);
CLF = V_phi+V_p;
D_CLF_x1 = diff(CLF,x1);
D_CLF_x2 = diff(CLF,x2);
D_CLF_x3 = diff(CLF,x3);