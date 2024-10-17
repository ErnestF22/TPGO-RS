function detLeibniz_test
dMax=8;
for d=1:dMax
    A=randn(d);
    disp([det(A) detLeibniz(A)])
end
