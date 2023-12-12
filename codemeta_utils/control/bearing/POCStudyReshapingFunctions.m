function POCStudyReshapingFunctions
syms c
f=acos(c)^2/2;
display(f)
figure(1)
subplot(2,1,1)
ezplot(f,[-1 1])

df=diff(f);
display(df)
subplot(2,1,2)
ezplot(df,[-1 1])

dfn=-acos(c);
dfd=sqrt(1-c^2);
disp('diff(dfn)')
disp(diff(dfn))
disp('diff(dfd)')
disp(diff(dfd))
disp('simplify(diff(dfn)/diff(dfd))')
disp(simplify(diff(dfn)/diff(dfd)))
% figure(2)
% subplot(2,1,1)
% ezplot(diff(dfn),[-1 1])
% % ezplot(-acos(c),[-1 1])
% subplot(2,1,2)
% ezplot(diff(dfd),[-1 1])
% % ezplot(sqrt(1-c^2))

bs=simplify((f+(1-c)*df)/acos(c)*2);
bsb=acos(c)-2*sqrt(1-c)/sqrt(1+c);
display(bs)
figure(3)
check_fun(@(x) subs(bs,c,x), @(x) subs(bsb,c,x),linspace(-1,1))
subs(bsb,c,1)
disp('diff(bs)')
display(simplify(diff(bs)))
