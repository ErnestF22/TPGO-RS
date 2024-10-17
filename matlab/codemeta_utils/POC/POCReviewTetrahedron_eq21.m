function POCReviewTetrahedron_eq21
%%
z=sym('z',[6 1]);
o=[0 z(4) z(5) z(2) z(3) 0;
    0 0 z(6) z(1) 0 z(3);
    0 0 0 0 z(1) z(2);
    0 0 0 0 z(6) z(5)
    0 0 0 0 0 z(4);
    0 0 0 0 0 0];
o=o+o.';
h=z*z.'+diag(diag(z*z.'))-2*eye(6)-o;
h=h-fliplr(diag(diag(fliplr(h))));
disp(h)

zeq=1/3*ones(6,1);
heq=subs(h,z,-1/3*ones(6,1));
disp(heq)
disp(heq*ones(6,1))

%%
A1=sym(zeros(6));
for iVar=1:6
    A1(:,iVar)=diff(h,z(iVar))*ones(6,1);
end
A1=subs(A1,z,zeq);
Lambda=ones(6)-fliplr(eye(6));
z0eq=-1/3;
A1b=(z0eq-1)/2*Lambda-(7*z0eq-1)*eye(6);
disp([A1-A1b])

A2b=(z0eq^2-z0eq)*Lambda+(z0eq^2+z0eq-2)*eye(6);
disp([heq-A2b])
