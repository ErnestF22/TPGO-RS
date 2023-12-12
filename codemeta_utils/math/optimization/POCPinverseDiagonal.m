function POCPinverseDiagonal
d=5;
u=cnormalize(ones(d,1));

H=orthComplementProjector(u);
D=diag(1:d);
DInv=inv(D);

X=DInv*H*DInv;

disp('[H-D*X*D]')
disp([H-D*X*D])

disp('[pinv(D*X*D) DInv*pinv(X)*DInv]')
disp([pinv(D*X*D) DInv*pinv(X)*DInv])

disp('difference E')
E=pinv(D*X*D)-DInv*pinv(X)*DInv;
disp(E)


PH=orthComplementProjector(null(H));
NH=null(H);
Id=eye(d);
In=eye(size(NH,2));

disp('PH*E*PH')
disp(PH*E*PH);




