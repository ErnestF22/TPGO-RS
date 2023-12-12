function POCRotationFollowVector
epsilon=0.1;
itMax=100;
k=3;
I=eye(3);
e3=I(:,k);
er=sphere_randn();

H=@(R) diag([-1 -1 1])*householderRotation(R'*er,k);
dRVec=@(R) -logrot(H(R));

f=@(R) sphere_dist(er,R*e3);

R=rot_randn();
fAll=zeros(1,itMax);
for it=1:itMax
    R=R*rot(epsilon*dRVec(R));
    fAll(it)=f(R);
end
plot(fAll)
