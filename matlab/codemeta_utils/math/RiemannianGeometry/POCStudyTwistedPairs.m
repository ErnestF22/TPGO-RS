function POCStudyTwistedPairs
T=cnormalize(randn(3,1));
R=rot_randn();
E=hat(T)*R;

e3=[0;0;1];
e3hat=hat(e3);
Pz=eye(3)-e3*e3';
Rzpihalf=rot(pi/2*e3);
Rzpi=rot(pi*e3);

Z=e3hat';
W=Rzpihalf;

[Ua,Sa,Va]=svd(E);
Ta1=vee(Ua*Z*Ua');
Ra1=Ua*W*Va';
Ta2=vee(Ua*Z*Ua');
Ra2=Ua*W'*Va';
disp([E hat(Ta1)*Ra1 hat(Ta2)*Ra2])


% R1a1=Ua*Rzpihalf';
% R2a1=Va;
% disp([E R1a1*e3hat*R2a1'])
