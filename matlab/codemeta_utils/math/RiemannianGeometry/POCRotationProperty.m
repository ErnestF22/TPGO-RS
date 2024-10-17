function POCRotationProperty
%R=rot_randn();
%R=sym('R',[3 3]);
%A=[0 -1 0; 1 0 0; 0 0 1];
%A=[0 1 0; 1 0 0; 0 0 1];

%t1=rand*pi;
%t2=rand*pi;
%t3=rand*pi;
syms t1 t2 t3
c1=cos(t1);
c2=cos(t2);
c3=cos(t3);
s1=sin(t1);
s2=sin(t2);
s3=sin(t3);
R=[ c2      -c3*s2          s2*s3;
    c1*s2   c1*c2*c3-s1*s3  -c3*s1-c1*c2*s3;
    s1*s2   c1*s3+c2*c3*s1  c1*c3-c2*s1*s3];

%disp(R.'*R)

ABig=zeros(9,9);
ABig(1,1)=1;
ABig(5,1)=1;
ABig(1,5)=1;
ABig(2,4)=-1;
ABig(4,2)=-1;
ABig(5,5)=1;
ABig(2,2)=1;
ABig(4,4)=1;
ABig(9,9)=-1;
c1=R(1,1)+R(2,2);
c2=R(1,2)-R(2,1);
m=R(1,1)^2+R(2,2)^2+2*R(1,1)*R(2,2)+R(1,2)^2+R(2,1)^2-2*R(1,2)*R(2,1);
c3=R(3,3);

p=R(1,1)^2+R(2,2)^2+2*R(1,1)*R(2,2)+R(1,2)^2+R(2,1)^2-2*R(1,2)*R(2,1)-R(3,3)^2-2*R(3,3);
%p2=trace(R.'*A*R);
p2=expand(trace(R(:).'*ABig*R(:)));
disp([p;p2;p-p2])
disp(ABig)
%disp(R*R.')
simplify(p)
