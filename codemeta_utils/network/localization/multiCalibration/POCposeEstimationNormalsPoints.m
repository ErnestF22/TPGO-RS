function POCposeEstimationNormalsPoints
load poseEstimationFromNormalsPoints_testData
flagSubsample=true;
flagReorder=false;

if flagReorder
    %reorder samples
    n=[n(:,1:2:end) n(:,2:2:end)];
    x=[x(:,1:2:end) x(:,2:2:end)];
    d=[d(1:2:end) d(2:2:end)];
end

if flagSubsample
    NPoints=9;
    n=n(:,1:NPoints);
    x=x(:,1:NPoints);
    d=d(1:NPoints);
end

M=n';
[UM,SM,VM]=svd(M,'econ');
NPoints=size(M,1);
nx=zeros(9,NPoints);
nx2=zeros(6,NPoints);
for inx=1:NPoints
    nx(:,inx)=reshape(n(:,inx)*x(:,inx)',9,1);
    nx2(:,inx)=reshape(n(:,inx)*x(2:3,inx)',6,1);
end
%b=@(R) (n.*(R*x))'*ones(3,1)-d';
%disp(n'*T+b(R))
%TfromR=@(R) -M'*M\(M'*b(R));
%disp(n'*TfromR(R)+b(R))
%disp(b(R)-U*U'*b(R))
%b=@(R) (n.*(R*x))'*ones(3,1);
p=d'-UM*(UM'*d');
%b=@(R) nx'*R(:);
%disp(b(R)-U*U'*b(R)-p)
%A=nx'-U*(U'*nx');
%disp(A*R(:)-p)
A2=nx2'-UM*(UM'*nx2');
disp(A2*R(4:9)'-p)
RPartEst=reshape(A2\p,3,2);
[UR,SR,VR]=svd([zeros(3,1) RPartEst]);
REst=UR*diag([1 1 det(UR*VR')])*VR';
disp('[R REst R-REst]')
disp([R REst R-REst])
%TEst=-(M'*M)\(M'*(nx2'*REst(4:9)'-d'));
TEst=-VM*diag(diag(SM).^-1)*UM'*(nx2'*REst(4:9)'-d');
disp([T TEst])