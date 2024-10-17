function POCHomFlow4ptBilinearRefine
resetRands()
[X,G,NVec,idxX]=homFlowDatasetStructure(1);
NX=size(X,2);
[x,dx,v,w]=homFlowDatasetFlow(X,G);

NVecinG=rigidTransformG(G,NVec,'references','planes','wc');
nVecinG=planeNVecToNScaled(NVecinG);
n=nVecinG;
nv=n*v';

preAv=[];
preAn=[];
A=[];
b=[];
for iX=1:NX
    xi=x(1,iX);
    yi=x(2,iX);
    %    1  2  3  4  5  6     7     8
%     Ai=[ 1 xi yi  0  0  0 xi*yi  xi^2;
%          0  0  0  1 xi yi  yi^2 xi*yi;
%          ]*VNToAlphaMat;
    A2i=[eye(2) -x(:,iX)]';
    A1i=[x(:,iX)' 1];
    Ai=kron(A2i',A1i);
    bi=dx(:,iX)-[-w(2)+yi*w(3)+xi*yi*w(1)-xi^2*w(2);
        w(1)-xi*w(3)+yi^2*w(1)-xi*yi*w(2)];
    
    %disp(A1i*nv*A2i-Mi)
    Ain=A2i'*v*A1i;
    Aiv=A2i'*(A1i*n);
    %disp(Ain*n-Mi)
    %disp(Aiv*v-Mi)
    
    preAn=[preAn;kron(A1i',A2i')];
    preAv=[preAv;kron(vec(A2i'),A1i)];
    A=[A;Ai];
    b=[b;bi];
end

AnVec=reshape(preAn*v,2,3,[]);
An=reshape(permute(AnVec,[1 3 2]),[],3);
AvVec=reshape(preAv*n,2,3,[]);
Av=reshape(permute(AvVec,[1 3 2]),[],3);
disp(max(abs(An*n-b)))
disp(max(abs(Av*v-b)))
disp(max(abs(A*vec(n*v')-b)))

dx=dx+0.05*randn(size(dx));
nTruth=n;
vTruth=v;
[n,v]=homFlowParametersEstimate4pt(x,dx,w);

E=energy(n,v,A,b);
for it=1:20
    AnVec=reshape(preAn*v,2,3,[]);
    An=reshape(permute(AnVec,[1 3 2]),[],3);
    n=An\b;
    E=[E energy(n,v,A,b)];
    AvVec=reshape(preAv*n,2,3,[]);
    Av=reshape(permute(AvVec,[1 3 2]),[],3);
    v=Av\b;
    E=[E energy(n,v,A,b)];
end

normv=norm(v);
v=v/normv;
n=n*normv;
semilogy(E)
disp([nTruth n; vTruth v])

function E=energy(n,v,A,M)
E=norm(A*vec(n*v')-M)^2;