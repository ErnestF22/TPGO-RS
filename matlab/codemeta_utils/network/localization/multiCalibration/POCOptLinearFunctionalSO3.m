function POCOptLinearFunctionalSO3
resetRands()
A=randn(20,9);
RTruth=rot_randn();
%b=A*RTruth(:);
b=zeros(size(A,1),1);

E=eye(3);
EHat=zeros(9,3);
for k=1:3
    EHat(:,k)=reshape(hat(E(:,k)),9,1);
end

switch 1
    case 1
        [R,~,~,~,vVec]=rot_randGeodFun(RTruth);

        RVec=@(R) R(:);
        f=@(R) A*R(:)-b;
        c=@(R) 0.5*f(R)'*f(R);
        df=@(R,v) A*kron(eye(3),R)*EHat*v;
        %dc=@(R,v) f(R)'*df(R,v);
        g=@(R) EHat'*kron(eye(3),R')*A'*f(R);
        dc=@(R,v) g(R)'*v;
        %H=@(R,v) EHat'*kron(eye(3),hat(v)'*R')*A'*f(R)+EHat'*kron(eye(3),R')*A'*A*kron(eye(3),R)*EHat*v;
        H=@(R,v) hessian(A,b,EHat,R)*v;
        dRVec=@(R,v) kron(eye(3),R)*EHat*v;
        ddc=@(R,v) H(R,v)'*v;

        %check_der(@(t) RVec(R(t)), @(t) dRVec(R(t),vVec))
        %check_der(@(t) f(R(t)), @(t) df(R(t),vVec))
        %check_der(@(t) c(R(t)), @(t) dc(R(t),vVec))
        check_der(@(t) dc(R(t),vVec), @(t) ddc(R(t),vVec))
    case 2
        load rot_minimize_badData
        %R=eye(3);
        R=R0;
        f=@(R) A*R(:)-b;
        c=@(R) 0.5*f(R)'*f(R);
        NIt=100;
        allc=zeros(NIt,1);
        allc(1)=c(R);
        for it=2:NIt
            v=-(A*kron(eye(3),R)*EHat)\f(R);
            R=R*rot(0.5*v);
            allc(it)=c(R);
            disp(norm(v))
            if norm(v)<1e-14
                break;
            end
        end
        semilogy(allc)
end

function H=hessian(A,b,EHat,R)
n=size(R,1);
d=size(EHat,2);
f=A*R(:)-b;
H=zeros(d,d);
for k=1:d
    H(:,k)=EHat'*kron(eye(3),reshape(EHat(:,k),n,n)'*R')*A'*f;
end
H=H+EHat'*kron(eye(3),R')*A'*A*kron(eye(3),R)*EHat;
