function POCTrifocalRotationResiduals

switch 1
    case 1
        %derivatives
        n0=sphere_randn();
        n1=sphere_randn();
        n2=sphere_randn();

        [R1,~,~,~,v1Vec]=rot_randGeodFun([],'speed',1);
        [R2,~,~,~,v2Vec]=rot_randGeodFun([],'speed',rand);

        % e=@(t) n0'*hat(R1(t)*n1)*R2(t)*n2;
        % de=@(t) -n0'*hat(R1(t)*n1)*R2(t)*hat(n2)*v2Vec+n0'*hat(R2(t)*n2)*R1(t)*hat(n1)*v1Vec;

        v12Vec=[v1Vec;v2Vec];
        e=@(t) residuals(n0,R1(t),n1,R2(t),n2);
        ge=@(t) gradResiduals(n0,R1(t),n1,R2(t),n2);
        de=@(t) ge(t)'*v12Vec;
        Dge=@(t) DgradResiduals(n0,R1(t),n1,R2(t),n2);
        dge=@(t) Dge(t)*v12Vec;

        funs=consensus_rot3_almostGlobal_functions('type','l1l2');
        
        fe=@(t) funs.f(e(t));
        gfe=@(t) funs.df(e(t))*ge(t);
        dfe=@(t) gfe(t)'*v12Vec;
        Dgfe=@(t) funs.ddf(e(t))*ge(t)*ge(t)'+funs.df(e(t))*Dge(t);
        dgfe=@(t) Dgfe(t)*v12Vec;
        
        
        %check_der(e,de)
        %check_der(ge,dge)
        %check_der(fe,dfe)
        check_der(gfe,dgfe)
    case 2
        %summation of squares
        n0a=sphere_randn();
        n1a=sphere_randn();
        n0b=sphere_randn();
        n1b=sphere_randn();
        R=rot_randn();
        
        e1=(n0a'*R*n1a)^2;
        e2=n0a'*R*n1a*n1a'*R'*n0a;

        Aa=buildQuadratic(n0a,n1a);
        e3=R(:)'*Aa*R(:);
        
        disp([e1 e2 e3])

        e1=e1+(n0b'*R*n1b)^2;
        e2=e2+n0b'*R*n1b*n1b'*R'*n0b;

        Ab=buildQuadratic(n0b,n1b);
        e3=e3+R(:)'*Ab*R(:);

        disp([e1 e2 e3])

        Aa2=buildQuadratic2(n0a,n1a*n1a',n0a);
        Ab2=buildQuadratic2(n0b,n1b*n1b',n0b);
        
        disp([Aa-Aa2])
        disp([Ab-Ab2])
        
end

function A=buildQuadratic(n0,n1)
A=zeros(9);
for idx1=1:9
    E1=zeros(3);
    E1(idx1)=1;
    for idx2=1:9
        E2=zeros(3);
        E2(idx2)=1;
        A(idx1,idx2)=n0'*E1*n1*n1'*E2'*n0;
    end
end

function A=buildQuadratic2(a,B,c)
A=zeros(9);
for idx1=1:9
    E1=zeros(3);
    E1(idx1)=1;
    for idx2=1:9
        E2=zeros(3);
        E2(idx2)=1;
        A(idx1,idx2)=a'*E1*B*E2'*c;
    end
end

function e=residuals(n0,R1,n1,R2,n2)
e=n0'*hat(R1*n1)*R2*n2;

function grade=gradResiduals(n0,R1,n1,R2,n2)
grade=[hat(n1)*R1'*hat(R2*n2)*n0; -hat(n2)*R2'*hat(R1*n1)*n0];

function Dgrade=DgradResiduals(n0,R1,n1,R2,n2)
Dgrade=[
     hat(n1)*hat(R1'*hat(R2*n2)*n0)  hat(n1)*R1'*hat(n0)*R2*hat(n2);
    -hat(n2)*R2'*hat(n0)*R1*hat(n1) -hat(n2)*hat(R2'*hat(R1*n1)*n0)];
% -hat(n2)*R2'*hat(R1*n1)*n0];
