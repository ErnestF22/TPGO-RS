function POCSolutionBoundsSecondOrder
testNumber=2;
switch testNumber
    case 1
        kv=1+rand;
        kx=1+rand;
        %a=0.1;
        a=2*kx*kv/(kv^2+4*kx);
    case 2
        syms kv kx a
end
   
A=[-kv -kx; 1 0];
P=[1 a; a kx];
f=@(t,z) A*z;

switch testNumber
    case 1
        z0=[0;1];
        [t,z]=ode45(f,[0,20],z0);

        %plot(t,x)
        z=z';
    case 2
        syms x dx
        z=[dx;x];
end

%E=0.5*sum(z.*(P*z));
Q=-(P*A+A.'*P.')/2;
disp('Q=')
disp(Q)

dE=sum(z.*(P*f(0,z)));
dEb=-sum(z.*(Q*z));

switch testNumber
    case 1
        plot(t,dE,'b',t,dEb,'rx')
        title('Check that derivative computed from Q is correct')
    case 2
        disp('Check that derivative computed from Q is correct')
        disp(simplify(dE-dEb))
end

Qb=[kv-a a*kv/2; a*kv/2 kx*a];

disp('Q-Qb')
disp(Q-Qb)

switch testNumber
    case 1
        disp('eig(Q)')
        disp(eig(Q))
    case 2
        rhCondition=simplify(4*det(Q)/a)>0;
        disp('Rough-Hurwitz conditions')
        disp('kx,kv,a > 0')
        disp(rhCondition)
end

Qneg=[kv-a -a*kv/2; -a*kv/2 kx*a];
disp('Compare evals of Q and Qneg')
disp(eig(Q)-eig(Qneg))

switch testNumber
    case 1
        l=@(a) eig([kv-a a*kv/2; a*kv/2 kx*a]);
        plotfun(l,linspace(0,4*kx*kv/(kv^2+4*kx)))
        aMax=(4*kv*kx + kv^2*(kv^2 + 4*kx)^(1/2) + 4*kv*kx^2 + kv^3*kx + kv^3 - kv^2*kx*(kv^2 + 4*kx)^(1/2))/(kv^4 + kv^2*kx^2 + 6*kv^2*kx + kv^2 + 4*kx^3 + 8*kx^2 + 4*kx);
        ax=axis();
        hold on
        plot([aMax aMax],[ax(3) ax(4)],'k')
        hold off
    case 2
        l=eig(Q);
        solve(simplify(diff(l(1),a)==0),a)
        
end

