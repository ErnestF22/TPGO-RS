function NotesPlots

plotnum=1;

t=linspace(1e-12,pi,100);
t2=linspace(1e-8,pi/2,100);

switch plotnum
    case 1

        f=@(t) t/sin(t)*(cos(t)+1);
        fe=@(t) 2;
        fh=@(t) t/sinh(t)*(cosh(t)+1);
        allf=@(t) [f(t);fe(t);fh(t)];

        figure(1)
        plotfun(allf,t)
        legend('\kappa=1','\kappa=0','\kappa=-1','Location','EastOutside')
        xlabel('l')
        ylabel('l(C_\kappa(l)+1)S^-^1_\kappa(l)','Interpreter','tex')
        %axis equal
        
        savefigure('../vision/cameranetwork/TAC11-manifolds/images/mumax-wide','epsc',[500 120])
        %savefigure('../vision/cameranetwork/cdc10-manifolds/mumax','epsc',[320 240])

    case 2
        figure(2)
        fn=@(t,k) t/(sinh(sqrt(abs(k))*t)/sqrt(abs(k)));
        fp=@(t,k) t/(sin(sqrt(k)*t)/sqrt(k));
        allf=@(t) [fp(t,2) fp(t,1) fp(t,0.5) fn(t,0) fn(t,-0.5) fn(t,-1) fn(t,-2)];
        plotfun(allf,t2);
        legend('k=2','k=1','k=0.5','k=0','k=-0.5','k=-1','k=-2')

        figure(3)
        fn=@(t,k) t*cos(sqrt(abs(k))*t)/(sinh(sqrt(abs(k))*t)/sqrt(abs(k)));
        fp=@(t,k) t*cos(sqrt(abs(k))*t)/(sin(sqrt(k)*t)/sqrt(k));
        allf=@(t) [fp(t,2) fp(t,1) fp(t,0.5) fn(t,0) fn(t,-0.5) fn(t,-1) fn(t,-2)];
        plotfun(allf,t2);
        legend('k=2','k=1','k=0.5','k=0','k=-0.5','k=-1','k=-2')
        
    case 3
        delta=0;
        Delta=2;
        dmax=pi/2;

        sd=@(t) sdelta(delta,t);
        sD=@(t) sdelta(Delta,t);
        cd=@(t) cdelta(delta,t);
        
        mumax=@(t) t*cd(t)/sd(t)+t/sD(t);
        plotfun(mumax,t2)
        
        disp(mumax(dmax));
        
        
end

function v=sdelta(delta,t)
if delta>0
    v=sin(sqrt(delta)*t)/sqrt(delta);
elseif delta==0
    v=t;
else
    v=sinh(sqrt(-delta)*t)/sqrt(-delta);
end

function v=cdelta(delta,t)
if delta>0
    v=cos(sqrt(delta)*t);
elseif delta==0
    v=1;
else
    v=cosh(sqrt(-delta)*t);
end
    
