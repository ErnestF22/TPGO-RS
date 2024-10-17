function essential_dist_test
switch 3
    case 1
        %test zero distance when only vertical motion
        plotfuntrials(@f1,1000)
    case 2
        %test commutativity
        plotfuntrials(@f2,1000)
    case 3
        %distance minimizing?
        %resetRands()
        NPoints=101;
        tMax=pi;

        t=linspace(0,tMax,NPoints)';
        d=zeros(NPoints,1);
        tMin=zeros(NPoints,1);
        tMinExp=zeros(NPoints,1);
        
        Q0=essential_randn();
        v0=essential_randTangentNormVector(Q0);
        
        for it=1:NPoints
            Q=essential_exp(Q0,t(it)*v0);
            [Q,tMinExp(it)]=essential_randomVerticalMotion(Q);
            [d(it),tMin(it)]=essential_dist(Q0,Q);
        end
        subplot(2,1,1)
        plot(t,d,t,t,'x')
        title('Distance vs. traveled distance')
        subplot(2,1,2)
        plot(t,tMin,t,-tMinExp,'x')
        title('Recovered rotation in equivalence rotation vs. negative of applied one')
        
end

function c=f1()
Q0=essential_randn();
Q=essential_randomVerticalMotion(Q0);
c=essential_dist(Q,Q0);

function c=f2()
Q1=essential_randn();
Q2=essential_randn();

c=essential_dist(Q1,Q2)-essential_dist(Q2,Q1);
