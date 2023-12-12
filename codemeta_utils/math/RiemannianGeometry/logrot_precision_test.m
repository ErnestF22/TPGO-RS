function logrot_precision_test

randn('seed',0)

Ntrials=1000;
t=logspace(-6,0,100);

totError=zeros(size(t));
for(itrial=1:Ntrials)
    w=randn(3,1);
    w=w/norm(w);
    for(itheta=1:length(t))
        wt=t(itheta)*w;
        wtRec=logrot(rot(wt));
        
        wt=wt/norm(wt);
        wtRec=wtRec/norm(wtRec);
        
        totError(itheta)=totError(itheta)+abs(1-wt'*wtRec);
    end
    if(mod(itrial,100)==0)
        disp(itrial)
    end
end
totError=totError/Ntrials;

semilogy(t,totError);