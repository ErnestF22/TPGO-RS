function [p,u,v,w] = ExpMapTS(p,u,v,w,L)

tol = 10^(-5);

if abs(norm(p,'fro')-1)>tol
    abs(norm(p)-1)
    error('Not a unit vector')
end

if nargin<5
L = 50;
end
 
stepsize = 1/L;

 for i=1:L
     
     pprev = p;
     uprev = u;
     vprev = v;
     wprev = w;
     
     p = ExpMapSphere(p,stepsize*v);
     u =  ParallelTranspot(pprev,p,uprev+stepsize*wprev);
     v =  ParallelTranspot(pprev,p,vprev-stepsize*CurvatureR(uprev,wprev,vprev));
     w = ParallelTranspot(pprev,p,wprev);
     

 end




end

