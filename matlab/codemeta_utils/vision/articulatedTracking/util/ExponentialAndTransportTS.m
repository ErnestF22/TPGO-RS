function [p,u,v,w,U] = ExponentialAndTransportTS(p,u,v,w,U,L)

Nvec = size(U,2);

tol = 10^(-5);

if abs(norm(p,'fro')-1)>tol
    abs(norm(p)-1)
    error('Not a unit vector')
end

 if nargin<6
     L = 10;
 end
 stepsize = 1/L;

Z = U(1:3,:);
W = U(4:6,:);

  
 for i=1:L
     
     
     pprev = p;
     uprev = u;
     vprev = v;
     wprev = w;
     
     Zprev = Z;
     Wprev = W;
     
     
     p = ExpMapSphere(p,stepsize*v);
     u =  ParallelTranspot(pprev,p,uprev+stepsize*wprev);
     v =  ParallelTranspot(pprev,p,vprev-stepsize*CurvatureR(uprev,wprev,vprev));
     w = ParallelTranspot(pprev,p,wprev);
     
     
     for j=1:Nvec

         Z(:,j)=  ParallelTranspot(pprev,p, Zprev(:,j)-...
                                .5*stepsize*CurvatureR(uprev,Wprev(:,j),vprev)-...
                                .5*stepsize*CurvatureR(uprev,wprev,Zprev(:,j)));
                            
         W(:,j) = ParallelTranspot(pprev,p, Wprev(:,j)+ .5*stepsize*CurvatureR(vprev,Zprev(:,j),uprev));
     end
     
 end


U = [Z;W];
 

end


