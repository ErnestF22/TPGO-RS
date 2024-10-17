function POCProcrustesContinuity
x0=[-1 0 1;
     0 1 0];
v=[  0  0 0;
     1 -1 1];
 
x=@(t) x0+t*v;
 
% f=@(t) [norm(center(x(t))-x0,'fro').^2; norm(-center(x(t))-x0,'fro').^2];
% 
% figure(1)
% funPlot(f,linspace(0,4))

f=@(t) procrustes(x0',x(t)','Scaling',false,'Reflection',false);
figure(2)
ezplot(f,[1.5 2.5])

end

function x=center(x)
x=x-mean(x,2)*ones(1,size(x,2));
end

function c=procrustesCost(x0,x)
c=procrustes(x0',x');
end
