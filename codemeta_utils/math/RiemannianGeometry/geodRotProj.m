%function [S,tfinal]=geodRotProj(R,R0,v,mode)
%Finds S, the projection of a rotation R along a geodesic passing through
%R0 and with tangent vector R0*v
function [S,tfinal]=geodRotProj(R,R0,v,mode)
v=v/norm(v);

if(exist('mode','var')==0)
    mode='closed';
end

switch(mode)
    case 'iter'
        gProj=Inf;
        t=0;
        while(abs(gProj)>1e-15)
            R1=R0*rot(t*v);
            g=logrot(R1'*R);
            gProj=v'*g;
            t=t+gProj;
        end

        tfinal=t;
        if(tfinal<0)
            tfinal=tfinal+2*pi;
        end
    case 'closed'
        tfinal=optimalt(R'*R0,v);
end
        
%S=R0*rot(tfinal*v);
S=R0*rot_exp(eye(3),tfinal*hat(v));

% T=linspace(0,2*pi,200);
% %T=linspace(5.3,5.7,1000);
% d=zeros(size(T));
% for(it=1:length(T))
%     t=T(it);
%     R1=rot(t*v);
%     d(it)=rot_dist(R0'*R,R1);
%     p(it)=logrot(R'*R1)'*v;
%     p2(it)=scaledTangentInnerProd(R'*R0,v,t);
% end
% 
% [mind,mindit]=min(d);
% tex=T(mindit)
% 
% plot(T,d,T,p,T,p2)
% grid on


function p=scaledTangentInnerProd(R,v,t)
vhat=hat(v);
Rd=R+sin(t)*R*vhat+(1-cos(t))*R*vhat^2;
RdRd=R-R'+sin(t)*(R*vhat-vhat'*R')+(1-cos(t))*(R*vhat^2-vhat'^2*R');
vhatRdRd=vhat'*(R-R'+R*vhat^2-vhat'^2*R')+sin(t)*vhat'*(R*vhat-vhat'*R')+cos(t)*(-vhat'*(R*vhat^2-vhat'^2*R'));
a=trace(vhat'*(R-R'+R*vhat^2-vhat'^2*R'));
b=trace(vhat'*(R*vhat-vhat'*R'));
c=trace(-vhat'*(R*vhat^2-vhat'^2*R'));
p=a+b*sin(t)+c*cos(t);
%p=(c^2+b^2)*cos(t)^2-2*a*c*cos(t)+(a^2-b^2);

function [topt,topt1,topt2]=optimalt(R,v)
vhat=hat(v);

a=trace(vhat'*(R-R'+R*vhat^2-vhat'^2*R'));
b=trace(vhat'*(R*vhat-vhat'*R'));
c=trace(-vhat'*(R*vhat^2-vhat'^2*R'));

cost1=(a*c-sqrt(a^2*c^2-(a^2-b^2)*(c^2+b^2)))/(c^2+b^2);
cost2=(a*c+sqrt(a^2*c^2-(a^2-b^2)*(c^2+b^2)))/(c^2+b^2);

if cost1>1 || cost1<-1 || cost2>1 || cost2<-1
    error(['Invalid values for the cosine: cost1 '...
        num2str(cost1) ' cost2 ' num2str(cost2)])
end

%We have 4 solutions.
%First find correct inverse of cosine by checking in the original equation
t1=acos(cost1);
t2=2*pi-t1;
if(abs(a+b*sin(t1)+c*cos(t1))<abs(a+b*sin(t2)+c*cos(t2)))
    topt1=t1;
else
    topt1=t2;
end

t3=acos(cost2);
t4=2*pi-t3;
if(abs(a+b*sin(t3)+c*cos(t3))<abs(a+b*sin(t4)+c*cos(t4)))
    topt2=t3;
else
    topt2=t4;
end

%Then check what is the closest and the furthest points
if(rot_dist(R',rot(topt1*v))<rot_dist(R',rot(topt2*v)))
    topt=topt1;
else
    topt=topt2;
end

