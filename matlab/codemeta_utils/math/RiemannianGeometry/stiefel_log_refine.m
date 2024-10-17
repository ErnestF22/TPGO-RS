function [h12,err]=stiefel_log_refine(y1,y2)
[n,p]=size(y1);

%make Y1 and Y2 square matrices in SO(n)
YY1=orthCompleteBasis(y1);
YY2=orthCompleteBasis(y2);
%move YY2 inside the equivalence class in order to minimize
%distance to YY1
[YY2]=stiefel_log_findOrthGrad(YY1,YY2,p);
%compute the log map
H=YY1*rot_log([],YY1'*YY2);   %projection of H on the vertical space should be zero here
H=projTangent(YY1,H);

for(it=1:30)
    g=YY2-rot_exp(YY1,H);
    geodCost=@(t) cost(YY1,YY2,projTangent(YY1,H+t*g));
    [t,costval]=fminbnd(geodCost,-0.1,0.1);

    H=projTangent(YY1,H+0.01*t*g);
    
    d(it)=cost(YY1,YY2,H);
    ng(it)=norm(projTangent(YY1,g));
    allt(it)=t;
end

h12=H(:,1:p);

% subplot(3,1,1)
% plot(d)
% title('||y2-stiefelexp(y1,h12)||')
% subplot(3,1,2)
% plot(ng)
% title('||proj(...)||')
% subplot(3,1,3)
% plot(allt)
% title(t)
err.d=d;
err.ng=ng;
err.allt=allt;

function c=cost(YY1,YY2,H)
[n,p]=size(H);
err=YY2-rot_exp(YY1,H);
c=sum(sum(err(:,1:p).^2));

function H=projTangent(YY1,H)
[n,p]=size(H);
H=rot_tangentProj(YY1,H);
E=[zeros(n-p,p) eye(n-p)]';
H=H-YY1*E*E'*YY1'*H*E*E';
