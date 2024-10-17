function yMean=sphere_mean(y,varargin)
yMean=lie_mean(y, sphere_funs(), varargin{:});

% maxIt=100;
% 
% N=size(y,2);
% ymean=y(:,1);
% fPrev=Inf;
% for(it=1:maxIt)
%     hy=sphere_log(ymean,y);
%     ymean=sphere_exp(ymean,sum(hy,2)/N);
%     fCurrent=cost(ymean,y);
% %    fPrev-fCurrent
%     if((fPrev-fCurrent)<1e-12)
%         break;
%     end
%     fPrev=fCurrent;
% %    allf(it)=cost(ymean,y);
%     
% end
% %semilogy(allf)
% 
% function f=cost(ymean,y)
% f=sum(sphere_dist(ymean,y).^2);