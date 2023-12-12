function y1=sphere_exp(y,h)
%flag to control tricks to improve numerical precision
flagNormalize=true;         %re-normalize result
flagTangentProject=true;    %re-project tangent vectors to the tangent space
flagRotMult=false;

if flagTangentProject
    h=sphere_tangentProj(y,h);
end

flagPermute=false;
if(size(h,2)==1)
    h=permute(h,[1 3 2]);
    flagPermute=true;
end

[D,N]=size(h);

y1=zeros(size(y));
for iN=1:N
    [hnorm,theta]=qr(h(:,iN),0);
    if flagRotMult
        R=[cos(theta) -sin(theta); sin(theta) cos(theta)];
        y1(:,iN)=([y(:,iN) hnorm]*R+[y(:,iN) hnorm]/R')*[0.5;0];
    else
        y1(:,iN)=y(:,iN)*cos(theta)+hnorm*sin(theta);
    end
end
% if(nargout>1)
%     h1=-y*sin(theta)*theta+h*cos(theta);
% end

if(flagPermute)
    y1=permute(y1,[1 3 2]);
%     if(nargout>1)
%         h1=permute(h1,[1 3 2]);
%     end
end

if flagNormalize
    y1=cnormalize(y1);
end
