function plotController(yt,y,k,style)
Xx=min(y(1,:)):1:max(y(1,:));
Yy=min(y(2,:)):1:max(y(2,:));
minn=100;
for i=1:size(Xx,2)
    for j=1:size(Yy,2)
        p=[Xx(i);Yy(j)];
        [flagPoints]=polygon_isCollision(y,p);
        if flagPoints==0
            Y=(yt-p*ones(1,size(yt,2)));
            for s=1:size(Y,2)
                M(2*s-1,1)=Y(1,s);
                M(2*s,1)=Y(2,s);
            end
            Y=M;
            u=k*Y;
            % for NotVectorized function comment line 11 to 16 and 
            % uncomment line 19
%             u=Y*k;%for notvectorized version
            if norm(u)<minn
                minn=norm(u);
            end
            quiver(p(1),p(2),u(1),u(2),'Color', [0.75 0.75 0.75])
            hold on
        end
    end
end