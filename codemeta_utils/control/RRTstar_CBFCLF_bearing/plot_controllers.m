function plot_controllers(tree,l1,l2,l3,l4)
X=0:10:100;
Y=0:10:100;
S = ["r*","b*","g*","m*","y*","c*","g*","r*"];
% y=[20 45;40 45;25 5;2 15]';
c=0;
for n=1:size(tree,2)
    if ~isempty(tree(n).parent) 
        y = tree(n).L;
        if y == 1
            y = l1;
        elseif y==2
            y = l2;
        elseif y==3
            y = l3;
        else
            y = l4;
        end
        k = tree(n).K;
        A = tree(n).convex.A;
        b = tree(n).convex.b;
        for i=1:size(X,2)
            for j=1:size(Y,2)
                x = [X(i);Y(j)];
                if all(A*x<=b)
                    v=(y-x*ones(1,size(y,2)));
                    Vy=(y);
                    Vx = x*ones(1,size(y,2));
                    M = zeros(2*size(y,2),1);
                    Ny=M;
                    Nx=M;
                    for s=1:size(v,2)
                        M(2*s-1,1)=v(1,s);
                        M(2*s,1)=v(2,s);
                        Ny(2*s-1,1)=Vy(1,s);
                        Ny(2*s,1)=Vy(2,s);
                        Nx(2*s-1,1)=Vx(1,s);
                        Nx(2*s,1)=Vx(2,s);
                    end
                    v=M;
                    u=k*v;
                    nky = k*Ny;
                    nkx = k*Nx;
%                     tree(n).CLF.xe
%                     KY = nky
                    
%                     c=c+1;
%                     subplot(1,2,1)
%                     plot(x(1),x(2),S(n))
               

%                     plot(u(1,1),u(2,1),'r*')
%                     hold on
%                     
%                     plot(nky(1,1),nky(2,1),sty)
%                     hold on
%                     plot(nkx(1,1),nkx(2,1),'b*')
%                     hold on
%                     subplot(1,2,2)
%                    [0.55 0.55 0.55]
% 'Color',[0.55 0.55 0.55] 
                    quiver(x(1),x(2),u(1),u(2),'Color',[0.55 0.55 0.55] )
                    hold on
                end
            end
        end
    end
end
end