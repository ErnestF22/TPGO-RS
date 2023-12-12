
function plot_controllers(A,b,k,k_added,y,e,flag,E,t)

X=0:e:100;
Y=0:e:100;
c=0;
for i=1:size(X,2)
    for j=1:size(Y,2)
        x = [X(i);Y(j)];
          if all(A*x<b)
            c=c+1;
            v=(y-x*ones(1,size(y,2)));
            if flag == 1
                v=t*(y-x*ones(1,size(y,2)));
            end
            v = reshape(v,[],1);
            u= k*v+k_added;
%             u = u/norm(u);
            quiver(x(1),x(2),u(1),u(2),'Color',[0.75 0.75 0.75]  )%
            hold on
          end
    end
end
end