clear
close all
x = randi(20,2,10);
a = [1:10];
plot(x(1,:),x(2,:),'ko','MarkerSize',12)
% text(x(1,:)+0.5,x(2,:)+0.5,string(a(:)))
hold on
s.A=[];
s.b=[];
for i=1:size(x,2)
    s(i).A =[];
    s(i).b = [];
    for j=1:size(x,2)
        if i~=j
            if (i ==8) && (j==10)
                a=10;
            end
            [A,B] = bisector_segment(x(:,i),x(:,j));
            [Ax,bx] = make_set_convex(A,B,x(:,i));
            s(i).A = [s(i).A; Ax];
            s(i).b = [s(i).b; bx];
        end
    end
end
plot_bisector_segment(s)