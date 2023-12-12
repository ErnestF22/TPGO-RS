clear all
close all



%location of vertices:
y = [1.7 1;1.5 8;11.5 9;13.4 2.5]'; %general -1<K<1 and E= [1;-0.1]; and
% Cl= 0.5 and Cb = 1;
% y = [2 2;1.5 6.5;11.5 16.5;11 11]'; %parallel
% y = [2 3.5;1.5 8.5;10 13;9 6]'; %open set
L = [12 8;13 5;16 10;18 1;16 10;18 1]' ;
% L = [12 5;13 5]'+[5;0];%16 10;18 1]'%+[0 7;0 7;0 5;0 5]' ;
L = y;


p = (y(:,2)+y(:,4))/2;
[A,b] = convexSet(y);
[Ax,bx] = LPSet(A,b,p);

Ah = -Ax;
bh = bx;


Ah(3,:)=[];
bh(3)=[];



exitDir =  -Ax(3,:)';
% exitDir = exitDir/norm(exitDir);

% xe = intersectionTwoLines(A(2,1),b(2),A(4,1),b(4));
xe = (y(:,3)+y(:,4))/2;
% xe = [5.5;5.5];

epsilon = 0.5;
Wb = 1*ones(size(Ah,1),1);
Wl = 1;

E = [1;-0.1];
% E = E/norm(E);
s = 1.6;
flag_trons_test = 0;
% 0: translation
% 1: bearing
flag_controller = 1;
flag_navigation = 1;

if flag_controller == 0% for in and KY cl=1 and cb=0.5
    Cl = 1;
    Cb = 1;
else
    % tuned by tron
    Cl = 0.01;
    Cb = 1000;
end

[K,k_added] = find_controller(exitDir,Wb,Wl,Ah,Ax,bx,bh,xe,y,L,Cl,Cb,E,flag_controller,s);
display(K)
display(k_added)
figure(2)
plot_controllers(Ax,bx,K,k_added,L,epsilon,flag_navigation,E,1/s)

plot( [y(1,end) y(1,1)],[y(2,end) y(2,1)], '.-r');
hold on
plot( [y(1,2) y(1,1)],[y(2,2) y(2,1)], '.-r');
hold on
plot( [y(1,2) y(1,3)],[y(2,2) y(2,3)], '.-r');
hold on
plot(xe(1),xe(2),'c*')


for i=1:size(y,2)
    plot(y(1,i),y(2,i),'b*')
    hold on
end

for i=1:size(L,2)
    plot(L(1,i),L(2,i),'g*')
    hold on
end

grid on


if flag_trons_test
    %exit direction
    quiver(11.5,5,exitDir(1),exitDir(2),'color',[0,0,1,0.5],'linewidth',1)
    
  %polygons of u for different S
    for i_corner=1:4
        xBase=y(:,i_corner);
        uBase=-K*vec(xBase-L)+k_added;
        sMin= 0.89
        sMax= 1.25
        uBaseS=@(s) -K*s*vec(xBase-L)+k_added;
        quiver(xBase(1),xBase(2),uBase(1),uBase(2),'b')
        uPolygon=[xBase, xBase+uBaseS(sMin), xBase+uBaseS(sMax),xBase];
        plot(uPolygon(1,:),uPolygon(2,:),'k','linewidth',1)
    end
    axis equal
    axis([-5 25 0 10])
    save('controller_data','K','k_added','sMin','sMax')
    
end