%this function load the Tree data and then deform the nevironment to test
%the contollers on the deformed environment. 
clear all
w = load('solved_tree.mat');
T = w.T;
L =[];
start = [45;60];
obstacle_list= [30 30 12;
    30 50 12;
    30 70 12;
    30 80 12;
    50 30 12;
    70 30 12;
    80 30 12]';

% 
l1 = [42 80;100 90;42 42;80 42]';
l2 = [0 0;30 18;80 18;80 0]';
l3 = [0 0;30 18;80 18;80 0]';
l4 = [0 0;30 18;18 80;0 90;18 30]';


figure(1)
for i=1:size(obstacle_list,2)
    plot_circle(obstacle_list(1,i),obstacle_list(2,i)...
        ,obstacle_list(3,i),'k')
end
plot(l1(1,:),l1(2,:),'b*','MarkerSize',5,'LineWidth',2)
hold on
plot(l2(1,:),l2(2,:),'b*','MarkerSize',5,'LineWidth',2)
hold on
plot(l3(1,:),l3(2,:),'b*','MarkerSize',5,'LineWidth',2)
hold on
plot(l4(1,:),l4(2,:),'b*','MarkerSize',5,'LineWidth',2)
axis([0 100 0 100])
% 
navigate(start,T,obstacle_list,1,'r*');
% axis([0 100 0 100])
% plotControllers(T,'g*')

obstacle_list2= [32 30 12;
    36 43 12;
    39 58 12;
    44 78 12;
    50 30 12;
    70 30 12;
    80 30 12]';
grid on
figure(1)
for i=1:size(obstacle_list2,2)
    plot_circle(obstacle_list2(1,i),obstacle_list2(2,i)...
        ,obstacle_list2(3,i),'r')
end

l1 = [56 78;114 88;52 47;80 42]';
l2 = [0 0;32 18;80 18;80 0]';
l3 = [0 0;32 18;80 18;80 0]';
l4 = [0 0;32 18;32 78;0 82;20 30]';

T = fitLandmarks(T,l1,l2,l3,l4);

plot(l1(1,:),l1(2,:),'g*','MarkerSize',5,'LineWidth',2)
hold on
plot(l2(1,:),l2(2,:),'g*','MarkerSize',5,'LineWidth',2)
hold on
plot(l3(1,:),l3(2,:),'g*','MarkerSize',5,'LineWidth',2)
hold on
plot(l4(1,:),l4(2,:),'g*','MarkerSize',5,'LineWidth',2)

axis([0 100 0 100])

figure(1)
% plotControllers(T,'r*')
navigate(start,T,obstacle_list,1,'b*');
axis([-5 105 -5 105])

grid on
