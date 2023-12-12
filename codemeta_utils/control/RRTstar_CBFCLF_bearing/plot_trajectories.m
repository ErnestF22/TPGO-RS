% clear 
% close all
X0 = 100*csvread('trajectory_sequence_0.csv',1,0);
% X0 = X0(1:end-3,:);
X1 = 100*csvread('trajectory_sequence_1.csv',1,0);
% X1 = X1(1:340,:);
X2 = 100*csvread('trajectory_sequence_2.csv',1,0);
% X2 = X2(21:end,:);
X3 = 100*csvread('trajectory_sequence_3.csv',1,0);

[T,obs,obs_offset,~,~,~,~] = experiment_tree_global;
s = [363,-51];
figure(2)
ll= deformed_l3;
ll(:,1) = ll(:,1)+[-30;-60];
plot_exp(deformed_l1,obsd_deformed,obs_offset2,deformed_l1,deformed_l2,ll,deformed_l4)
Ss = ["r","m","g","b","k"];
Sp = ["r.-","m.-","g.-","b.-","k.-"];
start = X0(1,:);
i = 1;
plot(start(1,1),start(1,2),'p','MarkerSize',6,'MarkerFaceColor',Ss(i),'LineWidth',1.5)
hold on
x = X0;
plot(x(:,1),x(:,2),Sp(i))

start = X1(1,:);
i = 2;
plot(start(1,1),start(1,2),'p','MarkerSize',6,'MarkerFaceColor',Ss(i),'LineWidth',1.5)
hold on
x = X1;
plot(x(:,1),x(:,2),Sp(i))

start = X2(1,:);
i = 3;
plot(start(1,1),start(1,2),'p','MarkerSize',6,'MarkerFaceColor',Ss(i),'LineWidth',1.5)
hold on
x = X2;
plot(x(:,1),x(:,2),Sp(i))

start = X3(1,:);
i = 4;
plot(start(1,1),start(1,2),'p','MarkerSize',6,'MarkerFaceColor',Ss(i),'LineWidth',1.5)
hold on
x = X3;
plot(x(:,1),x(:,2),Sp(i))

plot(s(1),s(2),'o','MarkerSize',6,'MarkerEdgeColor',[0 0.5 0.5],'MarkerFaceColor',[0 0.6 0.6],'LineWidth',1.5)
axis equal
grid on 
xlabel('x')
ylabel('y')

