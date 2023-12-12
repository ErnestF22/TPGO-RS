close all
clear

%network:
[A1,A2,A3,b1,b2,b3]= network_example;


%% first region:
z1 = [1 1 1 1 1];
Z1=generate_bit_neighbors(z1);
%plot ingredients:
 R1 = plot_neighbor_reigons(Z1,A1,A2,A3,b1,b2,b3,1,0.1);
 
%find vertices:
[V1,T1,table1] = solving_dual_simplex(R1);

%% second region:
z2 = [1 1 0 1 1];
Z2=generate_bit_neighbors(z2);

%plot ingredients:
 R2 = plot_neighbor_reigons(Z2,A1,A2,A3,b1,b2,b3,2,0.1);
 
%find vertices:
[V2,T2,table2] = solving_dual_simplex(R2);
