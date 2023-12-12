function Ab_set=ReluPlexNetwork(L)
A1 = [2 3;-3 8];
b1 = [-3;2];

A2 = [5 -1;3 -2];
b2 = [-10;21];

A3 = [-3 5];
b3 = 35;

A4 = -1;
b4 = 100;

Ab_set = struct;
for i=1:L
    Ab_set(i).A = eval(['A' num2str(i)]);
    Ab_set(i).b = eval(['b' num2str(i)]);
end

end
% t1 = -30;
% t2 = 45;
% 
% A1 = [cosd(t1) -sind(t1);sind(t1) cosd(t1)];
% b1 = [1;1];
% 
% A2 = [cosd(t2) -sind(t2);sind(t2) cosd(t2)];
% b2 = [1;1];
% 
% A3 = [1 1];
% b3 = 1;