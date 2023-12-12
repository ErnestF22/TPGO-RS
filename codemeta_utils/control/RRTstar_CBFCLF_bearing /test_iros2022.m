function [tree] = test_iros2022
tree(1).position = [10;10];
tree(1).parent = [];
tree(1).L = 3;

tree(2).position = [12.4;9.6];
tree(2).parent = 1;
tree(2).L = 4;

tree(3).position = [17.9;36.8];
tree(3).parent = 1;
tree(3).L = 4;

tree(4).position = [19.7;46];
tree(4).parent = 1;
tree(4).L = 4;

tree(5).position = [26.9;65.67];
tree(5).parent = 1;
tree(5).L = 4;

tree(6).position = [16.88;43.62];
tree(6).parent = 1;
tree(6).L = 4;

tree(7).position = [51.71;38];
tree(7).parent = 20;
tree(7).L = 3;

tree(8).position = [75.5;55.3];
tree(8).parent = 19;
tree(8).L = 1;

tree(9).position = [76.45;82.3];
tree(9).parent = 8;
tree(9).L = 1;

tree(10).position = [52.3;15.7];
tree(10).parent = 19;
tree(10).L = 3;

tree(11).position = [97.4;62.55];
tree(11).parent = 10;
tree(11).L = 1;

tree(12).position = [73.83;24];
tree(12).parent = 10;
tree(12).L = 3;

tree(13).position = [33.4;43.8];
tree(13).parent = 3;
tree(13).L = 4;

tree(14).position = [41.5;41.6];
tree(14).parent = 13;
tree(14).L = 4;

tree(15).position = [31.88;56.83];
tree(15).parent = 4;
tree(15).L = 4;

tree(16).position = [5.41;97.96];
tree(16).parent = 5;
tree(16).L = 2;

tree(17).position = [21.6;63.4];
tree(17).parent = 6;
tree(17).L = 2;

tree(18).position = [21.15;66.7];
tree(18).parent = 17;
tree(18).L = 2;

tree(19).position = [47.9;14.16];
tree(19).parent = 20;
tree(19).L = 3;

tree(20).position = [42.5;11];
tree(20).parent = 1;
tree(20).L = 3;
end