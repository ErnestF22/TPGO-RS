function [tree,obs1,obs2,l1,l2,l3,l4] = experiment_tree_global
tree(1).position = [364.54;-52];
tree(1).parent = [];
tree(1).L = [];

tree(2).position = [0;-75];
tree(2).parent = 1;
tree(2).L = 1;

tree(3).position = [132;-115];
tree(3).parent = 1;
tree(3).L = 2;

tree(4).position = [128;-200];
tree(4).parent = 3;
tree(4).L = 2;

tree(5).position = [50;-240];
tree(5).parent = 3;
tree(5).L = 2;

tree(6).position = [52;-304];
tree(6).parent = 5;
tree(6).L = 2;

tree(7).position = [133;-310];
tree(7).parent = 6;
tree(7).L = 3;

tree(8).position = [51;-350];
tree(8).parent = 6;
tree(8).L = 3;

tree(9).position = [102;-400];
tree(9).parent = 11;
tree(9).L = 3;

tree(10).position = [98;-475];
tree(10).parent = 11;
tree(10).L = 3;

tree(11).position = [160;-450];
tree(11).parent = 12;
tree(11).L = 1;
tree(11).L = 4;

tree(12).position = [290;-410];
tree(12).parent = 13;
tree(12).L = 1;
tree(12).L = 4;

tree(13).position = [330;-350];
tree(13).parent = 1;
tree(13).L = 1;

tree(14).position = [370;-500];
tree(14).parent = 1;
tree(14).L = 4;

obs1 = 100*[1.55	-1.36;
            1.85	-1.36;
            1.85	-3.36;
            2.78	-3.36;
            2.78	-3.66;
            1.85	-3.66;
            1.85	-3.95;
            1.55	-3.95;
            1.55	-2.84;
            0.9	    -2.84;
            0.9	    -2.4;
            1.55	-2.4]';
r=0.15;        
obs2 = 100*[1.55-r,	-1.36+r;
            1.85+r,	-1.36+r;
            1.85+r,	-3.36+r;
            2.78+r	,-3.36+r;
            2.78+r,	-3.66-r;
            1.85+r,	-3.66-r;
            1.85+r,	-3.95-r;
            1.55-r,	-3.95-r;
            1.55-r,	-2.84-r;
            0.9-r ,   -2.84-r;
            0.9-r ,   -2.4+r;
            1.55-r,	-2.4+r]';


l1=100*[4.187 0.165;1.8699	-1.923;1.8749	-3.023;4.186	-3.441]';
l2=100*[1.487	-0.227;1.522	-2.05;1.048	-2.375;-0.199	-0.524]';
l3=100*[-0.098	-2.951;0.6752	-4.5965;1.697	-3.979;2.438	-5.463]';
l4=100*[3.696	-4.456;2.7047	-3.6533;2.0606	-3.662]';
end
