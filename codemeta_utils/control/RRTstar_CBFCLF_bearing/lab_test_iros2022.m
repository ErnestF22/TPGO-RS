function tree = lab_test_iros2022
tree(1).position = [43;89];
tree(1).parent = [];
tree(1).L = [];

tree(2).position = [138.3;445];
tree(2).parent = 1;
tree(2).L = [22 23 36 35 33 27 28];

tree(3).position = [161.365;565.6];
tree(3).parent = 2;
tree(3).L = [23 22 6 11 10];

tree(4).position = [93.7;558];
tree(4).parent = 2;
tree(4).L = [10 6 22 23];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tree(5).position = [177;415];
tree(5).parent = 1;
tree(5).L = [22 23 21 33 35];

tree(6).position = [209;493];
tree(6).parent = 5;
tree(6).L = [6 11 22 23 12 21];

tree(7).position= [256;564];
tree(7).parent = 6;
tree(7).L = [9 34 6 11 20] ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tree(8).position = [156;350];
tree(8).parent = 1;
tree(8).L = [22 21 33 27 28];

tree(9).position = [185;414.3];
tree(9).parent= 8;
tree(9).L = [23 20 21 33];

tree(10).position = [286;528];
tree(10).parent = 9;
tree(10).L = [30 31 23 11 12 20];

tree(11).position = [198;391];
tree(11).parent = 8;
tree(11).L = [23 22 12 21 20] ;

tree(12).position = [280;476];
tree(12).parent= 11;
tree(12).L = [30 31 30 31 29] ;

tree(13).position = [325;384];
tree(13).parent = 11;
tree(13).L =[32 24 23 11 20] ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tree(14).position= [310;243];
tree(14).parent = 1;
tree(14).L = [2 39 35 37 27 28];

tree(15).position= [425;355];
tree(15).parent = 14;
tree(15).L = [32 24 2 1];

tree(16).position = [331;292];
tree(16).parent = 14;
tree(16).L = [24 1 2 38 37] ;

tree(17).position= [307;364];
tree(17).parent = 16;
tree(17).L = [24 32 1 2]  ;

tree(18).position = [425;432];
tree(18).parent= 16;
tree(18).L = [24 1 32 23] ;

tree(19).position = [442;553];
tree(19).parent= 18;
tree(19).L = [5 7 4 1];

end