function findComponents_test
%a graph with two components: [1 2 3] and [4 5 6 7]
E=[1 2 3 5 6 7 8;
   2 3 1 6 7 8 5]';

C=findComponents(E);

disp(C)

