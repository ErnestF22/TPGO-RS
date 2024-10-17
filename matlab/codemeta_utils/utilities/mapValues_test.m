function mapValues_test
v=[1 1 2 2 4 4 5];
m=[1 2 4; 2 1 3]';

vMapped=mapValues(v,m);
disp([v;vMapped])
