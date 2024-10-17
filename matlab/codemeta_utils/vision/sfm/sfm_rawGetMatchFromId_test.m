function sfm_rawGetMatchFromId_test
m=[1:4;2:5];

disp([1 2; 1 getMatchFromId(m,1)])
disp([3 4; 3 getMatchFromId(m,3)])
disp('This should be empty')
disp(['[' num2str(getMatchFromId(m,5)) ']'])
