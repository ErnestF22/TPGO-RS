function sfm_rawIntersectMatch_test
m1=[1:4;2:5];
m2=fliplr([3:8; 4:9]);

mTruth=[3:4;4:5];
disp([mTruth; intersectMatch(m1,m2)])

m3=[8:10; 9:11];
disp('This should be empty')
disp(['[' num2str(intersectMatch(m1,m3)) ']'])
