function POCRotateToMinimizeLastEntries
dLow=3;
d=5;
nbCols=10;

% Generate low-rank data
xLow=randn(dLow,nbCols);
xTrue=[xLow; zeros(d-dLow,nbCols)];
QTrue=rot_randn(eye(d));
x=QTrue*xTrue;



% Recover a Q for aligning the data
[Q,~,~]=svd(x*x');
disp('x=')
disp(x)
disp('Q''*x=')
disp(Q'*x)


