%Convert pairs of [3 x 3] matrices from side-by-side to stacked
function Q=essential_row2col(X)
Q=[X(:,1:3,:); X(:,4:6,:)];