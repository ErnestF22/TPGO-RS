%Convert pairs of [3 x 3] matrices from stacked to side-by-side
function X=essential_col2row(Q)
X=[Q(1:3,:,:) Q(4:6,:,:)];

