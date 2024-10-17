%Apply one of the flip ambiguity to a point in the QREM
%function Q=essential_flipAmbiguity(Q,k)
%Input
%   k   Type of flipping (from 1 to 4)
%       1   R1=R1, R2=R2 (no flipping)
%       2   R1=Rxpi*R1, R2=Rxpi*Rzpi*R2
%       3   R1=R1, R2=Rzpi*R2
%       4   R1=Rxpi*R1, R2=Rxpi*R2
%       If omitted, a flip is picked at random
%Note that k=1,2 preserves the sign of the essential matrix, k=3,4 flips it.
function Q=essential_flipAmbiguity(Q,k)
if ~exist('k','var')
    k=round(3*rand+1);
end

Q(1:3,:)=essential_flipAmbiguity_R1(Q(1:3,:),k);
Q(4:6,:)=essential_flipAmbiguity_R2(Q(4:6,:),k);
