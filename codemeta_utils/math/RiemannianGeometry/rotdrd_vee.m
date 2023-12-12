function v=rotdrd_vee(G,vHat)
d=size(g,1);
v=[rot_vee(G2R(G),vHat(1:d,1:d,:)); squeeze(vHat(1:d,d+1,:))];
