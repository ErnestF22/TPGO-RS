function v=rot3r3_vee(G,vHat)
v=[rot_vee(G2R(G),vHat(1:3,1:3,:)); squeeze(vHat(1:3,4,:))];
