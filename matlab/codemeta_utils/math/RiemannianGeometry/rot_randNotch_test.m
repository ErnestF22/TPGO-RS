function rot_randNotch_test
v=[45 90]*pi/180;
R=rot_randn([],[],300);
RNoise=rot_randNotch(R,v,[],'U',diag([1 1 0]));
d=rot_dist(R,RNoise,'vector');

cumDistPerc(d*180/pi)



