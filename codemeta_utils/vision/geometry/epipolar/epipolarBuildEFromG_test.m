function epipolarBuildEFromG_test
G1ref=RT2G(rot_randn(),randn(3,1));
G2ref=RT2G(rot_randn(),randn(3,1));

G1pose=invg(G1ref);
G2pose=invg(G2ref);

G12pose=computeRelativePoseFromG(G1pose,G2pose,'poses');
G12ref=computeRelativePoseFromG(G1ref,G2ref,'references');

E1=epipolarBuildEFromG(G12pose);
E2=epipolarBuildEFromG(G12ref);
E3=epipolarBuildEFromG(G1pose,G2pose,'poses');
E4=epipolarBuildEFromG(G1ref,G2ref,'references');

disp([E1 E1-E1;E2 E2-E1;E3 E3-E1; E4 E4-E1])
