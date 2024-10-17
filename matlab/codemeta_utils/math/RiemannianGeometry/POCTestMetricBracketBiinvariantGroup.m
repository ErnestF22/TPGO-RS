function POCTestMetricBracketBiinvariantGroup
R=eye(5);
X=rot_randTangentNormVector(R);
Y=rot_randTangentNormVector(R);
Z=rot_randTangentNormVector(R);

brXY=rot_bracket(R,X,Y);
brYZ=rot_bracket(R,Y,Z);

disp([rot_metric(R,brXY,Z) rot_metric(R,X,brYZ)])
