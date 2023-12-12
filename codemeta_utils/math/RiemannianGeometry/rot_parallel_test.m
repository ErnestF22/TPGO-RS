function rot_parallel_test
R=rot_exp(eye(3),rand*pi*rot_randTangentNormVector(eye(3)));
vatR=rot_randTangentNormVector(R);
S=rot_exp(R,rand*pi*vatR);

display('Geodesics parallel transport their own tangent vector')
[S*R'*vatR rot_parallel(R,S,vatR,'torotation')]

vatR2=rot_randTangentNormVector(R);
vatR2toS=rot_parallel(R,S,vatR2,'torotation');

sym=@(A) 0.5*(A+A');

display('Transported vector should be in tangent space at S')
sym(S'*vatR2toS)

display('Transitivity')
[vatR2 rot_parallel(S,R,vatR2toS,'torotation')]
