function X=homographyContinuousTriangulate(x,nVec)
[lambda,xHom]=homographyContinuousTriangulateDepts(x,nVec);
X=([1;1;1]*lambda).*xHom;
