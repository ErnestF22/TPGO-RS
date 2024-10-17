function writeTestData
resetRands();

writeMatrixFile('matrix0.txt',round(10*randn(3,3,2)));

N=1000;
N2=100;

writeMatrixFile('dim.txt',N,' %d');
writeMatrixFile('dim2.txt',N2,' %d');

A=round(10*randn(3,3,N));
writeMatrixFile('matrix1.txt',A);
%display(A)

B=round(10*randn(3,3,N));
writeMatrixFile('matrix2.txt',B);

C=round(10*randn(3,3,N));
writeMatrixFile('matrix3.txt',C);

v=round(10*randn(3,N));
writeMatrixFile('vector1.txt',v);

R=fastMult3x3MatMat(A,B);
writeMatrixFile('MatMat.txt',R);

R=fastMult3x3HatMat(v,A);
writeMatrixFile('HatMat.txt',R);

R=fastMult3x3MatHat(A,v);
writeMatrixFile('MatHat.txt',R);

R=fastMult3x3MatMatMat(A,B,C);
writeMatrixFile('MatMatMat.txt',R);

R=fastMult3x3HatNormSqMat(v,A);
writeMatrixFile('HatNormSqMat.txt',R);

R=fastMult3x3MatHatNormSq(A,v);
writeMatrixFile('MatHatNormSq.txt',R);

%---%

A=round(10*randn(3,3,N2));
writeMatrixFile('matrixShort1.txt',A);
%display(A)

C=round(10*randn(3,3,N2));
writeMatrixFile('matrixShort2.txt',C);

e1=randint(1,N,[1 N2]);
writeMatrixFile('idx1.txt',e1);

e2=randint(1,N,[1 N2]);
writeMatrixFile('idx2.txt',e2);

R=fastMult3x3MatIdxTranspMatMatIdx(A,e1,B,C,e2);
writeMatrixFile('MatIdxTranspMatMatIdx.txt',R);

%display(R)

function writeMatrixFile(fileName,M,s)
if ~exist('s','var')
    s=' %+e';
end
fid=fopen(fullfile('../testData/',fileName),'w');
fprintf(fid, s, M(:));
fclose(fid);
