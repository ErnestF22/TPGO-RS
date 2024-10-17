function AProj=rangeProjection(A,B)
B=orth(B);
AProj=orth(A-B*(B'*A));
