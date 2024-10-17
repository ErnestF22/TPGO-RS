function s=POCSubspaceAngles(A,B)
s=svd(orth(A)'*orth(B));
