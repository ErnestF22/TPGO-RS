function d=edges2maxDegree(E)
A=edges2adjmatrix(E);
d=max(sum(A));
