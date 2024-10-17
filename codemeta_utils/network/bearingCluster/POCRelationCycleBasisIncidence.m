function POCRelationCycleBasisIncidence
A=adjgallery(7,'banded',2);
[B,E]=adj2incmatrix(A,'oriented');
C=grOrientCycleBasis(grCycleBasis(E),E)';

PB=orthComplementProjector(B);
PC=C'*((C*C')\C);

disp(PB-PC)
