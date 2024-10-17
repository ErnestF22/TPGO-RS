function G=rotLocFrobCostPairMatrix(RRel)
I=eye(size(RRel,1));
G=[I -RRel; -RRel' I];
