function ETriangles=edges2triangles(E)
E=E(:,1:2);
EWalks=walksRemoveRedundant(edges2walks(E,3,'symmetrize'));
flagTriangle=EWalks(:,1)==EWalks(:,end);
ETriangles=EWalks(flagTriangle,1:end-1);
ETriangles=walksRemoveRedundant(ETriangles,'flagRemoveReverse',true);
