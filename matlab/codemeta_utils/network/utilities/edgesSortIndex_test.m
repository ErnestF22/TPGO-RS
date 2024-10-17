function edgesSortIndex_test
A=adjgallery(7,'banded',2);
E=adj2edges(A,'undirected');

C=grCycleBasis(E);

[idxc,signc]=edgesSortIndex(E,C(:,2));

disp([E(idxc,:),signc])
