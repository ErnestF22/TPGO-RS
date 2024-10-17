N=6;
%A=ones(N)-eye(N);
A=adjgallery(N,'banded',2);

figure(1)
gshow(A,'adj');
title(['Network diameter = ' num2str(findNetworkDiameter(A))])

