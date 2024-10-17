function localization_display(GiTruth,GiInit,GiEst,A)
optsDisplay={'references','adjmatrix',A};
subplot(1,2,1)
testNetworkDisplay(GiInit,optsDisplay{:},'Estimated')
hold on
testNetworkDisplay(GiTruth,optsDisplay{:})
hold off
title('Init')
subplot(1,2,2)
testNetworkDisplay(GiEst,optsDisplay{:},'Estimated')
hold on
testNetworkDisplay(GiTruth,optsDisplay{:})
hold off
