function locSegSampleNetworkGallery_test
iSubplot=0;
for sampleNum=1:4
    for dimData=2:3
        resetRands()
        [E,x]=locSegSampleNetworkGallery(sampleNum,dimData);
        iSubplot=iSubplot+1;
        subplot(4,2,iSubplot)
        graphShow(E,x);
    end
end
