function bearingNetworkPlot_cameras(R,x)
NCameras=size(R,3);
for iCamera=1:NCameras
    draw2dcameraFromAxesAndCenter(R(:,:,iCamera),x(:,iCamera));
end
axis equal
