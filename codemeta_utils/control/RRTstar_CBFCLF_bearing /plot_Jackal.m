figure(1)
I = imread('Lab_Cspace.JPG');
h = image([0 100],[0 100],I); 
uistack(h,'bottom')
hold on
axis equal
axis ([0 100 0 100])
plot(0,0,'r*')
