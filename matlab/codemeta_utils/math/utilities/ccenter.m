function xCentered=ccenter(x)
xCentered=x-mean(x,2)*ones(1,size(x,2));
