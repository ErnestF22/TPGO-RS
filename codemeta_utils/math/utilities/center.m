function xCentered=center(x)
warning('center is deprecated. Use ccenter instead.')
xCentered=x-mean(x,2)*ones(1,size(x,2));
