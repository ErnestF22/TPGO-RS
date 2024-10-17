function idxValid=convValidIdx(ny,nh,shape)
switch shape
    case 'full'
        idxValid=nh:ny;
    case 'same'
        idxValid=floor((nh+1)/2):(ny-floor(nh/2));
    case 'valid'
        idxValid=1:(ny-nh+1);
    otherwise
        error('Shape argument not valid.')
end
