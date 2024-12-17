p = 4;
d = 3;
n = 1;

% polar
Mpolar = stiefelfactory(p,d,n);
Mpolar.retr = Mpolar.retr_polar;

% QR
Mqr = stiefelfactory(p,d,n);


thr = 1e-4;

for ii=1:1000
    xi = make_rand_stiefel_3d_array(p,d,n);
    vi = stiefel_randTangentNormVector(xi);
    xpolar = Mpolar.retr(xi, vi);
    xqr = Mqr.retr(xi, vi);

    disp("ii")
    disp(ii)

    if max(abs(xpolar- xqr), [], "all") > thr
        disp("QR and POLAR retractions produce different outputs")
        disp("[xpolar, xqr]")
        disp([xpolar, xqr])
        break;
    end
end