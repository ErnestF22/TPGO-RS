%function B=fastMultHatMat(v,A)
%Same as fastMultMatHat, but with the order of the arguments reversed
%
%See also fastMultMatHat
function B=fastMult3x3HatMat(v,A)

B=[...
    v(3:3:end).*A(2:9:end)-v(2:3:end).*A(3:9:end);
    -v(3:3:end).*A(1:9:end)+v(1:3:end).*A(3:9:end);
    v(2:3:end).*A(1:9:end)-v(1:3:end).*A(2:9:end);
    v(3:3:end).*A(5:9:end)-v(2:3:end).*A(6:9:end);
    -v(3:3:end).*A(4:9:end)+v(1:3:end).*A(6:9:end);
    v(2:3:end).*A(4:9:end)-v(1:3:end).*A(5:9:end);
    v(3:3:end).*A(8:9:end)-v(2:3:end).*A(9:9:end);
    -v(3:3:end).*A(7:9:end)+v(1:3:end).*A(9:9:end);
    v(2:3:end).*A(7:9:end)-v(1:3:end).*A(8:9:end)
    ];
