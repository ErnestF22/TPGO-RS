function Asym = symm(A)
%SYMM Return symmetric part of A

d3 = size(A,3);
if (d3>1)
    Asym = zeros(size(A));
    for ii = 1:d3
        Asym(:,:,ii) = 0.5 .* (A(:,:,ii) + A(:,:,ii)')
    end
else
    Asym = 0.5 .* (A + A');
end


end

