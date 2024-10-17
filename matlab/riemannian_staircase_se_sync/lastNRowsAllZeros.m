function out = lastNRowsAllZeros(inputArg1, d, N)
%LASTNROWSALLZEROS Return true if for every ii, inputArg1(:,:,ii) last rows
% added by Riemannian Staircase are equal to 0

R_lastRowsAllZeros = matStack(any(multitransp(inputArg1)));
num_rows_to_check = size(inputArg1, 1) - d;

out = true;
for ii = 1:N
    for jj = 0:num_rows_to_check-1
        out_ij = boolean(~any(R_lastRowsAllZeros(:,end-jj)));
        if (out_ij == false)
            out = false;
            break;
        end
    end
end

end %function

