%Displays a matrix in row blocks separated by a blank line
%function displayMatrixInRowBlocks(M,sizeR)
function displayMatrixInRowBlocks(M,sizeR)
for iR=1:sizeR:size(M,1)
    disp(num2str(M(iR:iR+sizeR-1,:),' %+.3f'))
    fprintf('\n');
end
