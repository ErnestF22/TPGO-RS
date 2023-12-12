function [M] = extractComp(M,rowStart,rowEnd,colStart,colEnd)
% General function to extract components of a matrix. Useful for extracting
% components of matrix from function handles
M = M(rowStart:rowEnd,colStart:colEnd);
end