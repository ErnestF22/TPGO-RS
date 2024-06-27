function B = tail_array(A, k)
%TAIL Get last rows of array, table, or timetable
%   TAIL_ARRAY(A) displays the last eight rows of the array or table A in the
%   command window without storing a value.
%
%   TAIL_ARRAY(A,K) displays up to K rows from the end of the A. If A contains
%   fewer than K rows, then the entire array or table is displayed.
%
%   B = TAIL_ARRAY(A) or B = TAIL_ARRAY(A,K) returns the last eight rows, or up to K
%   rows, of the array or table A.
%
%   See also HEAD, TOPKROWS, SIZE

%   Copyright 2016-2022 The MathWorks, Inc.

h = size(A,1);
if k < 1
    error('tail_array() needs a vaild number of rows (integer > 0)');
end

if k > h 
    disp('WARNING: trying to get more tail rows than matrix has!');
    k = h;
end

if nargin<2
    k = min([8, h]);
    B = A(h-k+1:h,:);
else
    B = A(h-k+1:h,:);
end

end %tile function