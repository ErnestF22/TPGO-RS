%function plotfuntrials(f,N)
%Call the function f() N times and then plot the result
%f() can return a scalar or a vector
function rOut=funPlotTrials(f,N)
if ~exist('N','var')
    N=100;
end
r=f();
allr=zeros(length(r),N);
allr(:,1)=r;
for n=2:N
    allr(:,n)=f();
end
plot(allr')
if nargout>0
    rOut=allr;
end
disp('Maximum absolute value of f')
disp(max(abs(allr(:))))
disp('Minimum absolute value of f')
disp(min(abs(allr(:))))

