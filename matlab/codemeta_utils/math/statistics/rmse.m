%Compute Root Mean Square Error (RMSE)
%function e=rmse(X1,X2)
%If X2 is provided, compute RMSE of X1-X2, otherwise, compute rmse of X1
function e=rmse(X1,X2)
if exist('X2','var')
    X1=X1-X2;
end
X1=reshape(X1,size(X1,1),[]);
e=sqrt(mean(sum(X1.^2),2));
