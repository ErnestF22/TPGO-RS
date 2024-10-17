% Dehomogenize a vector

function [x]=dehom(u)

   % Get the size of the vector
   dim = size(u,1);

   % Now, divide by the last element
   x = u(1:dim-1,:)./(ones(dim-1,1)*u(dim,:));



