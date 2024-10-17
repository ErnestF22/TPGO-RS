function [Rgt,Tgt] = generate_poses(nV)

Rgt = zeros(3,3,nV);
Tgt = zeros(3,nV);

for i=1:nV,   Rgt(:,:,i) = randrot(3);  end
for i=1:nV,   Tgt(:,i) = rand(3,1);  end


