%Composes two rigid body motions in a numerically stable manner
function G=composeG(G,G1)
G=0.5*(G*G1+G/invg(G1));
