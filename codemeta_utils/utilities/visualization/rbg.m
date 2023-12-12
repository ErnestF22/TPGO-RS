%RGB colormap
%function c=rbg(N)
%Colormap that goes from red to blue to green. Useful for plots, as it
%excludes yellow.
function c=rbg(N)
u=ones(N,1);
c=hsv2rgb([mod(1-linspace(0,0.6249,N)'-0.0211,1) u u]);
