%Linear geodesic interpolation between rotations
%function RQuery=rot_interpolationLinearPair(t,R,tQuery)
%Returns the interpolated rotation obtained by using the geodesic between
%the two given rotations. It uses the same mechanisms also for future/past
%extrapolation.
%Inputs
%   t       [1 x 2] array of time instants
%   R       [3 x 3 x 2] array of rotations corresponding to t
%   tQuery  a scalar time instant
function RQuery=rot_interpolationLinearPair(t,R,tQuery)
if length(t)~=2 || size(R,3)~=2
    error('This funciton is restricted to computing the linear interpolation between pairs');
end
NTQuery=length(tQuery);
if NTQuery>1
    RQuery=zeros(3,3,NTQuery);
    for it=1:NTQuery
        RQuery(:,:,it)=rot_interpolationLinearPair(t,R,tQuery(it));
    end
else
    v=rot_log(R(:,:,1),R(:,:,2));
    if t(2)==t(1)
        tRatio=0;
    else
        tRatio=(tQuery-t(1))/(t(2)-t(1));
    end
    RQuery=rot_exp(R(:,:,1),tRatio*v);
end