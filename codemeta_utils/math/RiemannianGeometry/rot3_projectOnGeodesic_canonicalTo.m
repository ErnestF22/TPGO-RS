%Transform to canonical representation for projection on geodesic
%function [RCanonical,Re]=rot3_projectOnGeodesic_canonicalTo(R,Rz0,uz)
%Transform point to project to canonical representation
%(the geodesic is given by a rotation around z axis with origin at the
%identity)
function [RCanonical,Re]=rot3_projectOnGeodesic_canonicalTo(R,Rz0,uz)
Re=householderRotation(uz,3);
RCanonical=Re'*Rz0'*R*Re;
