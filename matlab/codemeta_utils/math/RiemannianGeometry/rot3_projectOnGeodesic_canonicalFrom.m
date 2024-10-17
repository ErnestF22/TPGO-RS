%Transform from canonical representation for projection on geodesic
%function RProj=rot3_projectOnGeodesic_canonicalFrom(RProjCanonical,Rz0,Re)
%Transform a projected point from the canonical to the original
%representation
function RProj=rot3_projectOnGeodesic_canonicalFrom(RProjCanonical,Rz0,Re)
RProj=Rz0*Re*RProjCanonical*Re';