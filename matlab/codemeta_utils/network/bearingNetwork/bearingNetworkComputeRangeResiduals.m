function q=bearingNetworkComputeRangeResiduals(y,yg,ny,nyg)
q=ny.*bearingNetworkComputeBearingsCosines(y,yg)-nyg;
