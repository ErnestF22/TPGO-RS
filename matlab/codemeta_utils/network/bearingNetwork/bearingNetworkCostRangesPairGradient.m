function gradPhi=bearingNetworkCostRangesPairGradient(y,yg,ny,nyg,funs)
gradPhiij=bearingCostGeneralRangesGradient(y,yg,ny,nyg,funs);
gradPhiji=bearingCostGeneralRangesGradient(-y,-yg,ny,nyg,funs);
gradPhi=[gradPhiij gradPhiji];
