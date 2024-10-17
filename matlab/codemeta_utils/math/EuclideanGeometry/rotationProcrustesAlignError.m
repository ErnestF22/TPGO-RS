function d=rotationProcrustesAlignError(R1,R2,varargin)
[~,R2Transformed]=rotationProcrustes(R1,R2,varargin{:});
d=rot_dist(R1,R2Transformed,'vector');

