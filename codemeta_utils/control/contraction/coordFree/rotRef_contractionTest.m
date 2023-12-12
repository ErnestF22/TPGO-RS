% rotRef_contractionTest.m
% Try to find some gains [kd;kv;kp] such that the augmented system with a
% reference trajectory convergences at some minimum exponential rate 
% beta>0 with metric m.

%% Init
% Clean work space
close all; clear all; clc;

% % Create all the constraining Ai matrices in trace(Ai*X) <= 0 
% syms kd kv kp beta m1 m2 m3 m4 m5 m6 mag_R mag_RRef mag_W real;
% A_all = rotRef_allBounds_Matrix(mag_R,mag_RRef,kd,kv,kp,mag_W,beta);
% 
% % Convert A_all to function handles
% A_all_func = cell(size(A_all));
% for ii = 1:length(A_all)
%     A_all_func{ii} = matlabFunction(A_all{ii},'Vars',[mag_R,mag_RRef,kd,kv,kp,mag_W,beta,m4]);
% end

% Default Parameters
maxIter = 5;
beta_LB = -10;
beta_UB = 10;
% init_gains = [100;80;1]; %[kd;kv;kp]; %[100;80;1] works with mag_W = 1
% kdList=linspace(1,100,10);
% [A,B]=meshgrid(kdList,kdList);
% c=cat(2,A',B');
% d=reshape(c,[],2);
% init_gains=[d';ones(1,100)];
kd_list = 100; % Define all kd gains to test
kv_list = 80; % Define all kv gains to test
kref_list = [.01, 5, 10, 15, 20];
% Generate all possible gains
init_gains = allcomb(kd_list,kv_list,kref_list)';
mag_W = 1; % Maximum angular speed
% init_gains = 100*rand(3,1);
% mag_W = 2;
mag_R = pi/4; % Maximum distance error from R to RRef, should be equal to norm(log(R,RRef),2)
mag_RRef = pi/4; % Maximum distance error from RRef to I, should be equal to norm(log(RRef,I),2)
% NOTE: {mag_R,mag_RRef} should cover the total space on SO3.

%% Batch jobs
% create a batch jobs to find local optimal gains for each initial gain
% (determine by the number of columns of init_gains)
for ii = 1:size(init_gains,2)
    tempJob = batch(@rotRef_gainSearch, 4, {mag_R,mag_RRef,...
        mag_W,...
        'gainsearch_maxiter',maxIter,...
        'init_gains',init_gains(:,ii),...
        'beta_ub',beta_UB,...
        'beta_lb',beta_LB,...
        'm4',1,...
        'maxgain',110});
end

%% Test (Comment for batch jobs)
% This is sample code to use rotRef_gainSearch incase you wanted to see the
% progress
for ii = 1:size(init_gains,2)
    [kd,kv,kp,beta] = rotRef_gainSearch(mag_R,mag_RRef,...
            mag_W,...
            'gainsearch_maxiter',maxIter,...
            'init_gains',init_gains(:,ii),...
            'beta_ub',beta_UB,...
            'beta_lb',beta_LB,...
            'm4',1);
end

% One of the opt conditions should be satisfied
% extraData = cell(4,1);
[M_nn,f,X,cvx_optval,set_A] = rotRef_contractionOpt(mag_R,mag_RRef,kd,kv,kp,mag_W,beta,'m4',1,'schur_comp',A_all_func);
w = cnormalize(randn(3,1))*mag_W;
M = rotRef_contractionMat2(beta,kd,kv,kp,M_nn);
RRef = rot_exp(eye(3),hat3(cnormalize(randn(3,1))*mag_RRef));
R = rot_exp(RRef,RRef*hat3(cnormalize(randn(3,1))*mag_R));
M_eval = M(R,w,RRef); M_symm = (M_eval + M_eval')/2;
fprintf('Contraction matrix eigenvalues:\n')
eig(M_symm)

% Check gershgorin disc bounds
[max_eVals] = rotRef_contractionMatrix_greshBound(mag_R,mag_RRef,kd,kv,kp,mag_W,beta,M_nn)

%% Function taken from https://www.mathworks.com/matlabcentral/fileexchange/10064-allcomb-varargin
function A = allcomb(varargin)
% ALLCOMB - All combinations
%    B = ALLCOMB(A1,A2,A3,...,AN) returns all combinations of the elements
%    in the arrays A1, A2, ..., and AN. B is P-by-N matrix where P is the product
%    of the number of elements of the N inputs. 
%    This functionality is also known as the Cartesian Product. The
%    arguments can be numerical and/or characters, or they can be cell arrays.
%
%    Examples:
%       allcomb([1 3 5],[-3 8],[0 1]) % numerical input:
%       % -> [ 1  -3   0
%       %      1  -3   1
%       %      1   8   0
%       %        ...
%       %      5  -3   1
%       %      5   8   1 ] ; % a 12-by-3 array
%
%       allcomb('abc','XY') % character arrays
%       % -> [ aX ; aY ; bX ; bY ; cX ; cY] % a 6-by-2 character array
%
%       allcomb('xy',[65 66]) % a combination -> character output
%       % -> ['xA' ; 'xB' ; 'yA' ; 'yB'] % a 4-by-2 character array
%
%       allcomb({'hello','Bye'},{'Joe', 10:12},{99999 []}) % all cell arrays
%       % -> {  'hello'  'Joe'        [99999]
%       %       'hello'  'Joe'             []
%       %       'hello'  [1x3 double] [99999]
%       %       'hello'  [1x3 double]      []
%       %       'Bye'    'Joe'        [99999]
%       %       'Bye'    'Joe'             []
%       %       'Bye'    [1x3 double] [99999]
%       %       'Bye'    [1x3 double]      [] } ; % a 8-by-3 cell array
%
%    ALLCOMB(..., 'matlab') causes the first column to change fastest which
%    is consistent with matlab indexing. Example: 
%      allcomb(1:2,3:4,5:6,'matlab') 
%      % -> [ 1 3 5 ; 1 4 5 ; 1 3 6 ; ... ; 2 4 6 ]
%
%    If one of the N arguments is empty, ALLCOMB returns a 0-by-N empty array.
%    
%    See also NCHOOSEK, PERMS, NDGRID
%         and NCHOOSE, COMBN, KTHCOMBN (Matlab Central FEX)
% Tested in Matlab R2015a and up
% version 4.2 (apr 2018)
% (c) Jos van der Geest
% email: samelinoa@gmail.com
% History
% 1.1 (feb 2006), removed minor bug when entering empty cell arrays;
%     added option to let the first input run fastest (suggestion by JD)
% 1.2 (jan 2010), using ii as an index on the left-hand for the multiple
%     output by NDGRID. Thanks to Jan Simon, for showing this little trick
% 2.0 (dec 2010). Bruno Luong convinced me that an empty input should
% return an empty output.
% 2.1 (feb 2011). A cell as input argument caused the check on the last
%      argument (specifying the order) to crash.
% 2.2 (jan 2012). removed a superfluous line of code (ischar(..))
% 3.0 (may 2012) removed check for doubles so character arrays are accepted
% 4.0 (feb 2014) added support for cell arrays
% 4.1 (feb 2016) fixed error for cell array input with last argument being
%     'matlab'. Thanks to Richard for pointing this out.
% 4.2 (apr 2018) fixed some grammar mistakes in the help and comments
narginchk(1,Inf) ;
NC = nargin ;
% check if we should flip the order
if ischar(varargin{end}) && (strcmpi(varargin{end}, 'matlab') || strcmpi(varargin{end}, 'john'))
    % based on a suggestion by JD on the FEX
    NC = NC-1 ;
    ii = 1:NC ; % now first argument will change fastest
else
    % default: enter arguments backwards, so last one (AN) is changing fastest
    ii = NC:-1:1 ;
end
args = varargin(1:NC) ;
if any(cellfun('isempty', args)) % check for empty inputs
    warning('ALLCOMB:EmptyInput','One of more empty inputs result in an empty output.') ;
    A = zeros(0, NC) ;
elseif NC == 0 % no inputs
    A = zeros(0,0) ; 
elseif NC == 1 % a single input, nothing to combine
    A = args{1}(:) ; 
else
    isCellInput = cellfun(@iscell, args) ;
    if any(isCellInput)
        if ~all(isCellInput)
            error('ALLCOMB:InvalidCellInput', ...
                'For cell input, all arguments should be cell arrays.') ;
        end
        % for cell input, we use to indices to get all combinations
        ix = cellfun(@(c) 1:numel(c), args, 'un', 0) ;
        
        % flip using ii if last column is changing fastest
        [ix{ii}] = ndgrid(ix{ii}) ;
        
        A = cell(numel(ix{1}), NC) ; % pre-allocate the output
        for k = 1:NC
            % combine
            A(:,k) = reshape(args{k}(ix{k}), [], 1) ;
        end
    else
        % non-cell input, assuming all numerical values or strings
        % flip using ii if last column is changing fastest
        [A{ii}] = ndgrid(args{ii}) ;
        % concatenate
        A = reshape(cat(NC+1,A{:}), [], NC) ;
    end
end
end
