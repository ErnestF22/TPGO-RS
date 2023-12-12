function medoids_centerCost_test
global x mu
resetRands()
x=sort(rand(1,8));
mu=linspace(0,1,200);

subplot(3,1,1)
run_test()

subplot(3,1,2)
weights=ones(size(x));
weights(end-2:end)=3;
run_test('weights',weights)

subplot(3,1,3)
weights=ones(size(x));
weights(end-2:end)=3;
run_test('weights',weights,...
    'bias',-2*ones(1,3))

function run_test(varargin)
global x mu
c=medoids_centerCost(x,'xEval',mu,varargin{:});
plot(mu,c)
hold on
plot(x,min(c)*ones(size(x)),'o')
