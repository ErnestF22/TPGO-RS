function POCLambdaLogRegularizer
a=2; 
fOld=@(l) relu_som(1-l)^2; 
lGrid=linspace(1/a+0.01,1.5);
%funCheckDer(@selfConcordant,@selfConcordantDer,linspace(0.1,1))

figure(1)
subplot(2,1,1)
funCheckDer(@(l) fRegularizer(l,a),@(l) fRegularizerDer(l,a),lGrid)
subplot(2,1,2)
funCheckDer(@(l) fRegularizerDer(l,a),@(l) fRegularizerDerDer(l,a),lGrid)
disp('Value of regularizer and first two derivatives should be zero at lambda=1')
disp([fRegularizer(1,a) fRegularizerDer(1,a) fRegularizerDerDer(1,a)])

figure(2)
funPlot(fOld,lGrid)
hold on
funPlot(@(l) fRegularizer(l,a),lGrid)
hold off
legend('ReLU^2','log-barrier')

function f=selfConcordant(x)
% see https://en.wikipedia.org/wiki/Self-concordant_function for options
p=0.5;
g=-x^p;
f=-log(g)-log(x);

function df=selfConcordantDer(x)
p=0.5; % need to match parameter in selfConcordant
g=-x^p;
dg=-p*x^(p-1);
df=-dg/g-1/x;


function f=fRegularizer(l,a)
b=-a/(a-1)^2;
if l<1
    f=-1/a*log(a*l-1)...
        +1/(a-1)*(l-1)...
        +b/2*(l-1)^2;
else
    f=0;
end

function f=fRegularizerDer(l,a)
b=-a/(a-1)^2;
if l<=1
    f=-1/(a*l-1)...
        +1/(a-1)...
        +b*(l-1);
else
    f=0;
end

function f=fRegularizerDerDer(l,a)
b=-a/(a-1)^2;
if l<=1
    f=a/(a*l-1)^2 ...
      +0 ...
      +b;
else
    f=0;
end
