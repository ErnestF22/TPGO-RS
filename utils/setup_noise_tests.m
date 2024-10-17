function input_data_nois = setup_noise_tests(input_data, sigma, mu)
%Returns noisy data according to sigma, mu arguments

% To compute which method performs best: graph with x-axis as input noise
% (e.g. $\sigma$-varying Gaussian noise), y-axis percentage error on result;
% always averaging multiple cases (around 30 cases per each $\sigma$) for each edge


if sigma==0 && mu == 0
    input_data_nois = input_data;
    return;
end

%input_data_nois = awgn(input_data); % to add white Gaussian noise
d = size(input_data, 1);
num_tijs = size(input_data, 2);

tmp = sigma * randn(1, num_tijs); %assume here: Gaussian, mean mu, std sigma, frequency continuous
input_data_nois = mu + (1-tmp/2).*input_data;


% setup plots


end

