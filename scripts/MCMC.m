function samples = MCMC(Avox, bvals, qhat, burn_in, sampling_interval, T)

% Set best-fit model parameters as starting point
x0 = findBestFit(Avox, bvals, qhat);

% Set the sigma to be used for calculating likelihood
Rsq = BallStickSSD_transformation(x0, Avox, bvals, qhat);
sigma = (1/(length(Avox) - length(x0))) * Rsq;

T = burn_in + T*sampling_interval;

samples = zeros(T, length(x0));
samples(1, :) = x0;

for t = 2:T
    
    x_c(1) = normrnd(samples(t-1, 1), 1e-01);
    x_c(2) = normrnd(samples(t-1, 2), 1e-05);
    while any(x_c(1:2) < 0)
        x_c(1) = normrnd(samples(t-1, 1), 1e-01);
        x_c(2) = normrnd(samples(t-1, 2), 1e-05);
    end
    x_c(3) = normrnd(samples(t-1, 3), 1e-04);
    while x_c(3) < 0 || x_c(3) > 1
        x_c(3) = normrnd(samples(t-1, 3), 1e-04);
    end
    x_c(4) = normrnd(samples(t-1, 4), 1e-05);
    x_c(5) = normrnd(samples(t-1, 5), 1e-05);

    % Compute the likelihood ratio alpha
    alpha = computeAlpha(x_c, samples(t-1,:), Avox, bvals, qhat, sigma);
    
    if alpha > rand
       samples(t, :) = x_c;
    else
       samples(t, :) = samples(t-1, :);
    end
    
end

% Discard examples for the burn_in and sampling_interval
samples(1:burn_in, :) = [];

sample_idxs = zeros(1, (T - burn_in)/sampling_interval);
i = 1;
for t=1:sampling_interval:(T - burn_in)
    sample_idxs(i) = t;
    i= i+1;
end

samples = samples(sample_idxs, :);

end

function alpha = computeAlpha(x_c, x_t, Avox, bvals, qhat, sigma)

% Compute model signals
S_t = BallStickModel(x_t, bvals, qhat);
S_c = BallStickModel(x_c, bvals, qhat);

% Compute SSDs
SSD_t = sum((Avox - S_t').^2);
SSD_c = sum((Avox - S_c').^2);

% Obtain theta
theta_c = x_c(4);
theta_t = x_t(4);

alpha = exp(log(sin(theta_c)) - log(sin(theta_t)) + (1/(2*sigma^2))*(SSD_t - SSD_c));

end

function parameter_hat = findBestFit(Avox, bvals, qhat)

% Define a starting point for the non-linear fit
startx = [sqrt(3.5e+00) sqrt(3e-03) sqrt(-log(2.5e-01)) 0 0];

% Define various options for the non-linear fitting
% algorithm.
h=optimset('MaxFunEvals',20000,...
 'Algorithm','quasi-newton',...
 'TolX',1e-10,...
 'TolFun',1e-10);

% Now run the fitting
[parameter_hat,~,~,~]=fminunc('BallStickSSD_transformation',startx,h,Avox,bvals,qhat);

% Transform parameters back to original values
parameter_hat(:, 1) = parameter_hat(:, 1)^2;
parameter_hat(:, 2) = parameter_hat(:, 2)^2;
parameter_hat(:, 3) = exp(-parameter_hat(:, 3)^2);

end

function sumRes = BallStickSSD_transformation(x, Avox, bvals, qhat)

% Transform parameters back to original values
S0 = x(1).^2;
diff = x(2).^2;
f = exp(-(x(3).^2));
theta = x(4);
phi = x(5);

% Synthesize the signals according to the model
fibdir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat) 1])');
S = S0*(f*exp(-bvals*diff.*(fibdotgrad.^2)) + (1-f)*exp(-bvals*diff));

% Compute the sum of square differences
sumRes = sum((Avox - S').^2);

end 

function S = BallStickModel(x, bvals, qhat)

% Transform parameters back to original values
S0 = x(1);
diff = x(2);
f = x(3);
theta = x(4);
phi = x(5);

% Synthesize the signals according to the model
fibdir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat) 1])');
S = S0*(f*exp(-bvals*diff.*(fibdotgrad.^2)) + (1-f)*exp(-bvals*diff));

end