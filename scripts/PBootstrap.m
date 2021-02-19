function samples = PBootstrap(Avox, bvals, qhat, n_runs_per_voxel, T)

% Find best-fit model parameters
[fitted_params, RESNORM] = findBestFit(n_runs_per_voxel, Avox, bvals, qhat);

% Computing model signals
S = BallStickModel(fitted_params, bvals, qhat);

% Compute sigma
sigma = sqrt((1/(length(Avox) - length(fitted_params))) * sum((Avox - S').^2));

% store samples
samples = zeros(length(fitted_params), T);

for t=1:T
    
    % Sample random noise values from N(0, sigma)
    E = normrnd(0, sigma, 1, 108);
    
    % Synthesize bootstrap dataset
    Ahat = S + E;
    Ahat = reshape(Ahat, 108, 1);
    
    % Estimate the best fit to A_hat
    [x, RESNORM] = findBestFit(n_runs_per_voxel, Ahat, bvals, qhat);
    
    % Transform parameters back to original values
    x(:, 1) = x(:, 1)^2;
    x(:, 2) = x(:, 2)^2;
    x(:, 3) = exp(-x(:, 3)^2);
    
    % Store sample values
    samples(:, t) = x;
    
end

end

function [best_params, min_error] = findBestFit(n_runs, Avox, bvals, qhat)

min_error = 1e+10;
run = 1;

while run <= n_runs
    
    fprintf('Run %i\n', run);
    
    % Define a set start point for the optimization
    startx = [3.5e+00 3e-03 2.5e-01 0 0];
    startx(:, 1) = startx(:, 1) + 1e+03*randn;
    startx(:, 2) = startx(:, 2) + 1e-03*randn;
    startx(:, 3) = startx(:, 3) + 1e-01*randn;
    startx(:, 4) = startx(:, 4) + 1e-01*randn;
    startx(:, 5) = startx(:, 5) + 1e-01*randn;
    
    % Perform inverse transform on the parameters
    startx(:, 1:2) = sqrt(startx(:, 1:2));
    startx(:, 3) = sqrt(-log(startx(:, 3)));

    % Define various options for the non-linear fitting
    % algorithm.
    h=optimset('MaxFunEvals',20000,...
     'Algorithm','quasi-newton',...
     'TolX',1e-10,...
     'TolFun',1e-10, 'Display', 'off');
    
    % Now run the fitting
    try
        [parameter_hat,RESNORM,~,~]=fminunc('BallStickSSD_transformation',startx,h,Avox,bvals,qhat);
    catch
        continue;
    end
    
    % Save the parameters if they provide the smallest RESNORM
    if RESNORM < min_error
        min_error = RESNORM;
        best_params = parameter_hat;
    end
    
    run = run + 1;
    
end

end

function sumRes = BallStickSSD_transformation(x, Avox, bvals, qhat)

S = BallStickModel(x, bvals, qhat);

% Compute the sum of square differences
sumRes = sum((Avox - S').^2);

end 

function S = BallStickModel(x, bvals, qhat)

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

end