function parameter_maps = createParameterMaps(image, bvals, qhat, runs_per_vox)

% Stores the parameter valuesye
parameter_maps = zeros(145, 174, 6);

for i=1:145
    for j=1:174
        % Select a single voxel
        Avox = image(:, i, j);
        
        if min(Avox) > 0
            % Perform multiple runs over voxel
            [parameter_hat, RESNORM] = multipleRuns(runs_per_vox, Avox, bvals, qhat);
            parameter_maps(i, j, 1:5) = parameter_hat;
            parameter_maps(i, j, 6) = RESNORM;
        end
    end
end

end

function [best_params, minError] = multipleRuns(n_runs, Avox, bvals, qhat)

% This will be updated based on whether a new minimum is found
minError = inf;
n = 1;

while n <= n_runs

    % Define a starting point for the non-linear fit
    startx = [3.5e+00 3e-03 2.5e-01 0 0];
    
    % Add gaussian noise to each of parameters
    startx_noisy(:, 1) = startx(:, 1) + 1e+03*randn;
    startx_noisy(:, 2) = startx(:, 2) + 1e-03*randn;
    while any(startx_noisy(:, 1:2) < 0)
        startx_noisy(:, 1) = startx(:, 1) + 1e+03*randn;
        startx_noisy(:, 2) = startx(:, 2) + 1e-03*randn;
    end
    startx_noisy(:, 3) = startx(:, 3) + 1e-01*randn;
    while startx_noisy(:, 3) > 1 || startx_noisy(:, 3) < 0
        startx_noisy(:, 3) = startx(:, 3) + 1e-01*randn;
    end
    startx_noisy(:, 4) = startx(:, 4) + 1e-01*randn;
    startx_noisy(:, 5) = startx(:, 5) + 1e-01*randn;
    
    % Perform inverse transformation on first three parameters
    startx_noisy(:, 1:2) = sqrt(startx_noisy(:, 1:2));
    startx_noisy(:, 3) = sqrt(-log(startx_noisy(:, 3)));
    
    % Define various options for the non-linear fitting
    % algorithm.
    h=optimset('MaxFunEvals',20000,...
     'Algorithm','quasi-newton',...
     'TolX',1e-10,...
     'TolFun',1e-10, 'Display', 'off');
    
    % Now run the fitting
    try
        [parameter_hat,RESNORM,~,~]=fminunc('BallStickSSD_transformation',startx_noisy,h,Avox,bvals,qhat); 
    catch
        continue;
    end
    
    % Perform transformation back on parameters
    parameter_hat(:, 1) = parameter_hat(:, 1)^2;
    parameter_hat(:, 2) = parameter_hat(:, 2)^2;
    parameter_hat(:, 3) = exp(-parameter_hat(:, 3)^2);

    if RESNORM < minError
        minError = RESNORM;
        best_params = parameter_hat;
    end
    
    n = n + 1;
end

end

function sumRes = BallStickSSD(x, bvals, qhat)

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