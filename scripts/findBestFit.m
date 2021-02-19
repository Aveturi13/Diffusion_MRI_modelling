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
[parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc('BallStickSSD_transformation',startx,h,Avox,bvals,qhat);

% Transform parameters back to original values
parameter_hat(:, 1) = parameter_hat(:, 1)^2;
parameter_hat(:, 2) = parameter_hat(:, 2)^2;
parameter_hat(:, 3) = exp(-1 * parameter_hat(:, 3)^2);

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