function sumRes = ZeppelinStickSSD(x, Avox, bvals, qhat)

% Extract the parameters
S0 = x(1)^2;
diff = x(1)^2*x(2)^2 ;
f = exp(-x(3)^2);
lambda1 = diff;
lambda2 = x(4)^2;
theta = x(5);
phi = x(6);

% Synthesize the signals according to the model
fibdir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat) 1])');
S = S0*(f*exp(-bvals*diff.*(fibdotgrad.^2)) + (1-f)*exp(-bvals.*(lambda2 + (lambda1 - lambda2).*fibdotgrad.^2)));

% Compute the sum of square differences
sumRes = sum((Avox - S').^2);

end