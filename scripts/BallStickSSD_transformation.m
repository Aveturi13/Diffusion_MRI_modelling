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