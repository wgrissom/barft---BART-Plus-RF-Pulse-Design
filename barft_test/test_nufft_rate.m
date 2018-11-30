% don't forget to run bart startup!
% also need a test for multiple shots/interleaves

% make a phantom
img = bart('phantom');

% make a non-Cartesian trajectory
traj = bart('traj -s -x 32 -a 3 -I 1');

% make a rate map (radians/s) and define T
[x,y] = ndgrid(-64:63);
rate = 2*pi*100*exp(-((x+32).^2 + y.^2)/100);
T = 10; % ms

% calculate dft in matlab
kspMat = exp(-1i*2*pi/128*(traj(1,:)'*x(:)' + traj(2,:)'*y(:)')) * img(:);

% calculate dft using bart
writecfl('img',img); % so we can test on cmd line
writecfl('spiralTraj',traj);
writecfl('rate',rate);
kspBart = bart('nufft -s',traj,img);  % -s tells nufft to do an nudft
printf('Bart DFT nrmse with no rate map: %0.1d',nrmse(kspMat(:),kspBart(:),1));

% calculate dft with rate in matlab
t = 0:T/(size(traj,2)-1):T;
kspMatRate = exp(-1i*2*pi/128*(traj(1,:)'*x(:)' + traj(2,:)'*y(:)') ...
    + 1i*t(:)*rate(:)') * img(:);

% calculate dft with rate using bart
