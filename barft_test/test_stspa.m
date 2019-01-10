% make a target pattern
N = 32;
[x,y] = meshgrid(-N/2:N/2-1);
target = x.^2 + y.^2 <= 5^2;
target = conv2(target,exp(-(x.^2+y.^2)/2),'same');
target = target./max(target(:));

% define a k-space sampling pattern
kmask = true(N);
kmask(1:2:end,:) = false;
Nt = sum(kmask(:));

% load the sensitivities
load fdtdsens
Nc = size(sens,3);
sens = sens(1:2:end,1:2:end,:);

% design the pulses (matlab)
A = exp(-1i*2*pi/N*(x(:)*x(:)' + y(:)*y(:)'));
A = A(:,kmask(:));
Abig = zeros(N*N,Nt*Nc);
for ii = 1:Nc
    senst = sens(:,:,ii);
    Abig(:,(ii-1)*Nt+1:ii*Nt) = bsxfun(@times,senst(:),A);
end
lambda = 1;
rfMat = (Abig'*Abig + lambda*eye(size(Abig,2)))\(Abig'*target(:));

%% design the pulses (barft)
for ii = 1:Nc
   target8(:,:,ii) = target;
   kmask8(:,:,ii) = kmask;
end

pulses = bart('stspa 1', sens, target8, double(kmask8));

for ii = 1:Nc
    holder = squeeze(pulses(:,:,ii));
    pulses2(:,:,ii) = holder(any(holder,2),:);
end

%% check differences
mMat = reshape(Abig*rfMat,[N N]);

mBarft = zeros(size(squeeze(pulses2(:,:,1))));
for ii = 1:Nc
    mBarft = mBarft + fftshift(ifft2(pulses2(:,:,ii)));
end
