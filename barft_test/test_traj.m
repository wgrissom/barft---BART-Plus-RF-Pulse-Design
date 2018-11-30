figure

% design a single leaf spiral k-space
traj = bart('traj -s -x 32 -a 3 -I 1');
subplot(121)
plot(traj(1,:)+1i*traj(2,:));
title 'One interleaf'
axis equal

% design a four-leaf spiral k-space
traj = bart('traj -s -x 32 -a 3 -I 4');
subplot(122)
hold on
plot(traj(1,:,1)+1i*traj(2,:,1));
plot(traj(1,:,2)+1i*traj(2,:,2));
plot(traj(1,:,3)+1i*traj(2,:,3));
plot(traj(1,:,4)+1i*traj(2,:,4));
title 'Four interleaves'
axis equal
