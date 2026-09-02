mask= abs(az)<2 & gs>1;
max_pot_load=max(PotRrLoad,max(PotRlLoad,max(PotFrLoad,PotFlLoad)));
max_accel_load=max(AccelRrLoad,max(AccelRlLoad,max(AccelFrLoad,AccelFlLoad)));
figure()
subplot(1,2,1);
plot(max_pot_load(mask),mu(mask),'.');
title('mu')
subplot(1,2,2);
plot(max_pot_load(mask),gs(mask),'.');
title('gs')

% figure()
% subplot(1,2,1);
% plot(max_accel_load(mask),mu(mask),'.');
% title('mu')
% subplot(1,2,2);
% plot(max_accel_load(mask),gs(mask),'.');
% title('gs')