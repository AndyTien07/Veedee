dt=0.01;
Fs = 1/dt;          
nfft = 2^11;    
window = hann(2^8);   
noverlap = 2^8 * 3/4; 
Start_Time=50*100;
End_Time=150*100;

% signal1=FL_pos+FR_pos;
% signal2=RL_pos+RR_pos;
tit="Pitch"

[Pxx, fp] = cpsd(signal1,signal2, window, noverlap, nfft, Fs);
[coherence, fc] = mscohere(signal1, signal2, window, noverlap, nfft, Fs);

phase_wrapped = rad2deg(angle(Pxx));
phase_unwrapped = rad2deg(unwrap(angle(Pxx)));

coh_min = 0.3;
phase_masked = phase_wrapped;
phase_masked(coherence < coh_min) = NaN;
% Plotting the wrapped and unwrapped phase
figure()
plot(fp, phase_wrapped, 'r', 'DisplayName', 'Wrapped Phase', 'LineWidth', 2);
hold on;
plot(fp, phase_unwrapped, 'b', 'DisplayName', 'Unwrapped Phase', 'LineWidth', 2);
xlim([0, 15]);
title([tit, ' Phase'], 'FontSize', 18);
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');
legend show;
grid on;

figure()
% Plotting in linear g^2/Hz
plot(fp, Pxx, 'DisplayName', 'csd','LineWidth',2);
xlim([0 15]);
title([tit,'csd'],FontSize=18)
xlabel('Frequency (Hz)');
xlim([0,15])
ylabel('PSD (Units are hard)'); % Updated Label
grid on;

figure()

plot(fc, coherence)
title([tit,'coherence'],FontSize=18)
xlabel('Frequency (Hz)')
xlim([0,15])
ylabel('Coherence (0 to 1)')
grid on

