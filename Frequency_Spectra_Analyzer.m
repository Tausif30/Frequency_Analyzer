% Laser Transfer Function Model
%  Input : Ch1 (Measured Photo detector pulse)
%  Output: Ch3 and Ch4 (measured output signals)
%  Method: H(f) = Y(f)/X(f) with Wiener regularization

% Load Data
data = readmatrix('Experiment_A.csv', 'NumHeaderLines', 2);
t    = data(:, 1);      % Time [s]
u_t  = data(:, 2);      % Ch1 = Input  [V]
y3_t = data(:, 3);      % Ch3 = Output [V]
y4_t = data(:, 4);      % Ch4 = Output [V]

% Center Vector
[~, idx_peak] = max(u_t);
t = t - t(idx_peak);

% Fast Fourier Transform Setup
dt = t(2) - t(1);           % Sampling interval [s]
Fs = 1 / dt;                % Sampling frequency [Hz]
N  = length(t);             % Number of samples

% Frequency axis (0 to Nyquist)
f_onesided = (0 : N-1) * (Fs / N);                  % for computation
f_shift    = (-N/2 : N/2-1) * (Fs / N);             % for plotting (centred)
fprintf('fs = %.2f GHz | Nyquist = %.2f GHz | df = %.4f GHz\n', ...
    Fs/1e9, Fs/2e9, (Fs/N)/1e9);

% Computing FFT
X_f  = fft(u_t);     % Input  spectrum (Ch1 measured)
Y3_f = fft(y3_t);    % Output spectrum (Ch3)
Y4_f = fft(y4_t);    % Output spectrum (Ch4)

% Transfer Function  H(s) = Y(s) / X()
%  Used Wiener-style regularization:
%  H(s) = X*(s) · Y(s) / ( |X(s)|² + ε )
%  This avoids divide-by-zero and smoothly tapers noise where X(f) is small.
epsilon = (max(abs(X_f)) * 0.01)^2;    % Wiener Regularization (1% power)
 
H3_f = (conj(X_f) .* Y3_f) ./ (abs(X_f).^2 + epsilon);
H4_f = (conj(X_f) .* Y4_f) ./ (abs(X_f).^2 + epsilon);

% Delay Correction (Automated Calculation)
[~, idx_input]    = max(u_t);          % index of Ch1 peak
[~, idx_ch3_peak] = max(abs(y3_t));    % index of Ch3 peak
[~, idx_ch4_peak] = max(abs(y4_t));    % index of Ch4 peak

delay_ch3 = t(idx_ch3_peak) - t(idx_input); 
delay_ch4 = t(idx_ch4_peak) - t(idx_input); 
fprintf('Calculated Delay Ch3 = %.3f ns | Ch4 = %.3f ns\n', delay_ch3*1e9, delay_ch4*1e9);

H3_f_corr = H3_f .* exp(+1j * 2*pi * f_onesided' * delay_ch3);
H4_f_corr = H4_f .* exp(+1j * 2*pi * f_onesided' * delay_ch4);

% Shift spectra for display
X_f_shift   = fftshift(X_f);
Y3_f_shift  = fftshift(Y3_f);
Y4_f_shift  = fftshift(Y4_f);
H3_f_shift  = fftshift(H3_f_corr);
H4_f_shift  = fftshift(H4_f_corr);

% H(s) Phase
H3_phase = unwrap(angle(H3_f_shift));
H4_phase = unwrap(angle(H4_f_shift));
df = f_shift(2) - f_shift(1);
gd3 = -diff(H3_phase) / (2*pi*df);     % Group delay Ch3 [s]
gd4 = -diff(H4_phase) / (2*pi*df);     % Group delay Ch4 [s]
f_gd = f_shift(1:end-1) + df/2;        % Midpoint frequency axis

% Plot
GHz = 1e9;
 
figure('Name','Laser Transfer Function','NumberTitle','off', ...
       'Position',[60 40 1200 860]);

% Panel 1: Time Domain (Input & Outputs)
subplot(3,2,1);
yyaxis left;
plot(t*1e9, u_t*1e3, 'b', 'LineWidth', 1.6, 'DisplayName','Ch1 Input [mV]');
ylabel('Ch1 Amplitude [mV]');

yyaxis right;
plot(t*1e9, y3_t*1e3, 'r', 'LineWidth', 1.3, 'DisplayName','Ch3 Output [mV]'); hold on;
plot(t*1e9, y4_t*1e3, 'g', 'LineWidth', 1.3, 'DisplayName','Ch4 Output [mV]');
ylabel('Ch3/Ch4 Amplitude [mV]');

xlabel('Time [ns]'); title('Time Domain');
legend('Location','best'); grid on;

% Panel 2: Input Pulse
subplot(3,2,2);
plot(t*1e9, u_t*1e3, 'b', 'LineWidth', 1.6, 'DisplayName','Ch1 Measured');
xlabel('Time [ns]'); ylabel('Amplitude [mV]');
title('Input Pulse');
legend('Location','best'); grid on;
xlim([-10 10]);

% Normalized Frequency Spectra
subplot(3,2,3);
plot(f_shift/GHz, abs(X_f_shift)/max(abs(X_f_shift)),  'b',  'LineWidth',1.5, 'DisplayName','Ch1 Input X(f)'); hold on;
plot(f_shift/GHz, abs(Y3_f_shift)/max(abs(Y3_f_shift)),'r',  'LineWidth',1.3, 'DisplayName','Ch3 Output Y(f)');
plot(f_shift/GHz, abs(Y4_f_shift)/max(abs(Y4_f_shift)),'g',  'LineWidth',1.3, 'DisplayName','Ch4 Output Y(f)');
xlabel('Frequency [GHz]'); ylabel('Normalized Magnitude');
title('Normalized Frequency Spectra');
legend('Location','best'); grid on; xlim([-5 5]);

% H(s) in dB
subplot(3,2,4);
H3_dB = 20*log10(abs(H3_f_shift) + eps);
H4_dB = 20*log10(abs(H4_f_shift) + eps);
plot(f_shift/GHz, H3_dB, 'r', 'LineWidth',1.5, 'DisplayName','|H_{Ch3}(f)| dB'); hold on;
plot(f_shift/GHz, H4_dB, 'g', 'LineWidth',1.5, 'DisplayName','|H_{Ch4}(f)| dB');
xlabel('Frequency [GHz]'); ylabel('|H(f)| [dB]');
title('Transfer Function — Magnitude (dB)');
legend('Location','best'); grid on; xlim([-5 5]);
top_dB = max([max(H3_dB) max(H4_dB)]);
ylim([top_dB-40, top_dB+5]);

% H(s) Phase
subplot(3,2,5);
plot(f_shift/GHz, H3_phase, 'r', 'LineWidth',1.5, 'DisplayName','Phase H_{Ch3}'); hold on;
plot(f_shift/GHz, H4_phase, 'g', 'LineWidth',1.5, 'DisplayName','Phase H_{Ch4}');
xlabel('Frequency [GHz]'); ylabel('Phase [rad]');
title('Transfer Function — Phase (Delay Corrected)');
legend('Location','best'); grid on; xlim([-5 5]);

% Group Delay
subplot(3,2,6);
plot(f_gd/GHz, gd3*1e9, 'r', 'LineWidth',1.5, 'DisplayName','GD Ch3'); hold on;
plot(f_gd/GHz, gd4*1e9, 'g', 'LineWidth',1.5, 'DisplayName','GD Ch4');
yline(0,'--k','LineWidth',0.8);
xlabel('Frequency [GHz]'); ylabel('Group Delay [ns]');
title('Group Delay \tau_g(f) = -d\phi/d\omega');
legend('Location','best'); grid on; xlim([-5 5]);
 
sgtitle('Laser Transfer Function — Experiment A  |  H(f) = Y(f)/X(f)', ...
        'FontSize',13,'FontWeight','bold');

% Save Results
pos_idx = f_onesided >= 0 & f_onesided <= 5e9;
tmp = f_onesided(pos_idx);   f_export = reshape(tmp, [], 1);
tmp = abs(H3_f(pos_idx));    H3_mag   = reshape(tmp, [], 1);
tmp = angle(H3_f(pos_idx));  H3_phase = reshape(tmp, [], 1);
tmp = abs(H4_f(pos_idx));    H4_mag   = reshape(tmp, [], 1);
tmp = angle(H4_f(pos_idx));  H4_phase = reshape(tmp, [], 1);
 
H_export = table(f_export, H3_mag, H3_phase, H4_mag, H4_phase, ...
    'VariableNames', {'Frequency_Hz', ...
                      'H_Ch3_magnitude','H_Ch3_phase_rad', ...
                      'H_Ch4_magnitude','H_Ch4_phase_rad'});
 
writetable(H_export, 'H_f_results.csv');
fprintf('Transfer function saved to H_f_results.csv\n');