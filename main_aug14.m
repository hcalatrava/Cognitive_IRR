%% Initialize scenario (paper 0)
close all,
clear all,
clc,
config = load_config(0); % Load configuration

format long e; % for displaying purposes

% Test for different number of subcarriers L
% Results are plotted on the bottom subfigure

% Configuration
rng('default') % Initialize seed

LOS_vector = get_LOS_vector(config); % get LOS vector
T_rangecell = sqrt(config.target_dist_east^2 + config.target_dist_north^2)*2/config.c; % roundtrip time

% Create figure
figure()
if config.clairvoyant
    sgtitle('ROC (clairvoyant detector)')
else
    sgtitle('ROC (X MLE)')
end

config.enable_awd = 0;

%% Different number of subcarriers    

L_vector = [1, 3, 5]; % values of L pqfor which a ROC curve is computed

subplot(2,1,2)
% config.N_mc = 50000;

% Loop through different number of subcarriers L
for idx_L = 1:length(L_vector)

    config.L = L_vector(idx_L); % update number of subcarriers

    [pfa, pd] = get_ROC_from_config(config, LOS_vector, T_rangecell)

    % Plot ROC for this configuration
    semilogx(pfa(:,1), pd(:,1), 'displayname', ['L = ', num2str(L_vector(idx_L))], ...
        'linewidth', 1.5, 'color', config.plot_color(idx_L), ...
        'linestyle', config.plot_linestyle(idx_L));
    hold on,
end

legend('show', 'location', 'southeast')
grid on
% title('ROC')
axis([0 1 0 1])
xlim([0.005, 1]) % same x axis limit as in paper #2
title(['SNR = ', num2str(config.SNR_predefined), ' dB'])
xlabel('FPR P(D=1 | H0)')
ylabel('TPR P(D=1 | H1)')

%% Different predefined SNR values

config.L = 2;
SNR_predefined_vector = [-15, -10, -5]; 

subplot(2,1,1)

% Loop through different number of predefined SNR values
for idx_SNR = 1:length(SNR_predefined_vector)

    % Generate OFDM measurements
    config.SNR_predefined = SNR_predefined_vector(idx_SNR); % update number of subcarriers

    [pfa, pd] = get_ROC_from_config(config, LOS_vector, T_rangecell);

    % Plot ROC for this configuration
    semilogx(pfa(:,1), pd(:,1), 'displayname', ['SNR = ', num2str(SNR_predefined_vector(idx_SNR))], ...
        'linewidth', 1.5, 'color', config.plot_color(idx_SNR), ...
        'linestyle', config.plot_linestyle(idx_SNR));
    hold on,
end

legend('show', 'location', 'southeast')
grid on
% title('ROC')
axis([0 1 0 1])
xlim([0.005, 1]) % same x axis limit as in paper #2
xlabel('FPR P(D=1 | H0)')
title(['L = ', num2str(config.L), ' subcarriers'])
ylabel('TPR P(D=1 | H1)')