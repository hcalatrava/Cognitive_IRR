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
% if config.clairvoyant
%     sgtitle('ROC (clairvoyant detector)')
% else
%     sgtitle('ROC (X MLE)')
% end

%% Adaptive Waveform Design (AWD) Module

% Parameter values for the test
config.L = 5;
config.SNR_predefined = -10;
config.enable_awd = 1;
config.N_mc = 50000;

% Get ROC curve from configuration, with and without enabling AWD

[pfa, pd, pfa_opt, pd_opt] = get_ROC_from_config(config, LOS_vector, T_rangecell)

% Plot ROC for this configuration
semilogx(pfa(:,1), pd(:,1), 'displayname', ['AWD disabled'], ...
    'linewidth', 1.5, 'color', config.plot_color(2), ...
    'linestyle', config.plot_linestyle(2));
hold on,
semilogx(pfa_opt(:,1), pd_opt(:,1), 'displayname', ['AWD enabled'], ...
    'linewidth', 1.5, 'color', config.plot_color(1), ...
    'linestyle', config.plot_linestyle(1));

legend('show', 'location', 'southeast')
grid on
axis([0 1 0 1])
xlim([0.005, 1]) % same x axis limit as in paper #2
title(['SNR = ', num2str(config.SNR_predefined), ' dB, L = ', num2str(config.L), ' subcarriers'])
xlabel('FPR P(D=1 | H0)')
ylabel('TPR P(D=1 | H1)')