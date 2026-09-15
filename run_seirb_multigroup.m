function run_seir_multigroup()
% SEIRb Multi-Group Model 
% This script simulates COVID-19 dynamics across 9 age-stratified groups
% incorporating risk-heterogeneity and behavioral feedback.

clear; clc; close all;

% Output directory for high-resolution figures
outdir = 'outputs_seirb_multigroup';
if ~exist(outdir,'dir'), mkdir(outdir); end

%% 1. Baseline Biological and Behavioral Parameters
R0_target = 3;          % Target basic reproduction number for calibration
alpha     = 400000;     % Behavioral feedback strength (Sensitivity to mortality)
tauE      = 3;          % Latent period (days)
tauI      = 10;         % Infectious period (days)
tauF      = 30;         % Time constant for perceived mortality risk (memory)
delta_b   = 0.005;      % Baseline Infection Fatality Ratio (IFR)
fs_all    = 20;         % Standardized font size
lw        = 1.8;        % Standardized line width


%% 2. Population and Contact Matrix Aggregation (C)

% Age cohorts aligned with CDC/Epidemiological reporting
age_labels = {'0-4', '5-17', '18-29', '30-39', '40-49', '50-64', '65-74', '75-84', '85+'};

% Raw 16-group US population data (Census-based)
pop16 = [18149, 19833, 21084, 21105, 20499, 21727, 21873, 22442, ...
         20556, 21327, 21495, 21301, 20121, 17565, 13745, 21200];

% Population weights to redistribute to our age ranges
w4_g2 = 0.6; w4_g3 = 0.4; 
w16_g8 = 0.7; w16_g9 = 0.3;

% Matrix to redistribute age groups
W_col = zeros(16, 9);
W_col(1, 1) = 1;                                        % G1: 0-4
W_col(2:3, 2) = 1; W_col(4, 2) = w4_g2;                 % G2: 5-17
W_col(4, 3) = w4_g3; W_col(5:6, 3) = 1;                 % G3: 18-29
W_col(7:8, 4) = 1;                                      % G4: 30-39
W_col(9:10, 5) = 1;                                     % G5: 40-49
W_col(11:13, 6) = 1;                                    % G6: 50-64
W_col(14:15, 7) = 1;                                    % G7: 65-74
W_col(16, 8) = w16_g8;                                  % G8: 75-84
W_col(16, 9) = w16_g9;                                  % G9: 85+

% Make population vector length 9
pop9 = NaN(1,9);
for j = 1:9
    pop9(j) = pop16*W_col(:,j);
end
m = pop9 / sum(pop9); % Normalized population shares (m_i)
N = numel(age_labels);


% Aggregates Prem et al. 16x16 matrix into 9x9 framework
C16 = readmatrix('Excel_1_(USA).xlsx');
% Matrix of number of contacts (not contact rate)
con_orig = C16.*repmat(pop16',1,16);

% Switch to 16x9 matrix of number of contacts
con_16_9 = NaN(16,9);
for i = 1:16 % loop through old age classes    
    for j = 1:9 % loop through new age classes
        con_16_9(i,j) = con_orig(i,:)*W_col(:,j); % combine contacts of original age group with others
    end
end

% Switch to 9x9 matrix of number of contacts
con_9_9 = NaN(9,9);
for i = 1:9 % loop through new age classes
    for j = 1:9 % loop through new age classes
        con_9_9(i,j) = con_16_9(:,j)'*W_col(:,i); % combine contacts of others 
    end
end

% Make the contact rate matrix
con_rate = NaN(9,9);
for i = 1:9
    for j = 1:9
        con_rate(i,j) = con_9_9(i,j)/pop9(i);
    end
end

C = con_rate;
fprintf('\n--- Contact Matrix ---\n');
disp(C);
max(eig(C));
fprintf('\n--- Max Eigenvalue---\n');
disp(max(eig(C)));

%% 3. IFR Normalization 

% Relative biological risks (Susceptibility rho_i and Mortality theta_i)
rel_inf_risk   = [0.5, 0.7, 1.0, 1.0, 0.9, 0.8, 0.6, 0.7, 0.8];
rel_death_risk = [0.3, 0.1, 1.0, 3.5, 10.0, 25.0, 60.0, 140.0, 360.0];

% Calculates group-specific IFR (delta_i) while maintaining aggregate delta_b
delta_vec = delta_b * (rel_death_risk ./ rel_inf_risk); 
delta_vec = delta_vec * (delta_b / sum(m .* delta_vec)); 


%% 4. R0 Calibration (Pure Spectral Radius Method)
% Anchors the multi-group system to a baseline R0=3 motivated by the NGM approach
NGM_base = zeros(N,N);
for i = 1:N
    for j = 1:N
        NGM_base(i,j) = rel_inf_risk(i) * C(i,j) * (m(i)/m(j)) * tauI;
    end
end
max_eig = max(eig(NGM_base));
beta0_calibrated = R0_target / max_eig; 
beta_vec = beta0_calibrated * rel_inf_risk;

%% 5. Initial Conditions and Simulation Protocol
I0_total = 1e-6; % Initial infectious seed share
I0_vec   = m(:) * I0_total;
E0_vec   = zeros(N, 1); % Consistent with two-group starting protocol
S0_vec   = m(:) - I0_vec;
y0       = [S0_vec; E0_vec; I0_vec; zeros(N,1); zeros(N,1); zeros(N,1)];

tspan    = 0:0.1:365;
params   = struct('beta_vec', beta_vec, 'alpha', alpha, 'tauE', tauE, ...
                  'tauI', tauI, 'tauF', tauF, 'delta', delta_vec, 'm', m, 'C', C);

opts     = odeset('RelTol',1e-10,'AbsTol',1e-12);
[t, y]   = ode45(@(t,y) rhs_paper(t, y, params, N), tspan, y0, opts);

% Parsing compartments for visualization
I = y(:, 2*N+1:3*N); 
R = y(:, 3*N+1:4*N);
D = y(:, 4*N+1:5*N);

%% 6. Visualization 

%% ========= Multi-group infection and death trajectories ===========
fprintf('\n=== Fig 8: Multi-group infection and death trajectories ===\n')

% Figure 8a: Infectious Prevalence Trajectories
f8a = figure;
plot(t, I, 'LineWidth', lw);
hold on; 
plot(t, sum(I,2), 'k--', 'LineWidth', 2.5);
hold off
xlabel('Days');
ylabel('Infectious prevalence'); 
xlim([0 365]);
set(gca, 'FontSize', fs_all);
legend([age_labels, {'Total Population'}], 'Location', 'northeast', 'FontSize', 16); 
save_all_formats(f8a, outdir, 'fig8_results_prevalence');
close

% Figure 8b: Daily Death Distribution
dailyDeaths = (I ./ tauI) .* delta_vec;
f8b = figure;
plot(t, dailyDeaths, 'LineWidth', lw); 
hold on;
plot(t, sum(dailyDeaths,2), 'k--', 'LineWidth', 2.5);
hold off
xlabel('Days');
ylabel('Daily death rate'); 
xlim([0 365]);
set(gca, 'FontSize', fs_all);
save_all_formats(f8b, outdir, 'fig8_results_daily_deaths');
close

%% ========= Multi-group summary outcomes ===========================
fprintf('\n=== Fig 9: Summary outcomes (cumulative deaths and final attack rate) ===\n')

% Figure 9: Summary Outcomes (Cumulative Deaths and Attack Rates)
f9a = figure;
%set(f9a,'Position', [10 10 1400 500]);
% Deaths per cohort
bar(1:N, D(end,:), 'FaceColor', [0.6 0.2 0.2]);
ylabel({'Cumulative deaths at day 365'});
xlabel('Age groups');
set(gca, 'FontSize', fs_all, 'xtick', 1:N, 'xticklabel', age_labels);
xtickangle(90); 
grid off;

f9b = figure;
% Final Attack Rate (R_i + D_i)
bar(1:N, (R(end,:)+D(end,:))./m, 'FaceColor', [0.2 0.6 0.2]);
ylabel('Final attack rate'); 
xlabel('Age groups');
set(gca, 'FontSize', fs_all, 'xtick', 1:N, 'xticklabel', age_labels)
xtickangle(90)
grid off

save_all_formats(f9a, outdir, 'fig9_multigroup_summary_deaths');
save_all_formats(f9b, outdir, 'fig9_multigroup_summary_attack');
close all

fprintf('\nSimulation finished. Publication figures saved in %s\n', outdir);

end

%% RHS Function: Multi-group SEIRb Equations 
function dy = rhs_paper(t, y, p, N)
    S = y(1:N); 
    E = y(N+1:2*N);
    I = y(2*N+1:3*N); 
    D = y(4*N+1:5*N);
    Ft = y(5*N+1:6*N);
    
    % Per-capita infectious pressure accounting for structured contact matrix
    mixing = p.C * (I ./ p.m(:)); 
    
    dS = zeros(N,1);
    dE = zeros(N,1);
    dI = zeros(N,1); 
    dR = zeros(N,1); 
    dD = zeros(N,1); 
    dFt = zeros(N,1);
    
    for i = 1:N
        % Behavioral Feedback Equation: transmission reduced by perceived risk
        beta_eff = p.beta_vec(i) * exp(-p.alpha * (Ft(i) / p.m(i)));
        
        dS(i) = - beta_eff * S(i) * mixing(i);
        dE(i) =   beta_eff * S(i) * mixing(i) - E(i)/p.tauE;
        dI(i) =   E(i)/p.tauE - I(i)/p.tauI;
        dR(i) =   (1 - p.delta(i)) * I(i) / p.tauI;
        dD(i) =   p.delta(i) * I(i) / p.tauI; 
        dFt(i) = (dD(i) - Ft(i)) / p.tauF; % Delay/memory in mortality awareness
    end
    dy = [dS; dE; dI; dR; dD; dFt];
end

%%  Helper: Export Figures 
function save_all_formats(hFig, folder, name)
        savefig(hFig, fullfile(folder, [name, '.fig']));
        saveas(hFig, fullfile(folder, [name, '.png']));
    
end