function run_seirb_twogroup()
% SEIRb Two-Group Model
% This script simulates dynamics across two groups
% incorporating risk-heterogeneity and behavioral feedback with
% differing assumptions in the two groups but

clearvars; clc; close all;

outdir = 'outputs_seirb_twogroup';
if ~exist(outdir,'dir'), mkdir(outdir); end

% Set the following to 1 for the outputs to run
type.fig3 = 1; % Vary mortality - individual scenarios
type.fig4 = 1; % Vary mortality (theta) across a range
type.fig5 = 1; % Vary susceptibility (rho) across a range
type.fig6 = 1; % Vary mortality (theta) and contact structure together
type.fig7 = 1; % Vary susceptibility (rho) and contact structure together
type.figA1 = 1; % Vary behavior (alpha) - exponential form
type.figB1 = 1; % Vary susceptibility (rho) and contact structure together
type.figC1 = 1;% Vary behavior (alpha) - fractional form
type.figC2=1;% Vary mortality - individual scenarios - fractional form
type.figC3=1; % Vary mortality (theta) across a range (Fractional Analysis)
type.figC4=1;% Vary susceptibility (rho) across a range (Fractional Analysis)
type.figC5=1;% Vary mortality (theta) and contact structure together (Fractional Analysis)
type.figC6=1;% Vary susceptibility (rho) and contact structure together (Fractional Analysis)
type.figD1 = 1; % Vary mortality (theta) across a range with aggregate mortality
% Base parameters
beta0  = 0.3;
alpha  = 200000;
tauE   = 3;
tauI   = 10;
tauF   = 30;
delta0 = 0.005;
k_y    = 1.0;      % two-group k_y (base)
k_o    = 1.0;      % two-group k_o (base)
C_base = [0.5 0.5; 0.5 0.5];  % homogeneous contacts
C_het  = [0.9 0.1; 0.1 0.9];

% Two-group baseline params
m  = 0.5; ry = m; ro = 1 - m;
p2 = struct('beta0',beta0,'alpha',alpha,'tauE',tauE,'tauI',tauI,'tauF',tauF, ...
    'delta_y',delta0,'delta_o',delta0,'k_y',k_y,'k_o',k_o, ...
    'ry',ry,'ro',ro,'C',C_base);
p2.beta_y = beta0;
p2.beta_o = beta0;


p2.use_fractional = 0; % Keep this 0 for baseline exponential. Set to 1 only to use fractional form of 

tgrid   = 0:0.1:365;
tspan   = tgrid;
optsBase= odeset('RelTol',1e-12,'AbsTol',1e-14);

% Initial conditions, two-group
I0 = 1e-6;
Iy0 = ry*I0; Io0 = ro*I0;
Sy0 = ry - Iy0; So0 = ro - Io0;
y0_2 = [Sy0; 0; Iy0; 0; 0;  So0; 0; Io0; 0; 0];   % [Sy Ey Iy Ry Fty  So Eo Io Ro Fto]
opts2 = odeset(optsBase,'NonNegative',1:numel(y0_2));

% Define fontsize and line width
fs_all   = 20;
fs_legend = 16;
lw = 1.8;

% Parameter vectors
Theta_list  = [1, 0.5, 0.25, 0.1, 0.05, 0.02, 0.01, 1e-3, 1e-4]; % mortality ratio
rho_list    = [1, 0.5, 0.25, 0.1, 0.05, 0.02, 0.01, 1e-3, 1e-4]; % susceptibility ratio
alpha_list = [0, 200000, 400000]; % behavioral sensitivity

%% ========= Delta scenarios (group mortality) ==========================
if type.fig3==1
    fprintf('\n=== Fig 3: Delta scenarios (vary group mortality) ===\n')

    m  = 0.5; ry = m; ro = 1 - m;
    Iy0 = ry*I0; Io0 = ro*I0;
    Sy0 = ry - Iy0; So0 = ro - Io0;
    y0_2 = [Sy0; 0; Iy0; 0; 0;  So0; 0; Io0; 0; 0];

    % Homogeneous: dy=do=delta0  (m=0.5)
    pA = p2; pA.ry=ry; pA.ro=ro; pA.delta_y=delta0; pA.delta_o=delta0;
    [tA, yA] = ode45(@(t,y) rhs_twogroup(t,y,pA), tspan, y0_2, opts2);
    % infectious prevalence per group
    IyA = yA(:,3);
    IoA = yA(:,8);
    ItotA = IyA+IoA; % total infectious prevalence
    % daily death proxy per group
    FyA = pA.delta_y .* IyA ./ tauI;
    FoA = pA.delta_o .* IoA ./ tauI;
    fA  = FyA + FoA;   % total daily deaths

    % Heterogeneous: dy=0, do=2*delta0
    pB = p2; pB.ry=ry; pB.ro=ro; pB.delta_y=0; pB.delta_o=2*delta0;
    [tB, yB] = ode45(@(t,y) rhs_twogroup(t,y,pB), tspan, y0_2, opts2);
    % infectious prevalence per group
    IyB = yB(:,3);
    IoB = yB(:,8);
    ItotB = IyB+IoB;
    % daily death proxy per group
    FyB = pB.delta_y .* IyB ./ tauI;
    FoB = pB.delta_o .* IoB ./ tauI;
    fB  = FyB + FoB;   % total daily deaths (Scenario B)
    
    % Plot 1: Vary mortality - prevalence (Upgraded to Tiledlayout)
    % -------------------------------------------------------------
    f_prev = figure('Name','Vary mortality - prevelance');
    set(f_prev, 'Position', [10 10 900 1100]); pause(.1)
    
    t_layout1 = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
    
    nexttile;
    plot(tA, IyA, '-',  'LineWidth', lw); 
    hold on;
    plot(tB, IyB, '--', 'LineWidth', lw);
    xlim([0 365]); 
    ylim([0 0.06]);
    yticks([0 0.02 0.04 0.06]);
    ylabel({'Prevalence','(young group)'});
    set(gca,'FontSize',fs_all);
    legend({'Homogeneous: \delta_y=\delta_o','Heterogeneous: \delta_y=0, \delta_o=2\delta_0'}, ...
        'Location','northeast','FontSize',fs_legend);
    
    nexttile;
    plot(tA, IoA, '-',  'LineWidth', lw);
    hold on;
    plot(tB, IoB, '--', 'LineWidth', lw);
    xlim([0 365]); 
    ylim([0 0.06]); 
    yticks([0 0.02 0.04 0.06]);
    ylabel({'Prevalence','(old group)'});
    set(gca,'FontSize',fs_all);
    
    nexttile;
    plot(tA, ItotA, '-',  'LineWidth', lw); 
    hold on;
    plot(tB, ItotB, '--', 'LineWidth', lw);
    xlabel('Days');
    ylabel({'Prevalence','(total)'});
    xlim([0 365]);
    ylim([0 0.06]); 
    yticks([0 0.02 0.04 0.06]);
    set(gca,'FontSize',fs_all);
    
    savepng(fullfile(outdir,'fig3_delta_scenarios_I.png'));
    savefig(f_prev, fullfile(outdir,'fig3_delta_scenarios_I.fig'));
    close(f_prev);
    % Plot 2: Vary mortality - daily deaths 
    % --------------------------------------
    f_deaths = figure('Name','Vary mortality - deaths');
    set(f_deaths, 'Position', [10 10 900 1100]); pause(.1)
    t_layout2 = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
    
    nexttile;
    plot(tA, FyA, '-',  'LineWidth', lw);
    hold on;
    plot(tB, FyB, '--', 'LineWidth', lw);
    xlim([0 365]); 
    ylim([0 1.5e-5]);
    ylabel({'Daily death rate','(young group)'});
    set(gca,'FontSize',fs_all);
    
    nexttile;
    plot(tA, FoA, '-',  'LineWidth', lw); 
    hold on;
    plot(tB, FoB, '--', 'LineWidth', lw);
    xlim([0 365]);
    ylim([0 1.5e-5]);
    ylabel({'Daily death rate','(old group)'});
    set(gca,'FontSize',fs_all);
    
    nexttile;
    plot(tA, fA, '-',  'LineWidth', lw); 
    hold on;
    plot(tB, fB, '--', 'LineWidth', lw);
    xlabel('Days');
    ylabel({'Daily death rate','(total)'});
    xlim([0 365]); 
    ylim([0 1.5e-5]);
    set(gca,'FontSize',fs_all);
    
    savepng(fullfile(outdir,'fig3_delta_scenarios_daily_deaths.png'));
    savefig(f_deaths, fullfile(outdir,'fig3_delta_scenarios_daily_deaths.fig'));
    close(f_deaths);
end

%% ========= Theta sweep + totals for infections & deaths, multi-alpha ===
if type.fig4 ==1
    % For each alpha in the list, and for each Theta:
    %   - keep r_y*delta_y + (1-r_y)*delta_o = delta0 using delta_from_theta()
    %   - simulate 2-group model to day 365
    %   - record R_y(365), R_o(365), R_total(365)
    %   - compute total deaths as time-integrals of F_y and F_o:
    %         F_y(t) = delta_y * I_y(t) / tauI,   F_o(t) = delta_o * I_o(t) / tauI
    %         Death_y = \int_0^{365} F_y(t) dt,   Death_o = \int_0^{365} F_o(t) dt

    fprintf('\n=== Fig 4: Theta sweep + totals for infections & deaths, multi-alpha ===\n')

    ry_list     = 0.5;   % We can extend to [0.3 0.5 0.7]

    for ry = ry_list
        ro = 1 - ry;

        for a = alpha_list
            % storage
            Ry365 = zeros(size(Theta_list));
            Ro365 = zeros(size(Theta_list));
            Rtot  = zeros(size(Theta_list));
            Dy    = zeros(size(Theta_list));
            Do    = zeros(size(Theta_list));
            Dtot  = zeros(size(Theta_list));

            for k = 1:numel(Theta_list)
                Theta = Theta_list(k);

                % fair normalization for deltas
                [dy, do] = delta_from_theta(Theta, ry, delta0);

                % parameters for this (alpha, Theta)
                pT = p2;                    % start from your base two-group struct
                pT.alpha    = a;
                pT.ry       = ry;
                pT.ro       = ro;
                pT.delta_y  = dy;
                pT.delta_o  = do;

                % ICs proportional to population shares (conserve total I0)
                Iy0 = ry*I0; Io0 = ro*I0;
                Sy0 = ry - Iy0; So0 = ro - Io0;
                y0_2 = [Sy0; 0; Iy0; 0; 0;   So0; 0; Io0; 0; 0];

                % simulate to the shared time grid
                [tT, yT] = ode45(@(t,y) rhs_twogroup(t,y,pT), tspan, y0_2, opts2);

                % cumulative removals (infections) at day 365
                Ry365(k) = yT(end,4);
                Ro365(k) = yT(end,9);
                Rtot(k)  = Ry365(k) + Ro365(k);

                % deaths over time: integrate F_y and F_o
                Iy_t = yT(:,3); Io_t = yT(:,8);
                Fy_t = dy .* Iy_t ./ tauI;      % F_y(t)
                Fo_t = do .* Io_t ./ tauI;      % F_o(t)
                Dy(k) = trapz(tT, Fy_t);
                Do(k) = trapz(tT, Fo_t);
                Dtot(k) = Dy(k) + Do(k);
            end

            % ---- plots (log-x)
            aTag = sprintf('alpha_%g_ry%02d', a, round(100*ry));  % unique tag for filenames

            % Infections
            figure;
            set(gcf, 'Position', [10 10 900 600]);pause(.1)
            plot(Theta_list, Ry365, 's--', 'LineWidth',1.8, 'MarkerSize',6, 'DisplayName','young'); hold on;
            plot(Theta_list, Ro365, 'o-',  'LineWidth',1.8, 'MarkerSize',6, 'DisplayName','old');
            plot(Theta_list, Rtot,  '^-',  'LineWidth',2.0, 'MarkerSize',6, 'DisplayName','total');
            set(gca,'XScale','log');
            xlabel('\theta = \delta_y / \delta_o (log scale)');
            ylabel('Cumulative infections at day 365');
            ylim([0 1]);
            set(gca,'FontSize',fs_all);
            title(sprintf('Infections vs \\theta  (\\alpha=%g)', a));
            if a==400000
                legend('Location','best');
            end

            savepng(fullfile(outdir, sprintf('fig4_cases_vs_Theta_%s.png', aTag)));
            savefig(fullfile(outdir, sprintf('fig4_cases_vs_Theta_%s.fig', aTag)));
            close;

            % Deaths
            figure;
            set(gcf, 'Position', [10 10 900 600]); pause(.1)
            plot(Theta_list, Dy, 's--', 'LineWidth',1.8, 'MarkerSize',6, 'DisplayName','young'); hold on;
            plot(Theta_list, Do, 'o-',  'LineWidth',1.8, 'MarkerSize',6, 'DisplayName','old');
            plot(Theta_list, Dtot,'^-',  'LineWidth',2.0, 'MarkerSize',6, 'DisplayName','total');
            set(gca,'XScale','log');
            ylim([0 5e-3])
            xlabel('\theta = \delta_y / \delta_o (log scale)');
            ylabel('Cumulative deaths at day 365')
            title(sprintf(['Deaths vs \' ...
                '\theta  (\\alpha=%g)'],  a));
            if a ==400000
                legend('Location','best');
            end
            set(gca,'FontSize',fs_all);
            savepng(fullfile(outdir, sprintf('fig4_deaths_vs_Theta_%s.png', aTag)));
            savefig(fullfile(outdir, sprintf('fig4_deaths_vs_Theta_%s.fig', aTag)));
            close;

        end
    end

end

%% ========= Susceptibility heterogeneity (beta_y != beta_o) =========
if type.fig5==1
    % theta = beta_y / beta_o
    % Keep equal-weight average fixed at beta0.
    % For each alpha and theta, simulate to day 365; compare prevalence dynamics,
    % total infections R_tot(365), and deaths D_tot = \int_0^{365} F_y + F_o dt.

    fprintf('\n=== Fig 5: Susceptibility heterogeneity (theta = beta_y/beta_o) ===\n');

    ry = 0.5; ro = 1 - ry;
    Iy0 = ry*I0; Io0 = ro*I0;
    Sy0 = ry - Iy0; So0 = ro - Io0;
    y0_2 = [Sy0; 0; Iy0; 0; 0;   So0; 0; Io0; 0; 0];

    rho_sus = [1, 0.5, 0.25, 0.1, 0.05, 0.02, 0.01, 1e-3, 1e-4, 0];
    alpha_sus = [0, 200000, 400000];

    for a = alpha_sus
        % storage
        Rtot = zeros(size(rho_sus));
        Ry365= zeros(size(rho_sus));
        Ro365= zeros(size(rho_sus));
        Dtot = zeros(size(rho_sus));
        Dy   = zeros(size(rho_sus));
        Do   = zeros(size(rho_sus));

        for k = 1:numel(rho_sus)
            rho = rho_sus(k);

            % Map theta -> (beta_y, beta_o) with average fixed at beta0
            [beta_y, beta_o] = beta_from_rho(rho, ry, beta0);

            % build parameter struct
            pS = p2;
            pS.alpha   = a;
            pS.ry      = ry;  pS.ro = ro;
            pS.beta_y  = beta_y;
            pS.beta_o  = beta_o;

            % simulate
            [tS, yS] = ode45(@(t,y) rhs_twogroup(t,y,pS), tspan, y0_2, opts2);

            % prevalence to plot
            IyS = yS(:,3); IoS = yS(:,8); ItS = IyS + IoS;

            % totals at day 365 (infections)
            Ry365(k) = yS(end,4);
            Ro365(k) = yS(end,9);
            Rtot(k)  = Ry365(k) + Ro365(k);

            % deaths over time: integrate F_y and F_o
            Iy_t = yS(:,3); Io_t = yS(:,8);
            Fy_t = pS.delta_y .* Iy_t ./ tauI;   % delta_y from p2 (unchanged)
            Fo_t = pS.delta_o .* Io_t ./ tauI;
            Dy(k) = trapz(tS, Fy_t);
            Do(k) = trapz(tS, Fo_t);
            Dtot(k) = Dy(k) + Do(k);
        end

        % Infections vs Theta
        figure;
        set(gcf, 'Position', [10 10 900 600]);pause(.1)
        plot(rho_sus, Ry365, 's--','LineWidth',1.8,'MarkerSize',6,'DisplayName','young'); hold on;
        plot(rho_sus, Ro365, 'o-','LineWidth',1.8,'MarkerSize',6,'DisplayName','old');
        plot(rho_sus, Rtot,  '^-','LineWidth',2.0,'MarkerSize',6,'DisplayName','total');
        set(gca,'XScale','log');
        xlabel('\rho= \beta_y / \beta_o (log scale)'); ylabel('Cumulative infections at day 365');ylim([0 1]);
        title(sprintf('Infections vs \\rho  (\\alpha=%g)', a));
        set(gca,'fontsize',fs_all)
        if a ==0
            legend('Location','best');
        end
        savepng(fullfile(outdir, sprintf('fig5_cases_vs_theta_alpha_%g.png', a)));
        savefig(fullfile(outdir, sprintf('fig5_cases_vs_theta_alpha_%g.fig', a)));
        close;

        % Deaths vs Theta
        figure;
        set(gcf, 'Position', [10 10 900 600]);pause(.1)
        plot(rho_sus, Dy, 's--','LineWidth',1.8,'MarkerSize',6,'DisplayName','young'); hold on;
        plot(rho_sus, Do, 'o-','LineWidth',1.8,'MarkerSize',6,'DisplayName','old');
        plot(rho_sus, Dtot, '^-','LineWidth',2.0,'MarkerSize',6,'DisplayName','total');
        set(gca,'XScale','log');
        ylim([0 5e-3])
        xlabel('\rho = \beta_y / \beta_o (log scale)'); ylabel('Cumulative deaths at day 365');
        title(sprintf('Deaths vs \\rho  (\\alpha=%g)', a));
        set(gca,'fontsize',fs_all)
        if a==0
            legend('Location','best');
        end
        savepng(fullfile(outdir, sprintf('fig5_deaths_vs_theta_alpha_%g.png', a)));
        savefig(fullfile(outdir, sprintf('fig5_deaths_vs_theta_alpha_%g.fig', a)));
        close;
    end
end

%% ========= IFR heterogeneity with heterogeneous C =================
% C_base = [0.5 0.5; 0.5 0.5] &  C_het  = [0.9 0.1; 0.1 0.9]
if type.fig6 ==1

    fprintf('\n=== Fig 6: IFR heterogeneity with heterogeneous C ===\n')

    C_alt  = C_het;
    % population split
    ry = 0.5;
    ro = 1 - ry;
    % initial conditions consistent with ry, ro
    Iy0 = ry*I0; Io0 = ro*I0;
    Sy0 = ry - Iy0; So0 = ro - Io0;
    y0_2 = [Sy0; 0; Iy0; 0; 0;   So0; 0; Io0; 0; 0];

    % storage: rows = alpha, cols = Theta
    Dtot_base = zeros(numel(alpha_list), numel(Theta_list));
    Dtot_alt  = zeros(numel(alpha_list), numel(Theta_list));
    for ia = 1:numel(alpha_list)
        a = alpha_list(ia);

        for k = 1:numel(Theta_list)
            Theta = Theta_list(k);

            % fair normalization for deltas
            [dy, do] = delta_from_theta(Theta, ry, delta0);

            % Run with C_base
            pB = p2;
            pB.alpha   = a;
            pB.C       = C_base;
            pB.ry      = ry;   pB.ro = ro;
            pB.delta_y = dy;   pB.delta_o = do;
            pB.beta_y  = beta0; pB.beta_o = beta0;

            [tB, yB] = ode45(@(t,y) rhs_twogroup(t,y,pB), tspan, y0_2, opts2);
            Rtot_base(ia,k) = yB(end,4) + yB(end,9);   % R_y + R_o

            IyB = yB(:,3);  IoB = yB(:,8);
            FyB = dy .* IyB ./ tauI;
            FoB = do .* IoB ./ tauI;
            Dtot_base(ia,k) = trapz(tB, FyB + FoB);

            % Run with C_alt
            pA = p2;
            pA.alpha   = a;
            pA.C       = C_alt;
            pA.ry      = ry;   pA.ro = ro;
            pA.delta_y = dy;   pA.delta_o = do;
            pA.beta_y  = beta0; pA.beta_o = beta0;

            [tA, yA] = ode45(@(t,y) rhs_twogroup(t,y,pA), tspan, y0_2, opts2);
            Rtot_alt(ia,k) = yA(end,4) + yA(end,9);

            IyA = yA(:,3);  IoA = yA(:,8);
            FyA = dy .* IyA ./ tauI;
            FoA = do .* IoA ./ tauI;
            Dtot_alt(ia,k) = trapz(tA, FyA + FoA);
        end
    end

    % Combined Plot: Infections: R_tot
    for ia = 1:numel(alpha_list)
        figure;
        set(gcf, 'Position', [10 10 900 600]); pause(.1)

        a = alpha_list(ia);
        plot(Theta_list, Rtot_base(ia,:), '-o', 'LineWidth',1.8, 'MarkerSize',5, ...
            'DisplayName','C base [0.5 0.5; 0.5 0.5]'); hold on;
        plot(Theta_list, Rtot_alt(ia,:),  '--s','LineWidth',1.8, 'MarkerSize',5, ...
            'DisplayName','C alt [0.9 0.1; 0.1 0.9]');
        set(gca,'XScale','log');
        ylim([0 1])
        ylabel('Cumulative infections at day 365')
        title(sprintf('Infections vs \\theta (\\alpha = %g)', a));
        set(gca,'fontsize',fs_all)

        xlabel('\theta = \delta_y/\delta_o (log scale)');

        if a ==0
            legend('Location','best');
        end

        savepng(fullfile(outdir,sprintf('fig6_Rtot_cases_%d.png',a)));
        savefig(gcf, fullfile(outdir,sprintf('fig6_Rtot_cases_%d.fig',a)));
        close;
    end

    % Combined Plot: Deaths: D_tot
    for ia = 1:numel(alpha_list)
        figure;
        set(gcf, 'Position', [10 10 900 600]); pause(.1)

        a = alpha_list(ia);
        plot(Theta_list, Dtot_base(ia,:), '-o', 'LineWidth',1.8, 'MarkerSize',5, ...
            'DisplayName','C base [0.5 0.5; 0.5 0.5]'); hold on;
        plot(Theta_list, Dtot_alt(ia,:),  '--s','LineWidth',1.8, 'MarkerSize',5, ...
            'DisplayName','C alt [0.9 0.1; 0.1 0.9]');
        set(gca,'XScale','log');
        ylim([0 5e-3])
        ylabel('Cumulative deaths at day 365')
        title(sprintf('Deaths vs \\theta (\\alpha = %g)', a));
        set(gca,'fontsize',fs_all)

        xlabel('\theta = \delta_y/\delta_o (log scale)');
        if a==0
            legend('Location','best');
        end

        savepng(fullfile(outdir,sprintf('fig6_Dtot_deaths_%d.png',a)));
        savefig(gcf, fullfile(outdir,sprintf('fig6_Dtot_deaths_%d.fig',a)));
        close;
    end
end

%% ========= Susceptibility heterogeneity with heterogeneous C =========
% C_base = [0.5 0.5; 0.5 0.5] & C_het  = [0.9 0.1; 0.1 0.9]
if type.fig7==1

    fprintf('\n=== Fig 7: Susceptibility heterogeneity with heterogeneous C  ===\n')
    C_alt  = C_het;

    ry = 0.5;
    ro = 1 - ry;
    % ICs consistent with ry, ro
    Iy0 = ry*I0; Io0 = ro*I0;
    Sy0 = ry - Iy0; So0 = ro - Io0;
    y0_2 = [Sy0; 0; Iy0; 0; 0;   So0; 0; Io0; 0; 0];
    Dtot_base = zeros(numel(alpha_list), numel(rho_list));
    Dtot_het  = zeros(numel(alpha_list), numel(rho_list));
    Rtot_base = zeros(numel(alpha_list), numel(rho_list));
    Rtot_het  = zeros(numel(alpha_list), numel(rho_list));
    for ia = 1:numel(alpha_list)
        a = alpha_list(ia);

        for k = 1:numel(rho_list)
            rho = rho_list(k);
            [beta_y, beta_o] = beta_from_rho(rho, ry, beta0);

            % base C
            pB = p2;
            pB.alpha  = a;
            pB.C      = C_base;
            pB.ry     = ry;  pB.ro = ro;
            pB.beta_y = beta_y;
            pB.beta_o = beta_o;

            [tB, yB] = ode45(@(t,y) rhs_twogroup(t,y,pB), tspan, y0_2, opts2);
            Rtot_base(ia,k) = yB(end,4) + yB(end,9);
            IyB = yB(:,3); IoB = yB(:,8);
            FyB = pB.delta_y .* IyB ./ tauI;
            FoB = pB.delta_o .* IoB ./ tauI;
            Dtot_base(ia,k) = trapz(tB, FyB + FoB);

            %C_het
            pH = p2;
            pH.alpha  = a;
            pH.C      = C_alt;
            pH.ry     = ry;  pH.ro = ro;
            pH.beta_y = beta_y;
            pH.beta_o = beta_o;

            [tH, yH] = ode45(@(t,y) rhs_twogroup(t,y,pH), tspan, y0_2, opts2);
            Rtot_het(ia,k) = yH(end,4) + yH(end,9);


            IyH = yH(:,3); IoH = yH(:,8);
            FyH = pH.delta_y .* IyH ./ tauI;
            FoH = pH.delta_o .* IoH ./ tauI;
            Dtot_het(ia,k) = trapz(tH, FyH + FoH);

        end
    end


    % Combined Plot: Infections: R_tot
    % Goal: Compare R_tot(365) vs theta = beta_y/beta_o under two contact matrices.
    % C_base = [0.5 0.5; 0.5 0.5] & C_het  = [0.9 0.1; 0.1 0.9]
    for ia = 1:numel(alpha_list)
        figure;
        set(gcf, 'Position', [10 10 900 600]); pause(.1)
        a = alpha_list(ia);
        plot(rho_list, Rtot_base(ia,:), '-o', 'LineWidth',1.8, 'MarkerSize',5, ...
            'DisplayName','C base [0.5 0.5; 0.5 0.5]'); hold on;
        plot(rho_list, Rtot_het(ia,:),  '--s','LineWidth',1.8, 'MarkerSize',5, ...
            'DisplayName','C alt [0.9 0.1; 0.1 0.9]');
        set(gca,'XScale','log');
        ylim([0 1])
        ylabel('Cumulatives infections at day 365')
        title(sprintf('Infections vs \\rho (\\alpha = %g)', a));
        set(gca,'fontsize',fs_all)
        xlabel('\rho = \beta_y/\beta_o (log scale)');

        if a ==0
            legend('Location','best');
        end

        savepng(fullfile(outdir,sprintf('fig7_Rtot_cases_%d.png',a)));
        savefig(gcf, fullfile(outdir,sprintf('fig7_Rtot_cases_%d.fig',a)));
        close;
    end

    % Combined Plot: Deaths: D_tot
    % Compare total deaths D_tot = ∫ (F_y + F_o) dt vs theta under two C matrices.
    % C_base  & C_het
    for ia = 1:numel(alpha_list)
        figure;
        set(gcf, 'Position', [10 10 900 600]);pause(.1)
        a = alpha_list(ia);

        plot(rho_list, Dtot_base(ia,:), '-o', 'LineWidth',1.8, 'MarkerSize',5, ...
            'DisplayName','C base [0.5 0.5; 0.5 0.5]'); hold on;
        plot(rho_list, Dtot_het(ia,:),  '--s','LineWidth',1.8, 'MarkerSize',5, ...
            'DisplayName','C alt  [0.9 0.1; 0.1 0.9]');
        set(gca,'XScale','log');
        ylabel('Cumulative deaths at day 365')
        title(sprintf('Deaths vs \\rho (\\alpha = %g)', a));
        ylim([0 5e-3])
        set(gca,'fontsize',fs_all)

        xlabel('\rho = \beta_y/\beta_o (log scale)');

        if a==0
            legend('Location','best');
        end

        savepng(fullfile(outdir,sprintf('fig7_Dtot_deaths_%d.png',a)));
        savefig(gcf, fullfile(outdir,sprintf('fig7_Dtot_deaths_%d.fig',a)));
        close;
    end
end
%% ========= Appendix Section ==================================
%% =========Appendix A: Continuous Alpha Sweep

if type.figA1 == 1
    fprintf('\n=== Fig A1 (Appendix A): Vary Behavior Responsiveness - alpha (0 to 10^6) ===\n')
    alpha_sweep = [linspace(0, 1e6, 100)]; 
    cases_365  = zeros(size(alpha_sweep));
    deaths_365 = zeros(size(alpha_sweep));
    ry = 0.5; 
    ro = 0.5;
    dy = delta0;
    do = delta0;
    
    for k = 1:numel(alpha_sweep)
        pA = p2;
        pA.alpha   = alpha_sweep(k);
        pA.ry      = ry; 
        pA.ro = ro; 
        pA.delta_y = dy; 
        pA.delta_o = do;
        Iy0 = ry*I0;
        Io0 = ro*I0; Sy0 = ry - Iy0; So0 = ro - Io0;
        y0_2 = [Sy0; 0; Iy0; 0; 0;  So0; 0; Io0; 0; 0];
        [tT, yT] = ode45(@(t,y) rhs_twogroup(t,y,pA), tspan, y0_2, opts2);
        cases_365(k) = yT(end,4) + yT(end,9);
        Iy_t = yT(:,3); Io_t = yT(:,8);
        Fy_t = dy .* Iy_t ./ tauI; Fo_t = do .* Io_t ./ tauI;
        deaths_365(k) = trapz(tT, Fy_t + Fo_t);
    end
    
    % Figure 1: Cumulative Infections vs Alpha
    figure('Name','Appendix_Alpha_Sweep_Cases');
    f_cases = figure(1); 
    set(f_cases, 'Position', [10 10 800 600]); pause(.1)
    
    plot(alpha_sweep, cases_365, '-', 'LineWidth', lw, 'Color', [0 0.4470 0.7410]); 
    hold on;
    ylabel('Cumulative infections at day 365'); 
    xlabel('\alpha (Behavioral responsiveness)'); 
    ylim([0 1]);
    grid off;
    set(gca, 'FontSize', fs_all);
    target_alphas = [0, 200000, 400000];
    colors = ['g', 'm', 'c']; 
    for idx = 1:numel(target_alphas)
        [~, closest_idx] = min(abs(alpha_sweep - target_alphas(idx)));
        plot(alpha_sweep(closest_idx), cases_365(closest_idx), 'o', 'MarkerFaceColor', colors(idx), 'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'LineWidth', 1.5);
        text(alpha_sweep(closest_idx), cases_365(closest_idx) + 0.03, sprintf(' \\alpha = %g', target_alphas(idx)), 'FontSize', fs_legend, 'FontWeight', 'bold');
    end
    
    savepng(fullfile(outdir, 'figA1_cases_vs_alpha.png')); 
    savefig(gcf, fullfile(outdir, 'figA1_cases_vs_alpha.fig'));
    close;
    
    % Figure 2: Cumulative Deaths vs Alpha
    figure('Name','Appendix_Alpha_Sweep_Deaths'); 
    f_deaths = figure(1); 
    set(f_deaths, 'Position', [10 10 800 600]); pause(.1)
    
    plot(alpha_sweep, deaths_365, '-', 'LineWidth', lw, 'Color', [0.8500 0.3250 0.0980]); 
    hold on;
    ylabel('Cumulative deaths at day 365'); 
    xlabel('\alpha (Behavioral responsiveness)'); 
    ylim([0 5e-3]);
    grid off; 
    set(gca, 'FontSize', fs_all);
    for idx = 1:numel(target_alphas)
        [~, closest_idx] = min(abs(alpha_sweep - target_alphas(idx)));
        plot(alpha_sweep(closest_idx), deaths_365(closest_idx), 's', 'MarkerFaceColor', colors(idx), 'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'LineWidth', 1.5);
        text(alpha_sweep(closest_idx), deaths_365(closest_idx) + 1.5e-4, sprintf(' \\alpha = %g', target_alphas(idx)), 'FontSize', fs_legend, 'FontWeight', 'bold');
    end
    savepng(fullfile(outdir, 'figA1_deaths_vs_alpha.png')); 
    savefig(gcf, fullfile(outdir, 'figA1_deaths_vs_alpha.fig')); 
    close;
end
%% ========= Appendix B: Susceptibility heterogeneity with heterogeneous C (multiple choices)
if type.figB1==1
    fprintf('\n=== Fig B1 (Appendix B): Susceptibility heterogeneity with heterogeneous C ===\n')

    p_vals = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
    alpha_list_test = [0, 200000, 400000];

    ry = 0.5; ro = 1 - ry;
    Iy0 = ry*I0; Io0 = ro*I0;
    Sy0 = ry - Iy0; So0 = ro - Io0;
    y0_2 = [Sy0; 0; Iy0; 0; 0;   So0; 0; Io0; 0; 0];


    for ia = 1:numel(alpha_list_test)
        a = alpha_list_test(ia);

        Ry_all_C = zeros(numel(p_vals), numel(rho_list)); % Infections Young
        Ro_all_C = zeros(numel(p_vals), numel(rho_list)); % Infections Old

        for ip = 1:numel(p_vals)
            p_val = p_vals(ip);
            C_current = [p_val, 1-p_val; 1-p_val, p_val];

            for k = 1:numel(rho_list)
                rho = rho_list(k);
                [beta_y, beta_o] = beta_from_rho(rho, ry, beta0);

                pS = p2;
                pS.alpha  = a;
                pS.C      = C_current;
                pS.ry     = ry; pS.ro = ro;
                pS.beta_y = beta_y; pS.beta_o = beta_o;

                [tS, yS] = ode45(@(t,y) rhs_twogroup(t,y,pS), tspan, y0_2, opts2);
                Ry_all_C(ip, k) = yS(end,4) + yS(end,5);
                Ro_all_C(ip, k) = yS(end,9) + yS(end,10);


                Iy_t = yS(:,3); Io_t = yS(:,8);
                Fy_t = pS.delta_y .* Iy_t ./ tauI;
                Fo_t = pS.delta_o .* Io_t ./ tauI;
            end
        end

        f1 = figure();
        set(gcf, 'Position', [10 10 900 600]); pause(.1)
        hold on;
        for ip = 1:numel(p_vals), plot(rho_list, Ry_all_C(ip,:), '-o', 'LineWidth', 1.8, 'MarkerSize', 5,'DisplayName', sprintf('p=%.1f', p_vals(ip))); end
        set(gca, 'XScale', 'log');
        ylim([0 0.5]);
        ylabel({'Cumulative infections at day 365','(young group)'});
        xlabel('\rho = \beta_y/\beta_o (log scale)');
        title(sprintf('Infections in young group vs \\rho (\\alpha = %d)', a)); grid off;
        box on;
        set(gca, 'FontSize', fs_all);
        if ia == 1, legend('Location', 'northwest'); end
        savepng(fullfile(outdir,sprintf('figB1_infections_YOUNG_vs_rho_%d.png',a)));
        savefig(f1,fullfile(outdir,sprintf('figB1_infections_YOUNG_vs_rho_%d.fig',a)));
        close;
        f2 = figure();
        set(gcf, 'Position', [10 10 900 600]); pause(.1)
        hold on;
        for ip = 1:numel(p_vals), plot(rho_list, Ro_all_C(ip,:), '-s', 'LineWidth', 1.8, 'MarkerSize', 5,'DisplayName', sprintf('p=%.1f', p_vals(ip))); end
        set(gca, 'XScale', 'log');
        ylabel({'Cumulative infections at day 365','(old group)'});
        xlabel('\rho = \beta_y/\beta_o (log scale)');
        ylim([0 0.5]); title(sprintf('Infections in old group vs \\rho (\\alpha = %d)', a)); grid off;
        box on;
        set(gca, 'FontSize', fs_all);

        savepng(fullfile(outdir,sprintf('figB1_infections_OLD_vs_rho_%d.png',a)));
        savefig(f2,fullfile(outdir,sprintf('figB1_infections_OLD_vs_rho_%d.fig',a)));
        close;
    end


end

%% ================Appendix C =============
%% ========= Appendix C1: Fractional Alpha Sweep 
if type.figC1 == 1
    fprintf('\n=== Fig C1 (Appendix C): Vary Behavior Responsiveness (alpha) with fractional form ===\n')
    
    alpha_sweep_frac = linspace(0, 1e6, 100); 
    cases_365_frac  = zeros(size(alpha_sweep_frac));
    deaths_365_frac = zeros(size(alpha_sweep_frac));
    
    ry = 0.5;
    ro = 0.5;
    dy = delta0;
    do = delta0;
    
    for k = 1:numel(alpha_sweep_frac)
        pA = p2;
        pA.alpha   = alpha_sweep_frac(k);
        pA.ry      = ry;
        pA.ro = ro; 
        pA.delta_y = dy; 
        pA.delta_o = do;
        pA.use_fractional = 1; % Activates the non-linear fractional behavioral feedback formulation (Equation C.1)
        
        Iy0 = ry*I0;
        Io0 = ro*I0;
        Sy0 = ry - Iy0;
        So0 = ro - Io0;
        y0_2 = [Sy0; 0; Iy0; 0; 0;  So0; 0; Io0; 0; 0];
       
        [tT, yT] = ode45(@(t,y) rhs_twogroup(t,y,pA), tspan, y0_2, opts2);
        cases_365_frac(k) = yT(end,4) + yT(end,9);
        Iy_t = yT(:,3);
        Io_t = yT(:,8);
        Fy_t = dy .* Iy_t ./ tauI; 
        Fo_t = do .* Io_t ./ tauI;
        deaths_365_frac(k) = trapz(tT, Fy_t + Fo_t);
    end
    
    target_alphas_frac = [0, 200000, 400000]; 
    colors_frac        = ['g', 'm', 'c']; 
    
  
    % Plot Separated Figure 1: Cases

    figure; set(gcf, 'Position', [10 10 800 600]); pause(.1)
    
    plot(alpha_sweep_frac, cases_365_frac, '-', 'LineWidth', lw, 'Color', [0 0.4470 0.7410]);
    hold on;
    ylabel('Cumulative infections at day 365'); 
    xlabel('\alpha (Fractional)');
    ylim([0 1]);
    set(gca, 'FontSize', fs_all);
    
    for idx = 1:numel(target_alphas_frac)
        [~, closest_idx] = min(abs(alpha_sweep_frac - target_alphas_frac(idx)));
        a_val = alpha_sweep_frac(closest_idx);
        c_val = cases_365_frac(closest_idx);
        
        plot(a_val, c_val, 'o', 'MarkerFaceColor', colors_frac(idx), 'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'LineWidth', 1.5);
        text(a_val, c_val + 0.03, sprintf(' \\alpha = %g', target_alphas_frac(idx)), 'FontSize', fs_legend, 'FontWeight', 'bold');
    end
    
    savepng(fullfile(outdir, 'figC1_fractional_cases_sweep.png'));
    savefig(gcf, fullfile(outdir, 'figC1_fractional_cases_sweep.fig'));
    close;
    
   
    % Plot Separated Figure 2: Deaths

    figure; set(gcf, 'Position', [10 10 800 600]);pause(.1)
    
    plot(alpha_sweep_frac, deaths_365_frac, '-', 'LineWidth', lw, 'Color', [0.8500 0.3250 0.0980]);
    hold on;
    ylabel('Cumulative deaths at day 365'); 
    xlabel('\alpha (Fractional)');
    ylim([0 5e-3]); 
    set(gca, 'FontSize', fs_all);
    
    for idx = 1:numel(target_alphas_frac)
        [~, closest_idx] = min(abs(alpha_sweep_frac - target_alphas_frac(idx)));
        a_val = alpha_sweep_frac(closest_idx);
        d_val = deaths_365_frac(closest_idx);
        
        plot(a_val, d_val, 's', 'MarkerFaceColor', colors_frac(idx), 'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'LineWidth', 1.5);
        text(a_val, d_val + 1.5e-4, sprintf(' \\alpha = %g', target_alphas_frac(idx)), 'FontSize', fs_legend, 'FontWeight', 'bold');
    end
    
    savepng(fullfile(outdir, 'figC1_fractional_deaths_sweep.png'));
    savefig(gcf, fullfile(outdir, 'figC1_fractional_deaths_sweep.fig'));
    close;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ========= Appendix C2: Delta scenarios (group mortality) ==========================
if type.figC2==1
    fprintf('\n=== Fig C2: Delta scenarios (vary group mortality) ===\n')

    m  = 0.5; ry = m; ro = 1 - m;
    Iy0 = ry*I0; Io0 = ro*I0;
    Sy0 = ry - Iy0; So0 = ro - Io0;
    y0_2 = [Sy0; 0; Iy0; 0; 0;  So0; 0; Io0; 0; 0];

    % Homogeneous: dy=do=delta0  (m=0.5)
    pA = p2; pA.ry=ry; pA.ro=ro; pA.delta_y=delta0; pA.delta_o=delta0;

    pA.use_fractional = 1; % Force active for the homogeneous scenario

    [tA, yA] = ode45(@(t,y) rhs_twogroup(t,y,pA), tspan, y0_2, opts2);
    % infectious prevalence per group
    IyA = yA(:,3);
    IoA = yA(:,8);
    ItotA = IyA+IoA; % total infectious prevalence
    % daily death proxy per group
    FyA = pA.delta_y .* IyA ./ tauI;
    FoA = pA.delta_o .* IoA ./ tauI;
    fA  = FyA + FoA;   % total daily deaths

    % Heterogeneous: dy=0, do=2*delta0
    pB = p2; pB.ry=ry; pB.ro=ro; pB.delta_y=0; pB.delta_o=2*delta0;
    pB.use_fractional = 1; % Force active for the heterogeneous scenario

    [tB, yB] = ode45(@(t,y) rhs_twogroup(t,y,pB), tspan, y0_2, opts2);
    % infectious prevalence per group
    IyB = yB(:,3);
    IoB = yB(:,8);
    ItotB = IyB+IoB;
    % daily death proxy per group
    FyB = pB.delta_y .* IyB ./ tauI;
    FoB = pB.delta_o .* IoB ./ tauI;
    fB  = FyB + FoB;   % total daily deaths (Scenario B)
    
    % Plot 1: Vary mortality - prevalence (Upgraded to Tiledlayout)
    % -------------------------------------------------------------
    f_prev = figure('Name','Vary mortality - prevelance');
    set(f_prev, 'Position', [10 10 900 1100]); pause(.1)
    
    t_layout1 = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
    
    nexttile;
    plot(tA, IyA, '-',  'LineWidth', lw); 
    hold on;
    plot(tB, IyB, '--', 'LineWidth', lw);
    xlim([0 365]); 
    ylim([0 0.06]);
    yticks([0 0.02 0.04 0.06]);
    ylabel({'Prevalence','(young group)'});
    set(gca,'FontSize',fs_all);
    legend({'Homogeneous: \delta_y=\delta_o','Heterogeneous: \delta_y=0, \delta_o=2\delta_0'}, ...
        'Location','northeast','FontSize',fs_legend);
    
    nexttile;
    plot(tA, IoA, '-',  'LineWidth', lw);
    hold on;
    plot(tB, IoB, '--', 'LineWidth', lw);
    xlim([0 365]); 
    ylim([0 0.06]); 
    yticks([0 0.02 0.04 0.06]);
    ylabel({'Prevalence','(old group)'});
    set(gca,'FontSize',fs_all);
    
    nexttile;
    plot(tA, ItotA, '-',  'LineWidth', lw); 
    hold on;
    plot(tB, ItotB, '--', 'LineWidth', lw);
    xlabel('Days');
    ylabel({'Prevalence','(total)'});
    xlim([0 365]);
    ylim([0 0.06]); 
    yticks([0 0.02 0.04 0.06]);
    set(gca,'FontSize',fs_all);
    
    savepng(fullfile(outdir,'figC2_delta_scenarios_I.png'));
    savefig(f_prev, fullfile(outdir,'figC2_delta_scenarios_I.fig'));
    close(f_prev);
    % Plot 2: Vary mortality - daily deaths 
    % --------------------------------------
    f_deaths = figure('Name','Vary mortality - deaths');
    set(f_deaths, 'Position', [10 10 900 1100]); pause(.1)
    t_layout2 = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
    
    nexttile;
    plot(tA, FyA, '-',  'LineWidth', lw);
    hold on;
    plot(tB, FyB, '--', 'LineWidth', lw);
    xlim([0 365]); 
    ylim([0 1.5e-5]);
    ylabel({'Daily death rate','(young group)'});
    set(gca,'FontSize',fs_all);
    
    nexttile;
    plot(tA, FoA, '-',  'LineWidth', lw); 
    hold on;
    plot(tB, FoB, '--', 'LineWidth', lw);
    xlim([0 365]);
    ylim([0 1.5e-5]);
    ylabel({'Daily death rate','(old group)'});
    set(gca,'FontSize',fs_all);
    
    nexttile;
    plot(tA, fA, '-',  'LineWidth', lw); 
    hold on;
    plot(tB, fB, '--', 'LineWidth', lw);
    xlabel('Days');
    ylabel({'Daily death rate','(total)'});
    xlim([0 365]); 
    ylim([0 1.5e-5]);
    set(gca,'FontSize',fs_all);
    
    savepng(fullfile(outdir,'figC2_delta_scenarios_daily_deaths.png'));
    savefig(f_deaths, fullfile(outdir,'figC2_delta_scenarios_daily_deaths.fig'));
    close(f_deaths);
end
%%================

%% ========= Appendix C3: Theta sweep + totals (Fractional Analysis)============
if type.figC3 ==1
   

    fprintf('\n=== Fig C3: Fractional Theta sweep + totals for infections & deaths ===\n')

    ry_list     = 0.5;   % We can extend to [0.3 0.5 0.7]

    for ry = ry_list
        ro = 1 - ry;

        for a = alpha_list
            % storage
            Ry365 = zeros(size(Theta_list));
            Ro365 = zeros(size(Theta_list));
            Rtot  = zeros(size(Theta_list));
            Dy    = zeros(size(Theta_list));
            Do    = zeros(size(Theta_list));
            Dtot  = zeros(size(Theta_list));

            for k = 1:numel(Theta_list)
                Theta = Theta_list(k);

                % fair normalization for deltas
                [dy, do] = delta_from_theta(Theta, ry, delta0);

                % parameters for this (alpha, Theta)
                pT = p2;                    % start from your base two-group struct
                pT.alpha    = a;
                pT.ry       = ry;
                pT.ro       = ro;
                pT.delta_y  = dy;
                pT.delta_o  = do;
                pT.use_fractional = 1;

                % ICs proportional to population shares (conserve total I0)
                Iy0 = ry*I0; Io0 = ro*I0;
                Sy0 = ry - Iy0; So0 = ro - Io0;
                y0_2 = [Sy0; 0; Iy0; 0; 0;   So0; 0; Io0; 0; 0];

                % simulate to the shared time grid
                [tT, yT] = ode45(@(t,y) rhs_twogroup(t,y,pT), tspan, y0_2, opts2);

                % cumulative removals (infections) at day 365
                Ry365(k) = yT(end,4);
                Ro365(k) = yT(end,9);
                Rtot(k)  = Ry365(k) + Ro365(k);

                % deaths over time: integrate F_y and F_o
                Iy_t = yT(:,3); Io_t = yT(:,8);
                Fy_t = dy .* Iy_t ./ tauI;      % F_y(t)
                Fo_t = do .* Io_t ./ tauI;      % F_o(t)
                Dy(k) = trapz(tT, Fy_t);
                Do(k) = trapz(tT, Fo_t);
                Dtot(k) = Dy(k) + Do(k);
            end

            % ---- plots (log-x)
            aTag = sprintf('alpha_%g_ry%02d', a, round(100*ry));  % unique tag for filenames

            % Infections
            figure;
            set(gcf, 'Position', [10 10 900 600]);pause(.1)
            plot(Theta_list, Ry365, 's--', 'LineWidth',1.8, 'MarkerSize',6, 'DisplayName','young'); hold on;
            plot(Theta_list, Ro365, 'o-',  'LineWidth',1.8, 'MarkerSize',6, 'DisplayName','old');
            plot(Theta_list, Rtot,  '^-',  'LineWidth',2.0, 'MarkerSize',6, 'DisplayName','total');
            set(gca,'XScale','log');
            xlabel('\theta = \delta_y / \delta_o (log scale)');
            ylabel('Cumulative infections at day 365');
            ylim([0 1]);
            set(gca,'FontSize',fs_all);
            title(sprintf('Infections vs \\theta  (\\alpha=%g)', a));
            if a==400000
                legend('Location','best');
            end

            savepng(fullfile(outdir, sprintf('figC3_cases_vs_Theta_%s.png', aTag)));
            savefig(fullfile(outdir, sprintf('figC3_cases_vs_Theta_%s.fig', aTag)));
            close;

            % Deaths plot
            figure;
            set(gcf, 'Position', [10 10 900 600]); pause(.1)
            plot(Theta_list, Dy, 's--', 'LineWidth',1.8, 'MarkerSize',6, 'DisplayName','young'); hold on;
            plot(Theta_list, Do, 'o-',  'LineWidth',1.8, 'MarkerSize',6, 'DisplayName','old');
            plot(Theta_list, Dtot,'^-',  'LineWidth',2.0, 'MarkerSize',6, 'DisplayName','total');
            set(gca,'XScale','log');
            ylim([0 5e-3])
            xlabel('\theta = \delta_y / \delta_o (log scale)');
            ylabel('Cumulative deaths at day 365')
            title(sprintf(['Deaths vs \' ...
                '\theta  (\\alpha=%g)'],  a));
            if a ==400000
                legend('Location','best');
            end
            set(gca,'FontSize',fs_all);
            savepng(fullfile(outdir, sprintf('figC3_deaths_vs_Theta_%s.png', aTag)));
            savefig(fullfile(outdir, sprintf('figC3_deaths_vs_Theta_%s.fig', aTag)));
            close;

        end
    end

end

%% ========= Appendix C4: Susceptibility heterogeneity (Fractional Analysis) =========
if type.figC4==1
    

    fprintf('\n=== Fig C4: Fractional Susceptibility heterogeneity (theta = beta_y/beta_o) ===\n');

    ry = 0.5; ro = 1 - ry;
    Iy0 = ry*I0; Io0 = ro*I0;
    Sy0 = ry - Iy0; So0 = ro - Io0;
    y0_2 = [Sy0; 0; Iy0; 0; 0;   So0; 0; Io0; 0; 0];

    rho_sus = [1, 0.5, 0.25, 0.1, 0.05, 0.02, 0.01, 1e-3, 1e-4, 0];
    alpha_sus = [0, 200000, 400000];

    for a = alpha_sus
        % storage
        Rtot = zeros(size(rho_sus));
        Ry365= zeros(size(rho_sus));
        Ro365= zeros(size(rho_sus));
        Dtot = zeros(size(rho_sus));
        Dy   = zeros(size(rho_sus));
        Do   = zeros(size(rho_sus));

        for k = 1:numel(rho_sus)
            rho = rho_sus(k);

            % Map theta -> (beta_y, beta_o) with average fixed at beta0
            [beta_y, beta_o] = beta_from_rho(rho, ry, beta0);

            % build parameter struct
            pS = p2;
            pS.alpha   = a;
            pS.ry      = ry;  pS.ro = ro;
            pS.beta_y  = beta_y;
            pS.beta_o  = beta_o;
            pS.use_fractional = 1; 

            % simulate
            [tS, yS] = ode45(@(t,y) rhs_twogroup(t,y,pS), tspan, y0_2, opts2);

            % prevalence to plot
            IyS = yS(:,3); IoS = yS(:,8); ItS = IyS + IoS;

            % totals at day 365 (infections)
            Ry365(k) = yS(end,4);
            Ro365(k) = yS(end,9);
            Rtot(k)  = Ry365(k) + Ro365(k);

            % deaths over time: integrate F_y and F_o
            Iy_t = yS(:,3); Io_t = yS(:,8);
            Fy_t = pS.delta_y .* Iy_t ./ tauI;   % delta_y from p2 (unchanged)
            Fo_t = pS.delta_o .* Io_t ./ tauI;
            Dy(k) = trapz(tS, Fy_t);
            Do(k) = trapz(tS, Fo_t);
            Dtot(k) = Dy(k) + Do(k);
        end

        % Infections vs Theta
        figure;
        set(gcf, 'Position', [10 10 900 600]);pause(.1)
        plot(rho_sus, Ry365, 's--','LineWidth',1.8,'MarkerSize',6,'DisplayName','young'); hold on;
        plot(rho_sus, Ro365, 'o-','LineWidth',1.8,'MarkerSize',6,'DisplayName','old');
        plot(rho_sus, Rtot,  '^-','LineWidth',2.0,'MarkerSize',6,'DisplayName','total');
        set(gca,'XScale','log');
        xlabel('\rho= \beta_y / \beta_o (log scale)'); ylabel('Cumulative infections at day 365');ylim([0 1]);
        title(sprintf('Infections vs \\rho  (\\alpha=%g)', a));
        set(gca,'fontsize',fs_all)
        if a ==0
            legend('Location','best');
        end
        savepng(fullfile(outdir, sprintf('figC4_cases_vs_theta_alpha_%g.png', a)));
        savefig(fullfile(outdir, sprintf('figC4_cases_vs_theta_alpha_%g.fig', a)));
        close;

        % Deaths vs Theta
        figure;
        set(gcf, 'Position', [10 10 900 600]);pause(.1)
        plot(rho_sus, Dy, 's--','LineWidth',1.8,'MarkerSize',6,'DisplayName','young'); hold on;
        plot(rho_sus, Do, 'o-','LineWidth',1.8,'MarkerSize',6,'DisplayName','old');
        plot(rho_sus, Dtot, '^-','LineWidth',2.0,'MarkerSize',6,'DisplayName','total');
        set(gca,'XScale','log');
        ylim([0 5e-3])
        xlabel('\rho = \beta_y / \beta_o (log scale)'); ylabel('Cumulative deaths at day 365');
        title(sprintf('Deaths vs \\rho  (\\alpha=%g)', a));
        set(gca,'fontsize',fs_all)
        if a==0
            legend('Location','best');
        end
        savepng(fullfile(outdir, sprintf('figC4_deaths_vs_theta_alpha_%g.png', a)));
        savefig(fullfile(outdir, sprintf('figC4_deaths_vs_theta_alpha_%g.fig', a)));
        close;
    end
end

%% ========= Appendix C5: IFR heterogeneity with heterogeneous C (Fractional Analysis) =================
% C_base = [0.5 0.5; 0.5 0.5] &  C_het  = [0.9 0.1; 0.1 0.9]
if type.figC5 ==1

    fprintf('\n=== Fig C5: Fractional IFR heterogeneity with heterogeneous C ===\n')

    C_alt  = C_het;
    % population split
    ry = 0.5;
    ro = 1 - ry;
    % initial conditions consistent with ry, ro
    Iy0 = ry*I0; Io0 = ro*I0;
    Sy0 = ry - Iy0; So0 = ro - Io0;
    y0_2 = [Sy0; 0; Iy0; 0; 0;   So0; 0; Io0; 0; 0];

    % storage: rows = alpha, cols = Theta
    Dtot_base = zeros(numel(alpha_list), numel(Theta_list));
    Dtot_alt  = zeros(numel(alpha_list), numel(Theta_list));
    for ia = 1:numel(alpha_list)
        a = alpha_list(ia);

        for k = 1:numel(Theta_list)
            Theta = Theta_list(k);

            % fair normalization for deltas
            [dy, do] = delta_from_theta(Theta, ry, delta0);

            % Run with C_base
            pB = p2;
            pB.alpha   = a;
            pB.C       = C_base;
            pB.ry      = ry;   pB.ro = ro;
            pB.delta_y = dy;   pB.delta_o = do;
            pB.beta_y  = beta0; pB.beta_o = beta0;
            pB.use_fractional = 1; 

            [tB, yB] = ode45(@(t,y) rhs_twogroup(t,y,pB), tspan, y0_2, opts2);
            Rtot_base(ia,k) = yB(end,4) + yB(end,9);   % R_y + R_o

            IyB = yB(:,3);  IoB = yB(:,8);
            FyB = dy .* IyB ./ tauI;
            FoB = do .* IoB ./ tauI;
            Dtot_base(ia,k) = trapz(tB, FyB + FoB);

            % Run with C_alt
            pA = p2;
            pA.alpha   = a;
            pA.C       = C_alt;
            pA.ry      = ry;   pA.ro = ro;
            pA.delta_y = dy;   pA.delta_o = do;
            pA.beta_y  = beta0; pA.beta_o = beta0;
            pA.use_fractional = 1;

            [tA, yA] = ode45(@(t,y) rhs_twogroup(t,y,pA), tspan, y0_2, opts2);
            Rtot_alt(ia,k) = yA(end,4) + yA(end,9);

            IyA = yA(:,3);  IoA = yA(:,8);
            FyA = dy .* IyA ./ tauI;
            FoA = do .* IoA ./ tauI;
            Dtot_alt(ia,k) = trapz(tA, FyA + FoA);
        end
    end

    % Combined Plot: Infections: R_tot
    for ia = 1:numel(alpha_list)
        figure;
        set(gcf, 'Position', [10 10 900 600]); pause(.1)

        a = alpha_list(ia);
        plot(Theta_list, Rtot_base(ia,:), '-o', 'LineWidth',1.8, 'MarkerSize',5, ...
            'DisplayName','C base [0.5 0.5; 0.5 0.5]'); hold on;
        plot(Theta_list, Rtot_alt(ia,:),  '--s','LineWidth',1.8, 'MarkerSize',5, ...
            'DisplayName','C alt [0.9 0.1; 0.1 0.9]');
        set(gca,'XScale','log');
        ylim([0 1])
        ylabel('Cumulative infections at day 365')
        title(sprintf('Infections vs \\theta (\\alpha = %g)', a));
        set(gca,'fontsize',fs_all)

        xlabel('\theta = \delta_y/\delta_o (log scale)');

        if a ==0
            legend('Location','best');
        end

        savepng(fullfile(outdir,sprintf('figC5_Rtot_cases_%d.png',a)));
        savefig(gcf, fullfile(outdir,sprintf('figC5_Rtot_cases_%d.fig',a)));
        close;
    end

    % Combined Plot: Deaths: D_tot
    for ia = 1:numel(alpha_list)
        figure;
        set(gcf, 'Position', [10 10 900 600]); pause(.1)

        a = alpha_list(ia);
        plot(Theta_list, Dtot_base(ia,:), '-o', 'LineWidth',1.8, 'MarkerSize',5, ...
            'DisplayName','C base [0.5 0.5; 0.5 0.5]'); hold on;
        plot(Theta_list, Dtot_alt(ia,:),  '--s','LineWidth',1.8, 'MarkerSize',5, ...
            'DisplayName','C alt [0.9 0.1; 0.1 0.9]');
        set(gca,'XScale','log');
        ylim([0 5e-3])
        ylabel('Cumulative deaths at day 365')
        title(sprintf('Deaths vs \\theta (\\alpha = %g)', a));
        set(gca,'fontsize',fs_all)

        xlabel('\theta = \delta_y/\delta_o (log scale)');
        if a==0
            legend('Location','best');
        end

        savepng(fullfile(outdir,sprintf('figC5_Dtot_deaths_%d.png',a)));
        savefig(gcf, fullfile(outdir,sprintf('figC5_Dtot_deaths_%d.fig',a)));
        close;
    end
end

%% ========= Appendix C6: Susceptibility heterogeneity with heterogeneous C (Fractional Analysis) ========
% C_base = [0.5 0.5; 0.5 0.5] & C_het  = [0.9 0.1; 0.1 0.9]
if type.figC6==1

    fprintf('\n=== Fig C6: Fractional Susceptibility heterogeneity with heterogeneous C  ===\n')
    C_alt  = C_het;

    ry = 0.5;
    ro = 1 - ry;
    % ICs consistent with ry, ro
    Iy0 = ry*I0; Io0 = ro*I0;
    Sy0 = ry - Iy0; So0 = ro - Io0;
    y0_2 = [Sy0; 0; Iy0; 0; 0;   So0; 0; Io0; 0; 0];
    Dtot_base = zeros(numel(alpha_list), numel(rho_list));
    Dtot_het  = zeros(numel(alpha_list), numel(rho_list));
    Rtot_base = zeros(numel(alpha_list), numel(rho_list));
    Rtot_het  = zeros(numel(alpha_list), numel(rho_list));
    for ia = 1:numel(alpha_list)
        a = alpha_list(ia);

        for k = 1:numel(rho_list)
            rho = rho_list(k);
            [beta_y, beta_o] = beta_from_rho(rho, ry, beta0);

            % base C
            pB = p2;
            pB.alpha  = a;
            pB.C      = C_base;
            pB.ry     = ry;  pB.ro = ro;
            pB.beta_y = beta_y;
            pB.beta_o = beta_o;
            pB.use_fractional = 1;

            [tB, yB] = ode45(@(t,y) rhs_twogroup(t,y,pB), tspan, y0_2, opts2);
            Rtot_base(ia,k) = yB(end,4) + yB(end,9);
            IyB = yB(:,3); IoB = yB(:,8);
            FyB = pB.delta_y .* IyB ./ tauI;
            FoB = pB.delta_o .* IoB ./ tauI;
            Dtot_base(ia,k) = trapz(tB, FyB + FoB);

            %C_het
            pH = p2;
            pH.alpha  = a;
            pH.C      = C_alt;
            pH.ry     = ry;  pH.ro = ro;
            pH.beta_y = beta_y;
            pH.beta_o = beta_o;
            pH.use_fractional = 1;

            [tH, yH] = ode45(@(t,y) rhs_twogroup(t,y,pH), tspan, y0_2, opts2);
            Rtot_het(ia,k) = yH(end,4) + yH(end,9);


            IyH = yH(:,3); IoH = yH(:,8);
            FyH = pH.delta_y .* IyH ./ tauI;
            FoH = pH.delta_o .* IoH ./ tauI;
            Dtot_het(ia,k) = trapz(tH, FyH + FoH);

        end
    end


    % Combined Plot: Infections: R_tot
    % Goal: Compare R_tot(365) vs theta = beta_y/beta_o under two contact matrices.
    % C_base = [0.5 0.5; 0.5 0.5] & C_het  = [0.9 0.1; 0.1 0.9]
    for ia = 1:numel(alpha_list)
        figure;
        set(gcf, 'Position', [10 10 900 600]); pause(.1)
        a = alpha_list(ia);
        plot(rho_list, Rtot_base(ia,:), '-o', 'LineWidth',1.8, 'MarkerSize',5, ...
            'DisplayName','C base [0.5 0.5; 0.5 0.5]'); hold on;
        plot(rho_list, Rtot_het(ia,:),  '--s','LineWidth',1.8, 'MarkerSize',5, ...
            'DisplayName','C alt [0.9 0.1; 0.1 0.9]');
        set(gca,'XScale','log');
        ylim([0 1])
        ylabel('Cumulatives infections at day 365')
        title(sprintf('Infections vs \\rho (\\alpha = %g)', a));
        set(gca,'fontsize',fs_all)
        xlabel('\rho = \beta_y/\beta_o (log scale)');

        if a ==0
            legend('Location','best');
        end

        savepng(fullfile(outdir,sprintf('figC6_Rtot_cases_%d.png',a)));
        savefig(gcf, fullfile(outdir,sprintf('figC6_Rtot_cases_%d.fig',a)));
        close;
    end

    % Combined Plot: Deaths: D_tot
    % Compare total deaths D_tot = ∫ (F_y + F_o) dt vs theta under two C matrices.
    % C_base  & C_het
    for ia = 1:numel(alpha_list)
        figure;
        set(gcf, 'Position', [10 10 900 600]);pause(.1)
        a = alpha_list(ia);

        plot(rho_list, Dtot_base(ia,:), '-o', 'LineWidth',1.8, 'MarkerSize',5, ...
            'DisplayName','C base [0.5 0.5; 0.5 0.5]'); hold on;
        plot(rho_list, Dtot_het(ia,:),  '--s','LineWidth',1.8, 'MarkerSize',5, ...
            'DisplayName','C alt  [0.9 0.1; 0.1 0.9]');
        set(gca,'XScale','log');
        ylabel('Cumulative deaths at day 365')
        title(sprintf('Deaths vs \\rho (\\alpha = %g)', a));
        ylim([0 5e-3])
        set(gca,'fontsize',fs_all)

        xlabel('\rho = \beta_y/\beta_o (log scale)');

        if a==0
            legend('Location','best');
        end

        savepng(fullfile(outdir,sprintf('figC6_Dtot_deaths_%d.png',a)));
        savefig(gcf, fullfile(outdir,sprintf('figC6_Dtot_deaths_%d.fig',a)));
        close;
    end
end

%% ========= Appendix D: Aggregate Mortality Signal Sweep 
if type.figD1 == 1
    fprintf('\n=== Fig D1 (Appendix D): Vary Mortality (theta) under Aggregate Mortality ===\n')
    
    % Defining the 5 alpha pairs (alpha_y, alpha_o)
    alpha_pairs = [
        200000, 200000;  
        100000, 300000;  
        0,      400000;  
        300000, 100000;  
        400000, 0        
    ];
    num_pairs = size(alpha_pairs, 1);
    num_thetas = numel(Theta_list);
    
    % Pre-allocate storage for Total results: rows = alpha pairs, cols = Theta
    Rtot_365 = zeros(num_pairs, num_thetas);
    Dtot_365 = zeros(num_pairs, num_thetas);
    
    ry = 0.5; 
    ro = 0.5;
    
    % Loop through each of the 5 alpha pairs
    for i_pair = 1:num_pairs
        ay = alpha_pairs(i_pair, 1);
        ao = alpha_pairs(i_pair, 2);
        
        % Sweep across the exact same Theta list (from 10^-4 to 10^0)
        for k = 1:num_thetas
            Theta = Theta_list(k);
            [dy, do] = delta_from_theta(Theta, ry, delta0);
            
            pA = p2;
            pA.alpha_y = ay; 
            pA.alpha_o = ao;
            pA.ry      = ry; 
            pA.ro      = ro; 
            pA.delta_y = dy; 
            pA.delta_o = do;
            pA.use_aggregate = 1; % Flags the RHS to use total (Fty + Fto)
            
            Iy0 = ry*I0; 
            Io0 = ro*I0; 
            Sy0 = ry - Iy0;
            So0 = ro - Io0;
            y0_2 = [Sy0; 0; Iy0; 0; 0;  So0; 0; Io0; 0; 0];
            
            [tT, yT] = ode45(@(t,y) rhs_twogroup(t,y,pA), tspan, y0_2, opts2);
            
            % Compute Total Cumulative Infections at day 365
            Rtot_365(i_pair, k) = yT(end,4) + yT(end,9);
            
            % Integrate total daily deaths to get Total Cumulative Deaths
            Iy_t = yT(:,3); 
            Io_t = yT(:,8);
            Fy_t = dy .* Iy_t ./ tauI; 
            Fo_t = do .* Io_t ./ tauI;
            Dtot_365(i_pair, k) = trapz(tT, Fy_t + Fo_t);
        end
    end
    
    % Markers and line styles for distinguishing the 5 pairs on the plot
    plot_styles = {'-o', '--s', '-.^', ':d', '-x'};
    plot_colors = [
        0      0.4470 0.7410;  % Blue
        0.8500 0.3250 0.0980;  % Orange
        0.9290 0.6940 0.1250;  % Yellow
        0.4940 0.1840 0.5560;  % Purple
        0.4660 0.6740 0.1880   % Green
    ];
    
    labels = {
        '(\alpha_young, \alpha_old) = (200k, 200k)', ...
        '(\alpha_young, \alpha_old) = (100k, 300k)', ...
        '(\alpha_young, \alpha_old) = (0, 400k)', ...
        '(\alpha_young, \alpha_old) = (300k, 100k)', ...
        '(\alpha_young, \alpha_old) = (400k, 0)'
    };

% Creating the Two-Panel Figure (Left: Infections, Right: Deaths)

    f_aggregate = figure('Name', 'Appendix_D_Aggregate_Signal_Robustness');
    set(f_aggregate, 'Position', [10 10 1400 600]); pause(.1)
    
    t_layout = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Minimal and clear legend labels inside the panel
    labels_inside = {
        '\alpha_y = 200,000, \alpha_o = 200,000', ...
        '\alpha_y = 100,000, \alpha_o = 300,000', ...
        '\alpha_y = 0, \alpha_o = 400,000', ...
        '\alpha_y = 300,000, \alpha_o = 100,000', ...
        '\alpha_y = 400,000, \alpha_o = 0'
    };
    
    % --- PANEL 1: Total Infections vs Theta
    nexttile;
    hold on;
    for i_pair = 1:num_pairs
        plot(Theta_list, Rtot_365(i_pair, :), plot_styles{i_pair}, ...
            'LineWidth', lw, 'Color', plot_colors(i_pair, :), ...
            'MarkerSize', 6);
    end
    set(gca, 'XScale', 'log');
    xlabel('\theta = \delta_y / \delta_o (log scale)');
    ylabel('Cumulative infections at day 365', 'FontSize', fs_all);
    ylim([0 1]);
    xticks([1e-4 1e-3 1e-2 1e-1 1e-0]);
    box on; 
    set(gca, 'FontSize', fs_all);
    legend(labels_inside, 'Location', 'northeast', 'FontSize', fs_legend,'box','on');
    
    % --- PANEL 2: Total Deaths vs Theta
    nexttile;
    hold on;
    for i_pair = 1:num_pairs
        plot(Theta_list, Dtot_365(i_pair, :), plot_styles{i_pair}, ...
            'LineWidth', lw, 'Color', plot_colors(i_pair, :), ...
            'MarkerSize', 6);
    end
    set(gca, 'XScale', 'log');
    xlabel('\theta = \delta_y / \delta_o (log scale)');
    ylabel('Cumulative deaths at day 365', 'FontSize', fs_all);
    ylim([0 5e-3]);
    xticks([1e-4 1e-3 1e-2 1e-1 1e-0]);
    box on; 
    set(gca, 'FontSize', fs_all);
    
    
    % Placing the legend inside the second panel (Deaths) on the Top-Right
    
    
    % Save outputs
    savepng(fullfile(outdir, 'figD1_aggregate_signal_sweep.png'));
    savefig(f_aggregate, fullfile(outdir, 'figD1_aggregate_signal_sweep.fig'));
    close(f_aggregate);
end
fprintf('\nSimulation finished. Publication figures saved in %s\n', outdir);
end

%% ================== Determine delta from theta ======================
function [dy, do] = delta_from_theta(theta, ry, delta0)
den = 1 + ry*theta - ry;
assert(den > 1e-12, 'Denominator too small: 1 + ry*theta - ry');
dy  = (theta/den) * delta0;
do  = (1/den)    * delta0;
end


%% ================== Determine beta from rho =========================
function [beta_y, beta_o] = beta_from_rho(rho, ry, bbar)
den = 1 +ry*rho-ry;
beta_y = bbar*rho/den;
beta_o = bbar/den;
end


%% =================== RHS FUNCTIONS ==================================

function dy = rhs_twogroup(t, y, p)
Sy=y(1); 
Ey=y(2);
Iy=y(3);
Ry=y(4);
Fty=y(5);
So=y(6);
Eo=y(7);
Io=y(8);
Ro=y(9);
Fto=y(10);
% per-capita infectious pressure (homogeneous mixing)
py = (p.C(1,1) * (Iy/p.ry) + p.C(1,2) * (Io/p.ro));
po = (p.C(2,1) * (Iy/p.ry) + p.C(2,2) * (Io/p.ro));

if isfield(p, 'use_aggregate') && p.use_aggregate == 1
    % Conceptually calculated as:
    % 1. Aggregate actual death rates: f = f_y + f_o
    % 2. Apply perception delay: d(f_tilde)/dt = (f - f_tilde)/tau_f
    % Due to identical tau_f, this is mathematically equivalent to pooling the group-specific lags:
    f_tilde_total = (Fty + Fto) / (p.ry + p.ro);
    qy = exp(- p.alpha_y * p.k_y * f_tilde_total);
    qo = exp(- p.alpha_o * p.k_o * f_tilde_total);
elseif isfield(p, 'use_fractional') && p.use_fractional == 1
    gamma_pow = 2; 
    qy = 1 / ((1 + p.alpha * p.k_y * (Fty / p.ry))^gamma_pow);
    qo = 1 / ((1 + p.alpha * p.k_o * (Fto / p.ro))^gamma_pow);
else
    qy = exp(- p.alpha * p.k_y * (Fty / p.ry));
    qo = exp(- p.alpha * p.k_o * (Fto / p.ro));
end

if isfield(p,'beta_y'), by = p.beta_y; else, by = p.beta0; end
if isfield(p,'beta_o'), bo = p.beta_o; else, bo = p.beta0; end

Fy = p.delta_y * Iy / p.tauI;
Fo = p.delta_o * Io / p.tauI;

dSy = - by * qy * Sy * py;
dEy =   by * qy * Sy * py - Ey/p.tauE;
dIy =   Ey/p.tauE - Iy/p.tauI;
dRy =   Iy/p.tauI;
dFty= (Fy - Fty)/p.tauF;
dSo = - bo * qo * So * po;
dEo =   bo * qo * So * po - Eo/p.tauE;
dIo =   Eo/p.tauE - Io/p.tauI;
dRo =   Io/p.tauI;
dFto= (Fo - Fto)/p.tauF;

dy = [dSy; dEy; dIy; dRy; dFty;  dSo; dEo; dIo; dRo; dFto];
end
% qy = exp(- p.alpha * p.k_y * (Fty / p.ry));
% qo = exp(- p.alpha * p.k_o * (Fto / p.ro));
% if isfield(p,'beta_y'), by = p.beta_y; else, by = p.beta0; end
% if isfield(p,'beta_o'), bo = p.beta_o; else, bo = p.beta0; end
% Fy = p.delta_y * Iy / p.tauI;
% Fo = p.delta_o * Io / p.tauI;
% dSy = - by * qy * Sy * py;
% dEy =   by * qy * Sy * py - Ey/p.tauE;
% dIy =   Ey/p.tauE - Iy/p.tauI;
% dRy =   Iy/p.tauI;
% dFty= (Fy - Fty)/p.tauF;
% dSo = - bo * qo * So * po;
% dEo =   bo * qo * So * po - Eo/p.tauE;
% dIo =   Eo/p.tauE - Io/p.tauI;
% dRo =   Io/p.tauI;
% dFto= (Fo - Fto)/p.tauF;
% dy = [dSy; dEy; dIy; dRy; dFty;  dSo; dEo; dIo; dRo; dFto];
% end


%% ===================  Save helper ==================================

function savepng(fname)
[d,~,ext]=fileparts(fname);
if ~exist(d,'dir'), mkdir(d); end
if isempty(ext), fname=[fname '.png']; end
try
    exportgraphics(gcf, fname, 'Resolution', 300, 'Padding', 'loose'); 
catch
    print(gcf, fname, '-dpng', '-r300'); 
end
end


