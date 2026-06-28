% CEST analysis using double-step multipool Lorentzian fitting (DMPLF)
% Jianpan Huang - jianpanhuang@outlook.com
clear all, close all, clc
addpath(genpath(pwd));

%% Load data
load('CEST_5pools_VariedConc_Noise0p002.mat'); % The mat file includes Z-spectrum data (zspec), frequency offsets (offs) and fractions of four exchange pools
if_load_result = 0; % 0 will run the fitting, 1 will load the fitted results

%% Seperate M0 and Z-spectrum data
zspec_m0 = zspec(1,:);
zspec(1,:) = [];
offs(1) = [];
zspec_norm = zspec./repmat(zspec_m0, [length(offs),1]);
zspec_num = size(zspec,2);

%% Z-spectrum fitting
if if_load_result == 0
    pn = 5; % pool number
    %             1. Water              2. Amide               3. rNOE                4. MT                  5. Guan
    %      Zi     A1    G1    dw1       A2     G2    dw2       A3     G3    dw3       A4     G4    dw4       A5     G5    dw5
    iv1 = [ 1     0.9   2.3   0         0.0   2.0   3.5        0.01   4.0   -3.5      0.1    50    -2.5      0.0    2.0   2.0];
    lb1 = [ 1     0.1   0.5   -1.5      0.0    0.5   3.5       0.0    1.0   -3.5      0.0025 10    -2.5      0.0    0.5   2.0];
    ub1 = [ 1     1     8.0   +1.5      0.0    12.5  3.5       0.4    35.0  -3.5      0.50   130   -2.5      0.0    15    2.0];

    fit_para = zeros(zspec_num, length(iv1));
    h = waitbar(0, 'Doing Lorentzian fitting voxel by voxel, please wait >>>>>>'); 
    set(h, 'Units', 'normalized', 'Position', [0.4, 0.2, 0.25, 0.08])
    count = 0;
    tic
    for nn = 1:zspec_num
        count = count+1;
        zspec = zspec_norm(:,nn);

        % step1: fit DS, rNOE, MT
        fit_center = [-0.6,0.6]; % in ppm
        fit_edge = [0.6,9]; % in ppm
        fit_ind_center = find(offs>fit_center(1)&offs<fit_center(2));
        fit_ind_neg = find(offs<fit_edge(1));
        fit_ind_pos = find(offs>fit_edge(2));
        fit_ind = [fit_ind_neg;fit_ind_center;fit_ind_pos];
        % options = optimset('MaxFunEvals', 1e3,'TolFun', 1e-6,'TolX', 1e-6, 'Display',  'off' );
        par1 = lsqcurvefit(@lorentzian5pool,iv1,offs(fit_ind),zspec(fit_ind),lb1,ub1);

        % step2: fit amide, guan
        %      Zi     A1    G1    dw1       A2     G2    dw2       A3     G3    dw3       A4     G4    dw4       A5     G5    dw5
        iv = [ 1      0.9   2.3   0         0.01   2.0   3.5       0.01   4.0   -3.5      0.1    50    -2.5      0.01   2.0   2.0];
        lb = [ 1      0.1   0.5   -1.5      0.0    0.5   3.5       0.0    1.0   -3.5      0.0025 10    -2.5      0.0    0.5   2.0];
        ub = [ 1      1     8.0   +1.5      0.25    12.5  3.5      0.4   35.0  -3.5      0.50   130   -2.5      0.25    15    2.0];
        iv([1:4,8:13]) = par1([1:4,8:13]); % fix parameters for DS, rNOE and MT
        lb([1:4,8:13]) = par1([1:4,8:13]); % fix parameters for DS, rNOE and MT
        ub([1:4,8:13]) = par1([1:4,8:13]); % fix parameters for DS, rNOE and MT
        [par, resnorm, residual, exitflag] = lsqcurvefit(@lorentzian5pool,iv,offs,zspec,lb,ub);

        cr(nn) = par(14);
        apt(nn) = par(5);
        rnoe(nn) = par(8);
        mt(nn) = par(11);
        fit_para(nn,:) = par;
        waitbar(count/zspec_num, h);
        % fitted curves
        cur = zspec_analysis(par, offs, pn);
        if mod(count, 200) == 0
            zspecfit = cur(:,1);
            water_cur = cur(:,2);
            apt_cur = cur(:,3);
            rnoe_cur = cur(:,5);
            mt_cur = cur(:,7);
            cr_cur = cur(:,9);
            figure(100);
            plot(offs, zspec, 'bo', offs, zspecfit, 'b-',...
            offs, water_cur, 'r-.', offs, cr_cur, 'k-.', offs, apt_cur, 'g-.',...
            offs, rnoe_cur, 'c-.', offs, mt_cur, 'm-.', 'LineWidth',1.5);
            axis([min(offs(:)),max(offs(:)),0,1.01]); 
            xlabel('Offset (ppm)'); ylabel('Z');
            title('Z-spectrum'); 
            legend('Z', 'Z_f_i_t', 'Water', 'CrCET', 'APT', 'rNOE', 'MT','Location', 'east');
            set(gca, 'Xdir', 'reverse', 'FontWeight', 'bold', 'FontSize', 14, 'LineWidth',3)
            set(gcf,'color','w');
        end   
    end
    toc
    close gcf;
    delete(h);
    save('Result_DMPLF.mat', 'apt','rnoe','mt','cr','offs','zspec','zspec_norm','zspec_m0','fit_para');
else
    load('Result_DMPLF.mat');
end

%% Plot correlation
params = {'APT','CrCEST','rNOE','MT'};
true_rows = [1,2,3,4];
meas_cells = {apt, cr, rnoe, mt};

screen_size = get(0, 'ScreenSize'); % [left, bottom, width, height]
fig_width = screen_size(3) * 1;
fig_height = screen_size(4) * 0.4;
fig_left = (screen_size(3) - fig_width) / 2; 
fig_bottom = (screen_size(4) - fig_height) / 2;
figure('Position', [fig_left, fig_bottom, fig_width, fig_height]);

for i = 1:4
    conc_vals = frac_apt_cr_rnoe_mt(true_rows(i), :);
    meas_vals = meas_cells{i}(:);
    conc_vals = conc_vals(:)*110000; % change to mM by normalizing to water proton concentration (110 M)

    % plot
    subplot(1, 4, i); 
    scatter(conc_vals, meas_vals, 30, 'filled', 'MarkerFaceColor', [0.2 0.4 0.8]);
    hold on;

    % fit
    coeffs = polyfit(conc_vals, meas_vals, 1);
    a = coeffs(1);
    b = coeffs(2);
    x_fit = [min(conc_vals), max(conc_vals)];
    y_fit = a * x_fit + b;
    plot(x_fit, y_fit, 'r-', 'LineWidth', 3); 

    % calculate Pearson r and p values
    [R, P] = corrcoef(conc_vals, meas_vals);
    r = R(1,2);
    p = P(1,2);
    R2 = r^2;  
    
    xlabel('Concentration (mM)', 'FontSize', 14);
    ylabel('CEST', 'FontSize', 14);
    title(['\bf' params{i}], 'FontSize', 16);
    text(0.05, 0.90, sprintf('r = %.3f', r), ...
        'Units', 'normalized', 'FontSize', 18, 'FontWeight', 'bold');

    grid on;
    hold off;

    ax = gca;
    ax.LineWidth = 3;          
    ax.FontSize = 14;          
    ax.FontWeight = 'bold';    
end