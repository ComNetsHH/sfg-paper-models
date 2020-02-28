err_good = 0.1;
err_bad = 1.0;
p12 = 0.25;
p21 = 0.5;
rtt = 4;
max_num_rtx = 4;
evaluate_num_timeslots = rtt * (max_num_rtx+1) * 10;
e_NA0 = 0.001;
e_NA1 = 0.001;
line_width = 1.0;
font_size = 14;
use_taylor_approx = 0;

%% Stop-and-Wait ARQ with at most N retransmissions.
[T_S_ideal, T_F_ideal] = msfg_max_num_rtx(err_good, err_bad, max_num_rtx, p12, p21, rtt);
% Construct full PGF and check for mathematical soundness (must sum to 1).
G_ideal = T_S_ideal + T_F_ideal;
[x_ideal,y_ideal] = evaluate_probability_generating_function(G_ideal, evaluate_num_timeslots, use_taylor_approx);
disp("Sum over all N timeslots in ideal feedback channel SFG's probability generating function is: " + string(sum(y_ideal)));
% Evaluate transfer function to success state only.
[x_ideal,y_ideal] = evaluate_probability_generating_function(T_S_ideal, evaluate_num_timeslots, use_taylor_approx);

%% Stop-and-Wait ARQ with at most N retransmissions and an erroneous feedback channel.
[T_S_erroneous, T_F_erroneous] = msfg_max_num_rtx_erroneous_feedback(err_good, err_bad, max_num_rtx, p12, p21, rtt, e_NA0, e_NA1);
% Construct full PGF and check for mathematical soundness (must sum to 1).
G_erroneous = T_S_erroneous + T_F_erroneous;
[x_erroneous,y_erroneous] = evaluate_probability_generating_function(G_erroneous, evaluate_num_timeslots, use_taylor_approx);
disp("Sum over all N timeslots in erroneous feedback channel SFG's probability generating function is: " + string(sum(y_erroneous)));
[x_erroneous,y_erroneous] = evaluate_probability_generating_function(T_S_erroneous, evaluate_num_timeslots, use_taylor_approx);

%% HARQ with an erroneous feedback channel.
scalarize = 1;
use_taylor_approx = 1;
[T_S_harq, T_F_harq] = msfg_max_num_rtx_erroneous_feedback_with_soft_combining(err_good, err_bad, max_num_rtx, p12, p21, rtt, e_NA0, e_NA1, scalarize);
% Construct full PGF and check for mathematical soundness (must sum to 1).
G_harq = T_S_harq + T_F_harq;
[x_harq_erroneous,y_harq_erroneous] = evaluate_probability_generating_function(G_harq, evaluate_num_timeslots, use_taylor_approx);
disp("Sum over all N timeslots in erroneous feedback channel HARQ SFG's probability generating function is: " + string(sum(y_harq_erroneous)));
% Evaluate transfer function to success state only.
[x_harq_erroneous,y_harq_erroneous] = evaluate_probability_generating_function(T_S_harq, evaluate_num_timeslots, use_taylor_approx);

%% Evaluate simuLTE.
evaluate_simulte

%% RLC-layer ARQ.
scalarize = 0; % Don't scalarize so that we can pass the matrix-form transfer functions into the RLC SFG.
use_taylor_approx = 1;
max_num_rtx_rlc = 3;
[T_S_harq_matrix, T_F_harq_matrix] = msfg_max_num_rtx_erroneous_feedback_with_soft_combining(err_good, err_bad, max_num_rtx, p12, p21, rtt, e_NA0, e_NA1, scalarize);
[T_S_rlc, T_F_rlc] = msfg_rlc(T_S_harq_matrix, T_F_harq_matrix, p12, p21, max_num_rtx_rlc, rtt);
% Construct full PGF and check for mathematical soundness (must sum to 1).
G_rlc = T_S_rlc + T_F_rlc;
[x_rlc,y_rlc] = evaluate_probability_generating_function(G_rlc, evaluate_num_timeslots, use_taylor_approx);
disp("Sum over all N timeslots in RLC Layer SFG's probability generating function is: " + string(sum(y_rlc)));
% Evaluate transfer function to success state only.
[x_rlc, y_rlc] = evaluate_probability_generating_function(T_S_rlc, evaluate_num_timeslots, use_taylor_approx);

%% Plot.
fig = figure;
hold on;
grid on;
set(gca, 'YScale', 'log');
xlabel('[ms] until successful reception', 'Interpreter', 'latex');
ylabel('$P(X>t)$', 'Interpreter','latex');

stairs(x_ideal, 1 - cumsum(y_ideal), '--', 'LineWidth', line_width);
stairs(x_erroneous, 1 - cumsum(y_erroneous), 'LineStyle', '-.', 'LineWidth', line_width);
stairs(x_harq_erroneous, 1 - cumsum(y_harq_erroneous), 'LineWidth', line_width);
stairs(x_sim*rtt, ccdfMeanVec, 'Color', 'black', 'LineWidth', line_width);
errorbar(x_sim*rtt, ccdfMeanVec, [0 confidence_intervals(2,:)], 'LineStyle', 'none', 'Color', 'black', 'LineWidth', 1.2*line_width, 'HandleVisibility','off');

legend('ideal feedback channel', 'erroneous feedback channel', 'erroneous feedback channel + HARQ', 'simulation', 'Location', 'southwest');
set(gca, 'FontSize', font_size);
ylim([10^(-2) 1]);
xlim([0 rtt * (max_num_rtx+1)]);

set(fig, 'Units', 'Inches');
pos = get(fig, 'Position');
set(fig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
print(fig, 'imgs/sfgs-comparison','-dpdf','-r0');
savefig(fig, 'imgs/sfgs-comparison');
hold off;

%% Print error between simulation and analytic approach.
errVec = zeros(size(ccdfMeanVec));
errVecRelative = zeros(size(ccdfMeanVec));
ccdfMeanVec_sfg = 1 - cumsum(y_harq_erroneous);
for i=x_sim(1:end-1)*rtt+1        
    errVec((i-1)/rtt+1) = abs(ccdfMeanVec_sfg(i) - ccdfMeanVec((i-1)/rtt+1));        
    disp(string(ccdfMeanVec_sfg(i)) + ' / ' + string(ccdfMeanVec((i-1)/rtt+1)) + ' = ' + string(ccdfMeanVec_sfg(i) / ccdfMeanVec((i-1)/rtt+1)))
    errVecRelative((i-1)/rtt+1) = ccdfMeanVec_sfg(i) / ccdfMeanVec((i-1)/rtt+1);
end
disp('Maximum absolute error between simulation and SFG is: ' + string(max(errVec)))
disp('Mean absolute error between simulation and SFG is: ' + string(mean(errVec)))
disp('Mean relative error (simulation/SFG) is: ' + string(mean(errVecRelative)))

%% Plot with RLC.
fig_rlc = figure;
hold on;
grid on;
set(gca, 'YScale', 'log');
xlabel('[ms] until successful reception', 'Interpreter', 'latex');
ylabel('$P(X>t)$', 'Interpreter','latex');

stairs(x_ideal, 1 - cumsum(y_ideal), '--', 'LineWidth', line_width);
stairs(x_erroneous, 1 - cumsum(y_erroneous), 'LineStyle', '-.', 'LineWidth', line_width);
stairs(x_harq_erroneous, 1 - cumsum(y_harq_erroneous), 'LineWidth', line_width);
stairs(x_sim*rtt, ccdfMeanVec, 'Color', 'black', 'LineWidth', line_width);
errorbar(x_sim*rtt, ccdfMeanVec, [0 confidence_intervals(2,:)], 'LineStyle', 'none', 'Color', 'black', 'LineWidth', 1.2*line_width, 'HandleVisibility','off');
stairs(x_rlc, 1 - cumsum(y_rlc), 'LineWidth', line_width);

legend('ideal feedback channel', 'erroneous feedback channel', 'erroneous feedback channel + HARQ', 'simulation', 'RLC Layer ARQ', 'Location', 'none', 'Position', [0.15 0.365 1 1]);
set(gca, 'FontSize', font_size);
ylim([10^(-6) 1]);

set(fig_rlc, 'Units', 'Inches');
pos = get(fig_rlc, 'Position');
set(fig_rlc, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
print(fig_rlc, 'imgs/sfgs-comparison_with-rlc', '-dpdf', '-r0');
savefig(fig_rlc, 'imgs/sfgs-comparison_with-rlc');
hold off;