data = csvread('mac_delay_simulte.csv');
%% Matrix structure: runs x num_packets.
num_runs = size(data, 2) / 2; % Columns are: time_run1, delay_run2, time_run2, delay_run2...
num_packets = size(data, 1); % Rows.
timeMat = zeros(num_runs, num_packets);
delayMat = zeros(num_runs, num_packets);
max_num_transmission_attempts = 5; % Initial transmission + up to four HARQ retransmissions.

%% Fill matrices.
for packet=1:num_packets % Packets go down the rows of the 'data' matrix.    
    run = 1;
    for j=1:size(data,2) % Go along the columns.
        if mod(j,2) == 1 % Odd numbers are time values.
            timeMat(run, packet) = data(packet, j);
        else % Even numbers are delay values.
            delayMat(run, packet) = data(packet, j);
            run = run + 1;
        end        
    end
end

%% Derive number of transmission attempts per run.
txAttemptMat = zeros(num_runs, max_num_transmission_attempts);
for run=1:num_runs    
    nonzero_delayMat = delayMat(delayMat(run,1:end)>0);
    binVec = discretize(nonzero_delayMat, max_num_transmission_attempts);
    for txAttempt = 1:max_num_transmission_attempts
        txAttemptMat(run, txAttempt) = sum(binVec==txAttempt);
    end
end

% %% Plot transmission attempts.
% figure;
% bar(0:max_num_transmission_attempts-1, sum(txAttemptMat))
% xlabel('\#{}Retransmissions', 'Interpreter', 'latex');
% ylabel('Count', 'Interpreter', 'latex');

%% Derive occurrence probabilities.
pMat = zeros(num_runs, max_num_transmission_attempts);
for run=1:num_runs
    sum_attempts = sum(txAttemptMat(run, 1:end));
    for attempt=1:max_num_transmission_attempts
        pMat(run, attempt) = txAttemptMat(run, attempt) / sum_attempts;
    end
end

ccdfMat = zeros(num_runs, max_num_transmission_attempts); % CCDF = P(X>x)
for run=1:num_runs    
    ccdfMat(run, :) = 1-cumsum(pMat(run, 1:end));
end

%% Calculate confidence intervals.
pdfMeanVec = mean(pMat);
ccdfMeanVec = 1 - cumsum(pdfMeanVec);
N = size(ccdfMeanVec, 2);
standard_error_of_mean = std(pMat)/sqrt(N);
t_dist_interval = tinv([0.05 0.95], N-1);
confidence_intervals = bsxfun(@times, standard_error_of_mean, t_dist_interval(:));
ccdfMeanVec = [1 ccdfMeanVec];
x_sim = 0:max_num_transmission_attempts;

%% Plot.
% figure 
% stairs(x_sim, ccdfMeanVec);
% hold on;
% errorbar(x_sim, ccdfMeanVec, [0 confidence_intervals(2,:)], [0 confidence_intervals(2,:)]);