function [T_S, T_F] = msfg_rlc(T_S_harq, T_F_harq, p12, p21, max_num_rtx, rtt)    
    N = max_num_rtx;    
    v = rtt;
        
    %% A status report is sent by the receiver when an out-of-order reception occurs. 
    % This informs the sender that a PDU has been lost, which starts a timer that expires when reception is possible at the latest.    
    timeslots_until_loss_detected = N+1;        
        
    %% Construct transfer function to success state S.
    syms z;    
    T_S = symfun(0, z);
    for i=0:N
        T_S = T_S + (T_F_harq * z^(v*timeslots_until_loss_detected))^i * T_S_harq;
    end    
    
    T_F = (T_F_harq * z^(v*timeslots_until_loss_detected))^N * T_F_harq;
    
    % Transform to scalar PGF.
    P = [1-p12 p12; p21 1-p21];
    pi1 = p21 / (p12 + p21); % Steady state distribution: good state.
    pi2 = p12 / (p21 + p12); % Steady state distribution: bad state.
    steady_state = [pi1 pi2]; % Steady state distribution.           
    T_S = (steady_state * T_S * ones(size(P, 1), 1)) / (steady_state * ones(size(P, 1), 1));        
    T_F = (steady_state * T_F * ones(size(P, 1), 1)) / (steady_state * ones(size(P, 1), 1));        
end

