function [T_S, T_F] = msfg_max_num_rtx_erroneous_feedback(err_good, err_bad, max_num_rtx, p12, p21, rtt, e_NA0, e_NA1)
    %% Set up variables.
    N = max_num_rtx;
    L = [err_good 0; 0 err_bad];
    pi1 = p21 / (p12 + p21); % Steady state distribution: good state.
    pi2 = p12 / (p21 + p12); % Steady state distribution: bad state.
    steady_state = [pi1 pi2]; % Steady state distribution.           
    P = [1-p12 p12; p21 1-p21];    
    I = eye(size(P));
    E_NA = [e_NA0 0; 0 e_NA1];
    E_NN = I - E_NA;    
    
    %% Construct transfer function to success state S.
    syms z;
    v = rtt;
    sum = 0;
    for i=0:N
        sum = sum + (z^v*L*E_NN*P)^i*z^v*(I-L);
    end        
    T_S = symfun(sum, z);
    % Transform to scalar PGF.
    T_S = (steady_state * T_S * ones(size(P, 1), 1)) / (steady_state * ones(size(P, 1), 1));    
    
    %% Construct transfer function to failure state F.
    sum = 0;
    for i=0:N-1
        sum = sum + (E_NN*P*z^v*L)^i;
    end
    T_F = symfun(z^v*L*sum*E_NA + (z^v*L*E_NN*P)^N*z^v*L, z);
    T_F = (steady_state * T_F * ones(size(P, 1), 1)) / (steady_state * ones(size(P, 1), 1));    
end

