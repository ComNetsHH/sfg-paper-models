function [T_S, T_F] = msfg_max_num_rtx_erroneous_feedback_with_soft_combining(err_good, err_bad, max_num_rtx, p12, p21, rtt, e_NA0, e_NA1, scalarize)
    %% Set up variables.
    N = max_num_rtx;
    L = [err_good 0; 0 err_bad];
    pi_A = p21 / (p12 + p21); % Steady state distribution: good state.
    pi_B = p12 / (p21 + p12); % Steady state distribution: bad state.
    steady_state = [pi_A pi_B]; % Steady state distribution.           
    P = [1-p12 p12; p21 1-p21];    
    I = eye(size(P));
    E_NA = [e_NA0 0; 0 e_NA1];
    E_NN = I - E_NA;
    syms z;
    v = rtt;        
    
%     expected_nack_ack_misunderstanding = steady_state * E_NA * ones(size(steady_state))';
     %disp('Gilbert-Elliot steady state distribution, good state=' + string(pi_A) + ' bad state=' + string(pi_B));    
%     disp('Expected NACK->ACK misinterpretation rate: ' + string(expected_nack_ack_misunderstanding));
    
    %% Construct transfer function to success state S.    
    sum = 0;   
    for i=0:N        
        product = 1;        
        for j=0:i-1
            product = product * a(j, z, v, pi_A, pi_B, L)*d(P, E_NN, v);
        end
        sum = sum + product*b(i, z, v, pi_A, pi_B, I, L);
    end        
    T_S = symfun(sum, z);     
    % Transform to scalar PGF.
    if scalarize > 0
        T_S = (steady_state * T_S * ones(size(P, 1), 1)) / (steady_state * ones(size(P, 1), 1));
    end
    
    %% Construct transfer function to failure state F.
    sum = 0;
    for i=0:N-1
        product = 1;
        for j=1:i
            product = product * d(P,E_NN,v)*a(j,z,v,pi_A,pi_B,L);
        end
        sum = sum + product*c(E_NA);
    end
    product = 1;
    for j=0:N-1
        product = product*a(j,z,v,pi_A,pi_B,L)*d(P,E_NN,v);
    end
    product = product*a(N,z,v,pi_A,pi_B,L);
    T_F = symfun(a(0,z,v,pi_A,pi_B,L)*sum + product, z);
    % Transform to scalar PGF.
    if scalarize > 0
        T_F = (steady_state * T_F * ones(size(P, 1), 1)) / (steady_state * ones(size(P, 1), 1));
    end
end

function a_i = a(i, z, v, pi_A, pi_B, L)
    sum = 0;
    for j=0:i
        k=i-j;
        sum = sum + nchoosek(i,j)*pi_A^j*pi_B^k*L^(max(1,j+1));
    end
    a_i = symfun(z^v * sum, z);
end

function b_i = b(i, z, v, pi_A, pi_B, I, L)
    sum = 0;
    for j=0:i
        k=i-j;
        sum = sum + nchoosek(i,j)*pi_A^j*pi_B^k*(I - L^(max(1,j+1)));
    end
    b_i = symfun(z^v * sum, z);
end

function c = c(E_NA)
    c = E_NA;
end

function d = d(P, E_NN, v)
    d = E_NN * P;
end