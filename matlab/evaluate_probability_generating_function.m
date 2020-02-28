function [x,y] = evaluate_probability_generating_function(G, evaluate_num_timeslots, use_taylor_approx)    
    %% Returns probability mass function (PMF) y that corresponds to probability generating function (PGF) G for x=0:evaluate_num_timeslots    
    x = 0:evaluate_num_timeslots;    
    %% Use retrieval of PMF through k-times differentiation.
    if use_taylor_approx == 0
        syms z;
        current_derivative = G;
        y = zeros(evaluate_num_timeslots+1, 1);           
        for timeslot = 0:evaluate_num_timeslots
            % Derive for z.
            current_derivative = diff(current_derivative, z);        
            % Transform to PMF through 1/t! * derivative(z=0).                        
            k = 1 + timeslot;           
            pmf = 1/factorial(k) * current_derivative(0);        
            y(1+timeslot) = pmf;         
        end      
        y = [0; y(1:end-1)];    
    %% Use Taylor series approximation (performs much faster than above solution).
    else
        taylor_series = sym2poly(taylor(G, 'Order', evaluate_num_timeslots));        
        taylor_series = fliplr(taylor_series);           
        taylor_series = taylor_series(2:end); % first derivative is in 2nd summand.                       
        y = [0 taylor_series zeros(1, evaluate_num_timeslots - length(taylor_series))];        
    end
end

