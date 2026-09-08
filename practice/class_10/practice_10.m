function practice_10()
% Calculations for the Class 10 Practice Assignment from REB, The Course
    % global constants available to all functions
    % given
    V = 4; % L
    k01 = 2.59e9; % min-1
    E1 = 16500; % cal/mol
    dH1 = -22200; % cal/mol
    Cp = 440; % cal/L/K
    A = 0.6; % ft2
    U = 5.65e3/60; % cal/ft2/min/K
    Cp_ex = 1000; % cal/kg/K
    CA_in = 2; % mol/L
    T_in = 60 + 273.15; % K
    Tex_in = 20 + 273.15; % K
    m_ex = 0.2; % kg/min
    tau_values = [50, 60]; % min
    % known
    R = 1.987; % cal/mol/K

    % global variable for the current value of tau
    g_tau = nan;

    % CSTR model function
    function [nA, nZ, T, Tex] = cstr_model_variables(tau, guess)
        % make Vdot_in available to the residuals function
        g_tau = tau;

        % solve the design equations
        [soln, success, message] = solve_ates(@cstr_residuals, guess);

        % check for solver issues
        if ~success
            disp(' ')
            disp(['    CSTR model function issue: ', message])
            disp(' ')
        end
        
        % extract and return the results
        nA = soln(1);
        nZ = soln(2);
        T = soln(3);
        Tex = soln(4);
    end

    % CSTR residuals function
    function epsilon = cstr_residuals(guess)
        % extract the individual guesses
        nA = guess(1);
        nZ = guess(2);
        T = guess(3);
        Tex = guess(4);

        % calculate the additional unknowns
        Vdot_in = V/g_tau;
        Vdot = Vdot_in;
        nA_in = CA_in*Vdot_in;
        nZ_in = 0;
        k1 = k01 * exp(-E1 / (R * T));
        CA = nA / Vdot;
        r1 = k1 * CA;
        Qdot = U * A * (Tex - T);
        
        % evaluate and return the CSTR design equation residuals
        epsilon_1 = nA_in - nA - r1 * V;
        epsilon_2 = nZ_in - nZ + r1 * V;
        epsilon_3 = Vdot_in*Cp*(T - T_in) + V*r1*dH1 - Qdot;
        epsilon_4 = -m_ex*Cp_ex*(Tex - Tex_in) - Qdot;

        % return the residuals as an array
        epsilon = [epsilon_1; epsilon_2; epsilon_3; epsilon_4];
    end

    % deliverables function
    function deliverables()
        % define an initial guess for the cstr model variables
        Vdot_in = V/tau_values(1);
        nA_in = CA_in*Vdot_in;
        guess = [nA_in, 0, T_in + 5, Tex_in + 5];

        % allocate storage for the conversions and temperatures corresponding to each space time
        fA = [nan, nan];
        T = [nan, nan];
        Tex = [nan, nan];

        % loop through the space times
        for i = 1:2
            % solve the cstr model for the current space time
            [nA, nZ, T(i), Tex(i)] = cstr_model_variables(tau_values(i), guess);

            % calculate the conversion
            Vdot_in = V/tau_values(i);
            nA_in = CA_in*Vdot_in;
            fA(i) = 100*(nA_in - nA)/nA_in;

            % use the results as the guess for the next iteration
            guess = [nA, nZ, T(i), Tex(i)];
        end

        % tabulate, display, and save the results
        results_table = table(tau_values', fA', T' - 273.155, Tex' - 273.155, 'VariableNames', {'Space Time (min)', 'Conversion (%)', 'T (°C)', 'Tex     (°C)'});
        disp(' ')
        disp(results_table);
        disp(' ')
        writetable(results_table, 'prac_10_results.csv');
    end

    % perform the calculations
    deliverables()
end