function practice_11()
% Calculations for Practice Assignment 11 of REB, The Course
    % global constants available to all functions 
    % given
    V = 500.0; % ml
    tau = 0.4; % min
    CA_in = 5.0E-3; % mol/ml
    T_in = 60 + 273.15; % K
    Cp = 1.0; % cal/ml/K
    dH1 = -30.0E3; % cal/mol
    k01 = 4.75E13; % /min
    E1 = 25000.; % cal/mol
    % known
    R = 1.987; % cal/mol/K
    % calculated
    Vdot_in = V/tau;
    Vdot = Vdot_in;
    nA_in = CA_in*Vdot_in;
    nZ_in = 0;

    % CSTR model function
    function [nA, nZ, T] = cstr_model_variables(init_guess)
        % solve the design equations
        [soln, success, message] = solve_ates(@cstr_residuals, init_guess);

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
    end

    % CSTR residuals function
    function epsilon = cstr_residuals(guess)
        % extract the individual guesses
        nA = guess(1);
        nZ = guess(2);
        T = guess(3);

        % calculate the additional unknowns
        k1 = k01*exp(-E1/(R*T));
        CA = nA/Vdot;
        r1 = k1*CA;

        % evaluate the residuals
        epsilon_1 = nA_in - nA - V*r1;
        epsilon_2 = nZ_in - nZ + V*r1;
        epsilon_3 = Vdot_in*Cp*(T - T_in) + V*r1*dH1;

        % return the residuals as an array
        epsilon = [epsilon_1; epsilon_2; epsilon_3];
    end

    % deliverables function
    function deliverables()
        % allocate storage for the requested results
        fA = ones(3,1)*nan;
        T = ones(3,1)*nan;

        % low temperature steady state
        % define an initial guess
        initial_guess = [nA_in, nZ_in, T_in + 5.0];
        % solve the CSTR design equations
        [nA, ~, T(1)] = cstr_model_variables(initial_guess);
        % calculate the conversion
        fA(1) = 100*(nA_in - nA)/nA_in;

        % repeat for the high temperature steady state
        % define an initial guess
        %initial_guess(3) = T_in + 100;       failed to converge
        initial_guess(3) = T_in + 200;
        % solve the CSTR design equations
        [nA, ~, T(3)] = cstr_model_variables(initial_guess);
        % calculate the conversion
        fA(3) = 100*(nA_in - nA)/nA_in;

        % repeat for the middle temperature steady state
        % define an initial guess
        %initial_guess(3) = 0.5*(T(1) + T(3));      % failed to converge
        initial_guess(3) = 0.5*(T(1) + T(3)) + 50;
        % solve the CSTR design equations
        [nA, ~, T(2)] = cstr_model_variables(initial_guess);
        % calculate the conversion
        fA(2) = 100*(nA_in - nA)/nA_in;

        % tabulate, show, and save the results
        resultsTable = table(fA, T - 273.15);
        resultsTable.Properties.VariableNames = ["Conversion (%)"
            "Temperature (°C)"];
        disp(' ')
        disp(resultsTable)
        disp(' ')
        writetable(resultsTable,'prac_11_results.csv')
    end

    % perform the calculations
    deliverables()
end