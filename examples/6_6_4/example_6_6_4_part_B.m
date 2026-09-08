function example_6_6_4_part_B()
% Calculations for Example 6.6.4 from REB, the Book
    % given constants available to all functions 
    V = 500; % cm^3
    Vdot_in = 1.0; % cm^3 /s
    CA_in = 0.015; % mol /cm^3
    CB_in = 0.015; % mol /cm^3
    Tin = 50 + 273.15; % K
    Cp = 0.35*4.184; % J /g /K
    rho = 0.93; % g /cm^3
    dH = -20000; % J /mol
    k0 = 3.24E12; % cm^3 /mol /s
    E = 105000; % J /mol
    % known
    R = 8.314; % J /mol /K
    % calculated
    nA_in = Vdot_in*CA_in;
    nB_in = Vdot_in*CB_in;

    % global variable for the current value of T
    T = nan;

    % CSTR model function
    function [nA, nB, nY, nZ, T] = cstr_model_variables(init_guess)
        % solve the design equations
        [soln, success, message] = solve_ates(@cstr_residuals, init_guess);

        % check for solver issues
        if ~success
            disp(' ')
            disp(['    CSTR model function issue: ', message])
        end
        
        % extract and return the results
        nA = soln(1);
        nB = soln(2);
        nY = soln(3);
        nZ = soln(4);
        T = soln(5);
    end

    % CSTR residuals function
    function epsilon = cstr_residuals(guess)
        % extract the individual guesses
        nA = guess(1);
        nB = guess(2);
        nY = guess(3);
        nZ = guess(4);
        T = guess(5);

        % calculate the additional unknowns
        k = k0*exp(-E/R/T);
        Vdot = Vdot_in;
        CA = nA/Vdot;
        CB = nB/Vdot;
        r = k*CA*CB;

        % evaluate the residuals
        epsilon_1 = nA_in - nA - V*r;
        epsilon_2 = nB_in - nB - V*r;
        epsilon_3 = -nY + V*r;
        epsilon_4 = -nZ + V*r;
        epsilon_5 = -Vdot_in*rho*Cp*(T - Tin) - V*r*dH;

        % return the residuals as an array
        epsilon = [epsilon_1; epsilon_2; epsilon_3; epsilon_4; epsilon_5];
    end

    % deliverables function
    function deliverables()
        % allocate storage for the results
        ss_fA = [nan; nan; nan];
        ss_T = [nan; nan; nan];

        % set an initial guess for the low T steady state
        init_guess = [nA_in, nB_in, 0.0, 0.0, 50.0 + 273.15];

        % solve the reactor design equations
        [nA, ~, ~, ~, ss_T(1)] = cstr_model_variables(init_guess);

        % calculate the conversion
        ss_fA(1) = 100*(nA_in - nA)/nA_in;

        % set an initial guess for the middle T steady state
        init_guess(5) = 140 + 273.15;

        % solve the reactor design equations
        [nA, ~, ~, ~, ss_T(2)] = cstr_model_variables(init_guess);

        % calculate the conversion
        ss_fA(2) = 100*(nA_in - nA)/nA_in;

        % set an initial guess for the high T steady state
        init_guess(5) = 270 + 273.15;

        % solve the reactor design equations
        [nA, ~, ~, ~, ss_T(3)] = cstr_model_variables(init_guess);

        % calculate the conversion
        ss_fA(3) = 100*(nA_in - nA)/nA_in;

        % tabulate, show, and save the results
        resultsTable = table(ss_fA, ss_T - 273.15);
        resultsTable.Properties.VariableNames = ["Conversion (%)"
            "Temperature (°C)"];
        disp(' ')
        disp(resultsTable)
        disp(' ')
        writetable(resultsTable,'ex_6_6_4_steady_states.csv')
    end

    % perform the calculations
    deliverables()
end