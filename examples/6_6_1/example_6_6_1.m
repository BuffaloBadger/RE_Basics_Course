function example_6_6_1()
% Calculations for Example 6.6.1 of REB, The Book.
    % global constants available to all functions
    yA_in = 0.1;
    yB_in = 0.65;
    yI_in = 0.25;
    T_in = 165 + 273.15; % K
    P = 5.0; % atm
    yA = 0.001;
    k0 = 1.37E5; % m^3 /mol /min
    E = 11100.0; % cal /mol
    dH = -7200.0; % cal /mol
    Cp_A = 7.6; % cal /mol /K
    Cp_B = 8.2; % cal /mol /K
    Cp_I = 4.3; % cal /mol /K
    % known
    Re = 1.987; % cal /mol /K
    Rw = 8.206E-5; % m^3 atm /mol /K
    % basis
    Vdot_in = 1.0; % m^3 /min
    % calculated
    nA_in = yA_in*P*Vdot_in/Rw/T_in;
    nB_in = yB_in*P*Vdot_in/Rw/T_in;
    nI_in = yI_in*P*Vdot_in/Rw/T_in;

    % CSTR model function
    function [V, nB, nI, nZ, T] = cstr_model_variables(init_guess)
        % solve the design equations
        [soln, success, message] = solve_ates(@cstr_residuals, init_guess);

        % check for solver issues
        if ~success
            disp(' ')
            disp(['    CSTR model function issue: ', message])
        end
        
        % extract and return the results
        V = soln(1);
        nB = soln(2);
        nI = soln(3);
        nZ = soln(4);
        T = soln(5);
    end

    % CSTR residuals function
    function epsilon = cstr_residuals(guess)
        % extract the individual guesses
        V = guess(1);
        nB = guess(2);
        nI = guess(3);
        nZ = guess(4);
        T = guess(5);

        % calculate the additional unknowns
        nA = yA*(nB + nI + nZ)/(1 - yA);
        k = k0*exp(-E/Re/T);
        CA = nA/(nA + nB + nI + nZ)*P/Rw/T;
        CB = nB/(nA + nB + nI + nZ)*P/Rw/T;
        r = k*CA*CB;
        epsilon_1 = nA_in - nA - V*r;
        epsilon_2 = nB_in - nB - V*r;
        epsilon_3 = nI_in - nI;
        epsilon_4 = -nZ + V*r;
        epsilon_5 = -(nA_in*Cp_A + nB_in*Cp_B + nI_in*Cp_I)*(T - T_in) - V*r*dH;

        % return the residuals as an array
        epsilon = [epsilon_1; epsilon_2; epsilon_3; epsilon_4; epsilon_5];
    end

    % deliverables function
    function deliverables()
        % set the initial guess
        init_guess = [1.0, nB_in, nI_in, 0.0, T_in + 5.0];

        % solve the reactor design equations
        [V, ~, ~, ~, T] = cstr_model_variables(init_guess);

        % calculate the space time
        tau = V/Vdot_in;

        % tabulate, show, and save the results
        item = ["tau";"T"];
        value = [tau;T - 273.15];
        units = ["min";"°C"];
        resultsTable = table(item,value,units);
        disp(' ')
        disp(resultsTable)
        disp(' ')
        writetable(resultsTable,'ex_6_6_1_results.csv')
    end

    % perform the calculations
    deliverables()
end
