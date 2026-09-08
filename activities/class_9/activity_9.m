function activity_9()
 % Calculations for the Class 9 Learning Activity in REB, The Course
    % global constants
    % given
    T_in = 150 + 273.15; % K
    CA_in = 2.0; % mol/L
    V = 500.0; % L
    Vdot_in = 250.0; % L/h
    T_ex = 180 + 273.15; % K
    A = 2.0; % m2
    U = 500E3; % cal/m2/h/K
    k01 = 1.14E9; % L/mol/h
    E1 = 16200.; % cal/mol
    Cp = 1.17E3; % cal/L/K
    dH1 = 18.2E3; % cal/mol
    % known
    R = 1.987; % cal/mol/K
    % calculated
    nA_in = CA_in*Vdot_in;

    % CSTR model function
    function [nA, nZ, T] = cstr_model_variables(init_guess)
        % solve the cstr design equations
        [soln, success, message] = solve_ates(@cstr_residuals, init_guess);

        % check for solver issues
        if ~success
            disp(' ')
            disp(['     CSTR model function issue: ', message])
            disp(' ')
        end

        % extract and return the results
        nA = soln(1);
        nZ = soln(2);
        T = soln(3);
    end

    % CSTR residuals function
    function epsilon = cstr_residuals(guess)
        % extract  the individual guesses
        nA = guess(1);
        nZ = guess(2);
        T = guess(3);

        % calculate the dditional unknowns
        Vdot = Vdot_in;
        k1 = k01*exp(-E1/(R*T));
        CA = nA/Vdot;
        r1 = k1*CA^2;
        Qdot = U*A*(T_ex - T);

        % evaluate and return the residuals
        epsilon_1 = nA_in - nA - r1*V;
        epsilon_2 = r1*V - nZ;
        epsilon_3 = Qdot - Vdot_in*Cp*(T - T_in) - r1*V*dH1;
        epsilon = [epsilon_1; epsilon_2; epsilon_3];
    end

    % deliverables function
    function deliverables()
        % define the initial guess for the cstr reactor variables
        % init_guess = [nA_in, 0.0, T_in - 100];   did not converge
        init_guess = [nA_in, 0.0, T_in - 100];

        % solve the cstr model equations
        [nA, ~, T] = cstr_model_variables(init_guess);

        % calculate the conversion
        fA = 100*(nA_in - nA)/nA_in;

        % tabulate, show, and save the results
        item = ["fA";"T"];
        value = [fA;T - 273.15];
        units = ["%";"°C"];
        resultsTable = table(item,value,units);
        disp(' ')
        disp(resultsTable)
        disp(' ')
        writetable(resultsTable,'act_9_results.csv')
    end

    % perform the calculations
    deliverables()
end