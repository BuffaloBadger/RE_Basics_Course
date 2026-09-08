function activity_11_steady_states()
% Calculations for the Class 11 Learning Activity in REB, The Course
    % global constants available to all functions 
    % given
    T_in = 260. + 273.15; % K
    P = 3.; % atm
    tau = 80.; % s
    k01 = 1.26e6; %L /mol /s
    E1 = 19000.; % cal /mol
    dH1 = -10500.; % cal /mol
    CpA = 42/4.184; % cal /mol /K
    CpB = 122/4.184; % cal /mol /K
    CpZ = 173/4.184; % cal /mol /K
    % known
    R = 1.987; % cal/mol
    Rpv = 0.08206; % L atm /mol /K
    % basis
    V = 1.; % L
    % calculated
    Vdot_in = V/tau;
    nA_in = 0.5*P*Vdot_in/Rpv/T_in;
    nB_in = 0.5*P*Vdot_in/Rpv/T_in;

    % CSTR model function
    function [nA, nB, nZ, T] = cstr_model_variables(init_guess)
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
        nZ = soln(3);
        T = soln(4);
    end

    % CSTR residuals function
    function epsilon = cstr_residuals(guess)
        % extract the individual guesses
        nA = guess(1);
        nB = guess(2);
        nZ = guess(3);
        T = guess(4);

        % calculate the additional unknowns
        vDot = (nA + nB + nZ)*Rpv*T/P;
        CA = nA/vDot;
        CB = nB/vDot;
        k1 = k01*exp(-E1/(R*T));
        r1 = k1*CA*CB;

        % evaluate the residuals
        epsilon_1 = nA_in - nA - V*r1;
        epsilon_2 = nB_in - nB - V*r1;
        epsilon_3 = -nZ + V*r1;
        epsilon_4 = (nA*CpA + nB*CpB + nZ*CpZ)*(T - T_in) + V*r1*dH1;

        % return the residuals as an array
        epsilon = [epsilon_1; epsilon_2; epsilon_3; epsilon_4];
    end

    % deliverables function
    function deliverables()
        % allocate storage for the requested results
        fA = ones(3,1)*nan;
        T = ones(3,1)*nan;

        % low temperature steady state
        % define an initial guess
        initial_guess = [nA_in, nB_in, 0.0, 280 + 273.15];
        % solve the CSTR design equations
        [nA, ~, ~, T(1)] = cstr_model_variables(initial_guess);
        % calculate the conversion
        fA(1) = 100*(nA_in - nA)/nA_in;

        % repeat for the middle temperature steady state
        % define an initial guess
        initial_guess(4) = 350 + 273.15;
        % solve the CSTR design equations
        [nA, ~, ~, T(2)] = cstr_model_variables(initial_guess);
        % calculate the conversion
        fA(2) = 100*(nA_in - nA)/nA_in;

        % repeat for the high temperature steady state
        % define an initial guess
        %initial_guess(4) = 425 + 273.15;     failed to converge
        initial_guess(4) = 500 + 273.15;
        % solve the CSTR design equations
        [nA, ~, ~, T(3)] = cstr_model_variables(initial_guess);
        % calculate the conversion
        fA(3) = 100*(nA_in - nA)/nA_in;

        % tabulate, show, and save the results
        resultsTable = table(fA, T - 273.15);
        resultsTable.Properties.VariableNames = ["Conversion (%)"
            "Temperature (°C)"];
        disp(' ')
        disp(resultsTable)
        disp(' ')
        writetable(resultsTable,'act_11_results.csv')
    end

    % perform the calculations
    deliverables()
end