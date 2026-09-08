function example_6_6_2()
% Calculations for Example 6.6.2 of REB, The Book
    % global constant available to all functions
    % given
    k0_1 = 10.2; % gal /mol /min
    k0_2 = 17.0; % gal /mol /min
    E_1 = 15300; % J /mol
    E_2 = 23700; % J /mol
    CA_in = 10; % mol /gal
    CB_in = 12; % mol /gal
    T_in = 350; % K
    V = 25; % gal
    Vdot_in = 12.5; % gal /min
    dH_1_298 = -12000; % J /mol
    dH_2_298 = -21300; % J /mol
    Cp_A = 85; % J /mol /K
    Cp_B = 125; % J /mol /K
    Cp_D = 200; % J /mol /K
    Cp_U = 170; % J /mol /K
    % known
    R = 8.314; % J /mol /K
    % calculated
    nA_in = CA_in*Vdot_in;
    nB_in = CB_in*Vdot_in;

    % CSTR model function
    function [nA, nB, nD, nU, T] = cstr_model_variables(init_guess)
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
        nD = soln(3);
        nU = soln(4);
        T = soln(5);
    end

    % CSTR residuals function
    function epsilon = cstr_residuals(guess)
        % extract the individual guesses
        nA = guess(1);
        nB = guess(2);
        nD = guess(3);
        nU = guess(4);
        T = guess(5);

        % calculate the additional unknowns
        k_1 = k0_1*exp(-E_1/R/T);
        k_2 = k0_2*exp(-E_2/R/T);
        CA = nA/Vdot_in;
        CB = nB/Vdot_in;
        r_1 = k_1*CA*CB;
        r_2 = k_2*CA*CB;
        dH_1 = dH_1_298 + (Cp_D - Cp_A - Cp_B)*(T - 298);
        dH_2 = dH_2_298 + (Cp_U - Cp_A - Cp_B)*(T - 298);
        epsilon_1 = nA_in - nA - V*(r_1 + r_2);
        epsilon_2 = nB_in - nB - V*(r_1 + r_2);
        epsilon_3 = -nD + V*r_1;
        epsilon_4 = -nU + V*r_2;
        epsilon_5 = (nA_in*Cp_A + nB_in*Cp_B)*(T-T_in) + V*(r_1*dH_1...
            + r_2*dH_2);

        % return the residuals as an array
        epsilon = [epsilon_1; epsilon_2; epsilon_3; epsilon_4; epsilon_5];
    end    

    % deliverables function
    function deliverables()
        % initial guess for the solver
        init_guess = [nA_in, nB_in, 0.0, 0.0, T_in + 5.0];

        % the next two lines were used to check for other steady states
%    init_guess[4] = T_in + 0.01
%    init_guess[4] = T_in + 100

        % solve the CSTR model
        [nA, ~, nD, nU, T] = cstr_model_variables(init_guess);

        % calculate the other quantities of interest
        fA = 100*(nA_in - nA)/nA_in;
        S_D_U = nD/nU;

        % tabulate, show, and save the results
        item = ["Conversion"; "Selectivity";"Temperature"];
        value = [fA; S_D_U; T];
        units = ["%"; "mol D per mol U"; "K"];
        
        resultsTable = table(item,value,units);
        disp(' ')
        disp(resultsTable)
        disp(' ')
        writetable(resultsTable,'ex_6_6_2_results.csv')
    end

    % perform the calculations
    deliverables()
end