function example_6_6_3()
% Calculations for Example 6.6.3 of REB, The Book
    % constants available to all functions
    % given
    V = 0.5; % m^3
    nA_in = 70; % mol /min
    nB_in = 1500; % mol /min
    Vdot_in = 40E-3; % m^3 /min
    k0 = 1.2E9; % m^3 /mol /min
    E = 25800*4.184; % J /mol
    K0 = 4.2E-18; % m^3 /mol
    dH = -22400*4.184; % J /mol
    Cp_A = 412; % J /mol /K
    Cp_B = 75.5; % J /mol /K
    % known
    R = 8.314; % J /mol /K

    % global variable for the current value of T_in
    g_T_in = nan;

    % CSTR model function
    function [nA, nB, nZ, T] = cstr_model_variables(T_in, init_guess)
        % make T_in available to the residuals function
        g_T_in = T_in;

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
        k = k0*exp(-E/R/T);
        K = K0*exp(-dH/R/T);
        CA = nA/Vdot_in;
        CB = nB/Vdot_in;
        CZ = nZ/Vdot_in;
        r = k*CA*CB*(1 - CZ/(K*CA*CB));
        
        % evaluate and return the CSTR design equation residuals
        epsilon_1 = nA_in - nA - V*r;
        epsilon_2 = nB_in - nB - V*r;
        epsilon_3 = -nZ + V*r;
        epsilon_4 = -(nA_in*Cp_A + nB_in*Cp_B)*(T-g_T_in) - V*r*dH;

        % return the residuals as an array
        epsilon = [epsilon_1; epsilon_2; epsilon_3; epsilon_4];
    end

    % deliverables function
    function deliverables()
        % set a range for T_in
        T_in_range = (linspace(75.,125.,100) + 273.15)';

        % allocate storage for conversion and temperature at each T_in value
        fA_range = ones(100,1)*nan;
        T_range = ones(100,1)*nan;

        % set the initial guess
        init_guess = [nA_in, nB_in, 0.0, T_in_range(1) + 5.0];

        % calculate fA and T for each T_in
        for iT = 1:100

            % solve the reactor design equations
            [nA, nB, nZ, T] = cstr_model_variables(T_in_range(iT), init_guess);

            % calculate the conversion
            fA_range(iT) = 100*(nA_in - nA)/nA_in;
            T_range(iT) = T;

            % use the result as the initial guess for the next T_in
            init_guess = [nA, nB, nZ, T];
        end

        % plot the results
        figure
        plot(T_in_range - 273.15 ,fA_range ,'k','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Inlet Temperature (°C)', 'FontSize', 14)
        ylabel('Conversion of A (%)', 'FontSize', 14)
        saveas(gcf,"ex_6_6_3_fA_vs_Tin.pdf")

        % find the maximum conversion
        [fA_max, i_max] = max(fA_range);

        % calculate the other quantities of interest
        T_in_opt = T_in_range(i_max) - 273.15;
        T_out = T_range(i_max) - 273.15;

        % tabulate, show, and save the results
        item = ["opt T_in"; "fA_max"; "T_out"];
        value = [T_in_opt; fA_max; T_out];
        units = ["°C"; "%"; "°C"];
        resultsTable = table(item, value, units);
        disp(' ')
        disp(resultsTable)
        disp(' ')
        writetable(resultsTable,'ex_6_6_3_results.csv')

        % check for a high temperature steady state at the lowest T_in
        %init_guess = [nA_in, nB_in, 0.0, T_in_range(1) + 400.0];
        init_guess = [nA_in, nB_in, 0.0, T_range(1) + 15.0];
        [~, ~, ~, T] = cstr_model_variables(T_in_range(1), init_guess);
        disp('')
        disp(['Outlet T at start of range with low guess: '...
        , num2str(T_range(1),3) , ' K']);

        disp(['Outlet T at start of range with high guess: '...
        , num2str(T,3) , ' K']);
    end

    % perform the calculations
    deliverables()
end