function practice_11_multiplicity()
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

    % global variable for the current value of T
    g_T = nan;

    % CSTR model function
    function [nA, nZ, Tin] = cstr_model_variables(T, init_guess)
        % make T available to the residuals function
        g_T = T;
        
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
        Tin = soln(3);
    end

    % CSTR residuals function
    function epsilon = cstr_residuals(guess)
        % extract the individual guesses
        nA = guess(1);
        nZ = guess(2);
        Tin = guess(3);

        % calculate the additional unknowns
        k1 = k01*exp(-E1/(R*g_T));
        CA = nA/Vdot;
        r1 = k1*CA;

        % evaluate the residuals
        epsilon_1 = nA_in - nA - V*r1;
        epsilon_2 = nZ_in - nZ + V*r1;
        epsilon_3 = Vdot_in*Cp*(g_T - Tin) + V*r1*dH1;

        % return the residuals as an array
        epsilon = [epsilon_1; epsilon_2; epsilon_3];
    end

    % deliverables function
    function deliverables()
        % set a range of outlet temperatures
        %T_range = linspace(25, 400, 100) + 273.15;
        T_range = linspace(25, 250, 100) + 273.15;

        % allocate storage for the corresponding inlet temperatures
        Tin_range = ones(1,100)*nan;

        % set an initial guess for the first outlet temperature in the range
        init_guess = [nA_in, nZ_in, 0.0, T_range(1) - 5.0];

        % calculate the inlet temperature for each outlet temperature
        for i = 1:100
            % solve the design equations
            [nA, nZ, Tin_range(i)] = cstr_model_variables(T_range(i)...
                , init_guess);

            % use the result as the next initial guess
            init_guess = [nA, nZ, Tin_range(i)];
        end

        % generate a multiplicity plot
            figure;
            plot(Tin_range - 273.15, T_range - 273.15,'k','LineWidth',2)
            xline(T_in - 273.15, 'b','Linewidth',2)
            set(gca, 'FontSize', 14);
            xlabel('Inlet Temperature (°C)','FontSize', 14)
            ylabel('Outlet Temperature (°C)','FontSize', 14)
            saveas(gcf,"prac_11_multiplicity_plot.pdf")
    end

    % perform the calculations
    deliverables()
end