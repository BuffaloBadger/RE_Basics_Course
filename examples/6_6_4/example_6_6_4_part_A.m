function example_6_6_4_part_A()
% Calculations for Example 6.6.4 from REB, the Book
    % given constants available to all functions 
    V = 500; % cm^3
    Vdot_in = 1.0; % cm^3 /s
    CA_in = 0.015; % mol /cm^3
    CB_in = 0.015; % mol /cm^3
    T_in = 50 + 273.15; % K
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
    function [nA, nB, nY, nZ, Tin] = cstr_model_variables(init_guess)
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
        Tin = soln(5);
    end

    % CSTR residuals function
    function epsilon = cstr_residuals(guess)
        % extract the individual guesses
        nA = guess(1);
        nB = guess(2);
        nY = guess(3);
        nZ = guess(4);
        Tin = guess(5);

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
        % set a range of outlet temperatures
        T_range = (linspace(-25.0, 325.0, 350) + 273.15)';

        % set an initial guess
        init_guess = [nA_in, nB_in, 0.0, 0.0, T_range(1) - 1.0];

        % allocate storage for the corresponding inlet temperatures
        Tin_range = ones(100,1)*nan;

        for i = 1:350
            % set the current outlet temperature
            T = T_range(i);

            % solve the design equations
            [nA, nB, nY, nZ, Tin_range(i)] = cstr_model_variables(init_guess);

            % use the result as the next initial guess
            init_guess = [nA, nB, nY, nZ, Tin_range(i)];
        end

        % convert to °C
        T_range = T_range - 273.15;
        Tin_range = Tin_range - 273.15;

        % generate a multiplicity plot
            figure;
            plot(Tin_range, T_range,'k','LineWidth',2)
            xline(T_in - 273.15, 'b','Linewidth',2)
            set(gca, 'FontSize', 14);
            xlabel('Inlet Temperature (°C)','FontSize', 14)
            ylabel('Outlet Temperature (°C)','FontSize', 14)
            saveas(gcf,"multiplicity_plot.png")
    end

    % perform the calculations
    deliverables()
end