function activity_11_multiplicity()
% Calculations for the Class 11 Learning Activity from REB, the Course
    % global constants available to all functions 
    % given
    tau = 80.; % s
    P = 3.; % atm
    T_in = 260. + 273.15; % K
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

    % global variable for the current value of T
    g_T = nan;

    % CSTR model function
    function [nA, nB, nZ, Tin] = cstr_model_variables(T, init_guess)
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
        nB = soln(2);
        nZ = soln(3);
        Tin = soln(4);
    end

    % CSTR residuals function
    function epsilon = cstr_residuals(guess)
        % extract the individual guesses
        nA = guess(1);
        nB = guess(2);
        nZ = guess(3);
        Tin = guess(4);

        % calculate the additional unknowns
        nA_in = 0.5*P*Vdot_in/Rpv/Tin;
        nB_in = 0.5*P*Vdot_in/Rpv/Tin;
        vDot = (nA + nB + nZ)*Rpv*g_T/P;
        CA = nA/vDot;
        CB = nB/vDot;
        k1 = k01*exp(-E1/(R*g_T));
        r1 = k1*CA*CB;

        % evaluate the residuals
        epsilon_1 = nA_in - nA - V*r1;
        epsilon_2 = nB_in - nB - V*r1;
        epsilon_3 = -nZ + V*r1;
        epsilon_4 = (nA*CpA + nB*CpB + nZ*CpZ)*(g_T - Tin) + V*r1*dH1;

        % return the residuals as an array
        epsilon = [epsilon_1; epsilon_2; epsilon_3; epsilon_4];
    end

    % deliverables function
    function deliverables()
        % set a range of outlet temperatures
        %T_range = linspace(260, 600, 100) + 273.15;
        T_range = linspace(275, 450, 100) + 273.15;

        % set an initial guess
        T_in_guess = T_range(1) - 5;
        nA_in_guess = 0.5*P*Vdot_in/Rpv/T_in;
        nB_in_guess = 0.5*P*Vdot_in/Rpv/T_in;
        init_guess = [nA_in_guess, nB_in_guess, 0.0, T_in_guess];

        % allocate storage for the corresponding inlet temperatures
        Tin_range = ones(1,100)*nan;

        % calculate the inlet temperature for each outlet temperature
        for i = 1:100
            % solve the design equations
            [nA, nB, nZ, Tin_range(i)] = cstr_model_variables(T_range(i), init_guess);

            % use the result as the next initial guess
            init_guess = [nA, nB, nZ, Tin_range(i)];
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
            saveas(gcf,"act_11_multiplicity_plot.png")
    end

    % perform the calculations
    deliverables()
end