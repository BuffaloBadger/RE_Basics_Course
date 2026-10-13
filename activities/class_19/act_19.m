function act_19()
% Calculations for the Class 19 Learning Activity from REB, The Course
    % global constants available to all functions
    % given
    dH_1 = -184000.0; % J/mol
    k0_1 = 4.866E14; % L/mol/min
    E_1 = 74100.0; % J/mol
    P = 1.0; % atm
    V0 = 4.0; % L
    CA0 = 10.0; % mol /L
    T0 = 20 + 273.15; % K
    Tex0 = 20 + 273.15; % K
    Vex = 500.0/1.0E3; % L
    Aex = 550.0/1.0E4; % m^2
    U = 8480.0; % J/m^2/min/K
    Vdot_ex = 1.0; % L/min
    Tex0 = 20 + 273.15; % K
    Cp = 4.184; % J/g/K
    rho = 1.0E3; % g /L
    CBin = 10.0; % mol/L
    Tex_in = 20 + 273.15; % K
    Tin = 20 + 273.15; % K
    Tmax = 80.0 + 273.15; % K
    % known
    R = 8.3136; % J/mol/K
    Rpv = 8.206E-2; % L-atm/mol/K
    % calculated
    nA0 = CA0*V0;
    nAf = 0.01*nA0;

    % global variable for the current feed rate
    g_Vdot = nan;

    % SBSTR reactor function
    function [t, nA, nB, nZ, T, Tex, V] = sbstr_model_variables(Vdot)
        % make Vdot available to the derivatives function
        g_Vdot = Vdot;

        % define the initial values
        ind_0 = 0;
        dep_0 = [nA0, 0, 0, T0, Tex0, V0];

        % define the stopping criterion
        stopVar = 1;
        stopVal = nAf;

        % solve the SBSTR design equations
        odes_are_stiff = true;
        [t, dep, success, message] = solve_ivodes(ind_0, dep_0, stopVar...
            , stopVal, @sbstr_derivatives, odes_are_stiff);
        
        % check for solver issues
        if ~success
            disp('')
            disp(["SBSTR model function issue: ", message])
            disp('Press any key to continue or Ctrl-C to exit.')
            pause
        end
        
        % return the sbstr model variables as separate arrays
        nA = dep(:,1);
        nB = dep(:,2);
        nZ = dep(:,3);
        T = dep(:,4);
        Tex = dep(:,5);
        V = dep(:,6);
    end

    % SBSTR derivatives function
    function ddt = sbstr_derivatives(~, dep)
        % extract the dependent variables
        nA = dep(1);
        nB = dep(2);
        T = dep(4);
        Tex = dep(5);
        V = dep(6);

        % calculate the additional unknowns
        nBin = g_Vdot*CBin;
        k1 = k0_1*exp(-E_1/R/T);
        CA = nA/V;
        CB = nB/V;
        r1 = k1*CA*CB;
        Qdot = U*Aex*(Tex - T);
        mDot_ex = Vdot_ex*rho;

        % evaluate the derivatives
        dnAdt = -r1*V;
        dnBdt = nBin - r1*V;
        dnZdt = 2*r1*V;
        dTdt = (Qdot - g_Vdot*Cp*rho*(T - Tin) - r1*V*dH_1...
             + P*g_Vdot*R/Rpv)/(V*Cp*rho);
        dTexdt = -(Qdot + mDot_ex*Cp*(Tex - Tex_in))/(rho*Vex*Cp);
        dVdt = g_Vdot;

        % return the derivatives as a column vector
        ddt = [dnAdt; dnBdt; dnZdt; dTdt; dTexdt; dVdt];
    end

    % coupled unknown residual function
    function epsilon = coupled_unknown_residual(guess)
        % solve the sbstr design equations using the guess
        [~, ~, ~, ~, T, ~, ~] = sbstr_model_variables(guess);

        % evaluate and return the residual
        epsilon = max(T) - Tmax;
    end

    % items of interest function
    function items_of_interest()
        % guess the optimum flow rate
        % Vdot_guess = 1.0; % failed to converge
        Vdot_guess = 0.01;

        % solve the implicit equation for the optimum flow rate
        [soln, success, message] = solve_ates(@coupled_unknown_residual...
            , Vdot_guess);
        Vdot_opt = soln;

        % check for solver issues
        if ~success
            disp('')
            disp(["Issue solving the implicit ATE: ", message])
            disp('Press any key to continue or Ctrl-C to exit.')
            pause
        end
        
        % solve the sbstr design equations using the optimum Vdot
        [t, nA, nB, ~, T, Tex, V] = sbstr_model_variables(Vdot_opt);

        % find the maximum T and concentration of B
        maxT = max(T) - 273.15;
        max_CB = max(nB./V);

        % tabulate, show, and save the results
        item = ["Optimum Feed Rate"; "Maximum Temperature"...
            ; "Maximum Concentration of B"];
        value = [Vdot_opt; maxT; max_CB];
        units = ["L/min"; "°C"; "M"];
        resultsTable = table(item, value, units);
        disp('')
        disp(resultsTable)
        disp('')
        writetable(resultsTable, 'act_19_results.csv')

        % for discussion, plot fA and T vs t
        figure
        plot(t,100*(nA0 - nA)/nA0, 'LineWidth', 2)
        xlabel('Time (min)','FontSize', 18)
        ylabel('Conversion of A (%)','FontSize', 18)
        set(gca,'FontSize', 18)
        saveas(gcf, 'act_19_fA_vs_t.png')

        figure
        hold on
        plot(t, T-273.15, 'LineWidth', 2)
        plot(t,Tex-273.15, 'LineWidth', 2)
        hold off
        xlabel('Time (min)', 'FontSize', 18)
        ylabel('Temperature (°C)','FontSize', 18)
        legend({'Reacting Fluid', 'Cooling Water'}, Location='best')
        set(gca,'FontSize', 18)
        saveas(gcf, 'act_19_T_vs_t.png')
    end

    % perform the calculations
    items_of_interest()  
end