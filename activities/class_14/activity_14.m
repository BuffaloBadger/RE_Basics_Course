function activity_14()
% Calculations for the Class 14 Learning Activity from REB, The Course

    % global constants available to all functions
    % given
    V = 1200; % L
    nA_0 = 2500; % mol
    nB_0 = 50000; % mol
    nZ_0 = 0; % mol
    T_0 = 300.; % K
    fA_f = 0.5;
    k0_1 = 6e5; % /min
    E_1 = 42000; % J/mol
    dH1_0 = -16500.; % J/mol
    CpA = 180.; % J/mol/K
    CpB = 70.; % J/mol/K
    CpZ = 225.; % J/mol/K
    % known
    R = 8.3144; % J/mol/K
    % calculated
    nA_f = nA_0*(1-fA_f);

    % BSTR reactor function
    function [t, nA, nB, nZ, T] = bstr_model_variables()
        % set the initial values
        ind0 = 0;
        dep0 = [nA_0, nB_0, nZ_0, T_0];

        % set the stopping criterion
        StopVar = 1;
        StopVal = nA_f;

        % solve the design equations
        odes_are_stiff = false;
        [t, dep, success, message] = solve_ivodes(ind0, dep0, StopVar...
            , StopVal, @bstr_derivatives, odes_are_stiff);
        
        % check for solver issues
        if ~success
            disp('')
            disp(["BSTR model function issue: ", message])
            disp('')
        end

        % extract the bstr model variables to return
        nA = dep(:,1);
        nB = dep(:,2);
        nZ = dep(:,3);
        T = dep(:,4);
    end

    % BSTR derivatives function
    function ddt = bstr_derivatives(~, dep)
        % extract the dependent variables
        nA = dep(1);
        nB = dep(2);
        nZ = dep(3);
        T = dep(4);

        % calculate the additional unknowns
        k_1 = k0_1*exp(-E_1/(R*T));
        CA = nA/V;
        r_1 = k_1*CA;
        dH1 = dH1_0 + (CpZ - CpA - CpB)*(T - 298);

        % evaluate the derivatives
        dnAdt = -V*r_1;
        dnBdt = -V*r_1;
        dnZdt = V*r_1;
        dTdt = -V*r_1*dH1/(nA*CpA + nB*CpB + nZ*CpZ);

        % return the derivatives as a column vector
        ddt = [dnAdt; dnBdt; dnZdt; dTdt];
    end

    % deliverables function
    function deliverables()
        % solve the design equations
        [t, nA, ~, ~, T] = bstr_model_variables();

        % tabulate, show and save the deliverables
        item = ["Final Time"; "Final Temperature"];
        value = [t(end); T(end)];
        units = ["min"; "K"];
        resultsTable = table(item, value, units);
        disp(' ')
        disp(resultsTable)
        disp(' ')
        writetable(resultsTable,'results.csv')

        % calculate, plot, show, and save the instantaneous rate vs time
        r = k0_1*exp(-E_1./(R*T)) .* nA/V;
        figure
        plot(t,r, 'LineWidth', 2)
        xlabel('Time (min)', 'FontSize', 14)
        ylabel('Instantaneous Rate (mol/min/L)', 'FontSize', 14)
        set(gca, 'FontSize', 14)
        saveas(gcf,'r_vs_t.png')
    end

    % perform the calculations
    deliverables()
end