function practice_12()
% Calculations for the Class 12 Practice Assignment from REB, The Course
    % global constants available to all functions
    % given
    V = 1.5; % L
    T_in = 50 + 273.15; % K
    CA_in = 1.5; % mol/L
    CB_in = 2.5; % mol/L
    k01 = 1.4E14; % /min
    E1 = 22000; % cal/mol
    dH1 = -211000*0.239; % cal/mol
    Cp = 1.3*1000; % cal/L/K
    tau = 10; % min
    T_0 = T_in;
    % known
    R = 1.987; % cal/mol/K
    % calculated
    Vdot_in = V/tau;
    Vdot = Vdot_in;
    nA_in = CA_in*Vdot_in;
    nB_in = CB_in*Vdot_in;

    % cstr model function
    function [t, nA, nB, nY, nZ, T] = cstr_model_variables()
        % set the initial values
        ind_0 = 0.0;
        dep_0 = [0, 0, 0, 0, T_0];

        % set the stopping criterion
        stopVar = 0;
        stopVal = 60; % min

        % solve the cstr design equations
        odes_are_stiff = false;
        [t, dep, success, message] = solve_ivodes(ind_0, dep_0, stopVar...
            , stopVal,@cstr_derivatives, odes_are_stiff);
        
        % check for solver issues
        if ~success
            disp('')
            disp(['CSTR model function issue: ', message])
            disp('')
        end

        % return the individual cstr model variables
        nA = dep(:,1);
        nB = dep(:,2);
        nY = dep(:,3);
        nZ = dep(:,4);
        T = dep(:,5);
    end

    % cstr derivatives function
    function ddt = cstr_derivatives(~, dep)
        % extract the dependent variables
        nA = dep(1);
        nB = dep(2);
        nY = dep(3);
        nZ = dep(4);
        T = dep(5);

        % calculate the additional unknowns
        k1 = k01*exp(-E1/(R*T));
        CA = nA/Vdot;
        r1 = k1*CA;

        % evaluate the derivatives
        dnAdt = Vdot/V*(nA_in - nA - 2*r1*V);
        dnBdt = Vdot/V*(nB_in - nB - r1*V);
        dnYdt = Vdot/V*(-nY + 2*r1*V);
        dnZdt = Vdot/V*(-nZ + 2*r1*V);
        dTdt = -(Vdot*Cp*(T - T_in) + r1*V*dH1)/(V*Cp);

        % return the derivatives as a column vector
        ddt = [dnAdt; dnBdt; dnYdt; dnZdt; dTdt];
    end

    % deliverables function
    function deliverables()
        % solve the cstr design equations
        [t, ~, ~, nY, ~, T] = cstr_model_variables();

        % generate, show, and save the requested graphs
        figure;
        plot(t, nY, 'k', 'LineWidth', 2)
        xlabel('Time (min)','FontSize', 14)
        ylabel('Outlet Molar Flow of Y (mol /min)','FontSize', 14)
        set(gca, 'FontSize', 14)
        saveas(gcf, 'nY_vs_t.png')

        figure;
        plot(t, T - 273.15, 'k', 'LineWidth', 2)
        xlabel('Time (min)','FontSize', 14)
        ylabel('Temperature (°C)','FontSize', 14)
        set(gca, 'FontSize', 14);
        saveas(gcf, 'T_vs_t.png')
    end

    % perform the calculations
    deliverables();
end
