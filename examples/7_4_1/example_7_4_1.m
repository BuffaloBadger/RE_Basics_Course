function example_7_4_1()
% Calculations for Example 7.4.1 from REB, The Book

    % global constants available to all functions
    % given
    dH_1 = -101.2E3; % J/mol
    k0_1 = 5.11e4 * 3600; % L/mol/h
    E_1 = 74.8e3; % J/mol
    T_0 = 180 + 273.15; % K
    V = 1900.0; % L
    CA_0 = 2.9; % mol/L
    CB_0 = 3.2; % mol/L
    Cp = 1.23 * 4.184; % J/g/K
    rho = 1.02 * 1000.0; % g/L
    t_f = 2.0; % h
    % known
    R = 8.314; % J/mol/K
    % calculated
    nA_0 = CA_0*V;
    nB_0 = CB_0*V;

    % BSTR reactor function
    function [t, nA, nB, nY, nZ, T] = bstr_model_variables()
        % set the initial values
        ind0 = 0.0;
        dep0 = [nA_0, nB_0, 0.0, 0.0, T_0];

        % define the stopping criterion
        stopVar = 0;
        stopVal = t_f;
        
        % solve the IVODEs
        odes_are_stiff=false;
        [t, dep, success, message] = solve_ivodes(ind0, dep0, stopVar, stopVal...
                , @bstr_derivatives, odes_are_stiff);

        % check for solver issues
        if ~success
            disp('')
            disp(["BSTR function issue: ", message])
            disp('')
        end

        % extract the dependent variable profiles
        nA = dep(:,1);
        nB = dep(:,2);
        nY = dep(:,3);
        nZ = dep(:,4);
        T = dep(:,5);
    end

    % BSTR derivatives function
    function ddt = bstr_derivatives(~, dep)
        % extract necessary dependent variables for this integration step
        nA = dep(1);
        nB = dep(2);
        T = dep(5);

        % calculate the additional unknowns
        CA = nA/V;
        CB = nB/V;
        k = k0_1*exp(-E_1/R/T);
        r = k*CA*CB;

        % evaluate the derivatives
        dnAdt = -V*r;
        dnBdt = -V*r;
        dnYdt = V*r;
        dnZdt = V*r;
        dTdt = -r*dH_1/rho/Cp;

        % return the derivatives as a column vector
        ddt = [dnAdt; dnBdt; dnYdt; dnZdt; dTdt];
    end

    % deliverables function
    function deliverables()
        % solve the reactor design equations
        [t, nA, nB, nY, nZ, T] = bstr_model_variables();

        % calculate the other quantities of interest
        CA = nA/V;
        CB = nB/V;
        CY = nY/V;
        CZ = nZ/V;
        k = k0_1*exp(-E_1/R/T);
        r = k*CA*CB;
        T_C = T - 273.15;
        
        % display and save the graphs
        figure
        plot(t,CA,'m',t,CB,'b',t,CY,'g',t,CZ,'r:','LineWidth',2)
        xlabel("Time (h)", "FontSize", 14)
        ylabel("Concentration (mol /L)", "FontSize", 14)
        legend({"A","B","Y","Z"})
        set(gca,'FontSize',14)
        saveas(gcf,'Ci_vs_t.png')

        figure % temperature profile
        plot(t,T_C,'k', 'LineWidth', 2)
        xlabel("Time (h)","FontSize", 14)
        ylabel("Temperature (°C)","FontSize",14)
        set(gca,'FontSize',14)
        saveas(gcf,'T_vs_t.png')

        figure % rate profile
        plot(t,r, 'k', 'LineWidth', 2)
        xlabel("Time (h)","FontSize", 14)
        ylabel("Rate (mol/L/h)","FontSize",14)
        set(gca,'FontSize',14)
        saveas(gcf,'r_vs_t.png')

        end 

    % perform the calculations
    deliverables()
end