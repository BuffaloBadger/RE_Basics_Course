function activity_15()
% Calculations for the Class 15 Learning Activity in REB, The Course"""
    % global constants available to all functions
    % given
    CA_0 = 3.0/1000; % mol/cc
    CB_0 = 3.0/1000; % mol/cc
    T_0 = 50 + 273.15; % K
    T_max = 90 + 273.15; % K
    A = 66; % cm^2
    V_ex = 40; % cc
    U = 35*252.1*0.001076/60*1.8; % cal/cm^2/min/K
    T_ex_in = 40 + 273.15; % K
    t_turn = 30; % min
    Cp = 0.35; % cal/g/K
    rho = 0.93; % g/cc
    k0_1 = 1.24E13*60; % cc/mol/min
    E_1 = 20000; % cal/mol
    dH_1 = -80000/4.184; % cal/mol
    rho_ex = 1; % g/cc
    Cp_ex = 1; % cal/g/K
    VA_0 = 250; % cc
    T_ex_0 = 40 + 273.15; % K
    VB_0 = 250; % cc
    fA = 0.9;
    % known;
    R = 1.987; % cal/mol/K
    % calculated
    V = VA_0 + VB_0;
    nA_0 = CA_0*VA_0;
    nB_0 = CB_0*VB_0;
    nA_f = nA_0*(1 - fA);

    % define a global variable for the coolant flow rate
    g_m_ex = nan;

    % BSTR reactor model function
    function [t, nA, nB, nY, nZ, T, Tex] = bstr_model_variables(mDot_ex)
        % make the coolant flow rate available to the derivatives function
        g_m_ex = mDot_ex;

        % set the initial values
        ind_0 = 0;
        dep_0 = [nA_0, nB_0, 0, 0, T_0, T_ex_0];

        % define the stopping criterion
        stopVar = 1;
        stopVal = nA_f;

        % solve the design equations
        odes_are_stiff = false;
        [t, dep, success, message] = solve_ivodes(ind_0, dep_0, stopVar...
                , stopVal, @bstr_derivatives, odes_are_stiff);
        
        % check for solver issues
        if ~success
            disp('')
            disp(["BSTR model function issue: ", message])
            disp('')
        end
        
        % extract the dependent variables to be returned
        nA = dep(:,1);
        nB = dep(:,2);
        nY = dep(:,3);
        nZ = dep(:,4);
        T = dep(:,5);
        Tex = dep(:,6);
    end

    % BSTR derivatives function
    function ddt = bstr_derivatives(~, dep)
        % extract the dependent variables
        nA = dep(1);
        nB = dep(2);
        T = dep(5);
        Tex = dep(6);
        
        % calculate the additional unknowns
        k = k0_1*exp(-E_1/R/T);
        CA = nA/V;
        CB = nB/V;
        r = k*CA*CB;
        Q = U*A*(Tex - T);

        % calc the derivatives
        dnAdt = -r*V;
        dnBdt = -r*V;
        dnYdt = r*V;
        dnZdt = r*V;
        dTdt = (Q - r*dH_1*V)/(rho*V*Cp);
        dTexdt = (-Q - g_m_ex*Cp_ex*(Tex - T_ex_in))/(rho_ex*V_ex*Cp_ex);

        % return the derivatives as a column vector
        ddt = [dnAdt; dnBdt; dnYdt; dnZdt; dTdt; dTexdt];
    end

    % coupled unknown residual function
    function epsilon = coupled_unknown_residual(m_ex_guess)
        % solve the BSTR design equations using the guess
        [~, ~, ~, ~, ~, T, ~] = bstr_model_variables(m_ex_guess);

        % evaluate and return the residual
        epsilon = max(T) - T_max;
    end

    % deliverables function
    function deliverables()
        % guess the coolant flow rate
        m_ex_guess = 100; % g/min

        % calculate the coolant flow rate
        [soln, success, message] = solve_ates(@coupled_unknown_residual...
            , m_ex_guess);
        mDot_ex = soln(1);

        % check for solver issues
        if ~success
            disp('')
            disp(["Issue solving for the coupled unknown: ", message])
            disp('')
        end
        
        % solve the BSTR design equations
        [t, ~, ~, ~, nZ, T, ~] = bstr_model_variables(mDot_ex);

        % calculate the net rate
        net_rate = nZ/(t + t_turn);
        % generate, show, and save the requested graph
        figure
        plot(t, net_rate, 'k', 'LineWidth', 2)
        xlabel('Time (min)', 'FontSize', 14)
        ylabel('Net Rate (mol /min)','FontSize', 14)
        set(gca,'FontSize',14)
        saveas(gcf,'r_vs_t.png')

        % for discussion, solve using 90% and 110% of the base flow rate
        [t90, ~, ~, ~, nZ, T90, ~] = bstr_model_variables(0.9*mDot_ex);
        net_rate_90 = nZ/(t90 + t_turn);
        [t110, ~, ~, ~, nZ, T110, ~] = bstr_model_variables(1.1*mDot_ex);
        net_rate_110 = nZ/(t110 + t_turn);

        % plot, show and save net rate vs t
        figure
    
        l1 = [num2str(0.9*mDot_ex,3), 'g/min'];
        l2 = [num2str(mDot_ex,3), 'g/min'];
        l3 = [num2str(1.19*mDot_ex,3), 'g/min'];
        plot(t90, net_rate_90,'r',t, net_rate,'b',t110,net_rate_110,'g'...
            ,'LineWidth', 2)
        xlabel('Time (min)','FontSize', 14)
        ylabel('Net Rate (mol/min)','FontSize',14)
        set(gca,'FontSize',14)
        legend(l1, l2, l3)
        saveas(gcf,'r_vs_t_mex.png')

        % plot, show and save T vs t
        figure
        plot(t90, T90 - 273.15,'r',t, T-273.15, 'b', t110,T110 - 273.15, 'g'...
            ,'LineWidth', 2)
        xlabel('Time (min)','FontSize', 14)
        ylabel('Temperature (°C)','FontSize', 14)
        set(gca,'FontSize',14)
        legend({l1, l2, l3})
        saveas(gcf,'T_vs_t_mex.png')
    end

    % perform the calculations
    deliverables();
end    