function example_7_4_5()
% Calculations for Example 7.4.5 from REB, The Book
    % global constants available to all functions
    % given
    k01 = 1.6E18; % gal /lbmol /h
    E1 = 46000.; % BTU /lbmol
    k02 = 4.5E18; % gal /lbmol /h
    E2 = 48000.; % BTU /lbmol
    Cp = 65.0 * 0.1337; % BTU /gal /degR
    dH1 = 45000.; % BTU /lbmol
    dH2 = 39500.; % BTU /lbmol
    V = 50; % gal
    CA0 = 0.014; % lbmolA /gal
    CB0 = 0.020; % lbmolB /gal
    T0 = 70 + 459.67; % degR
    U = 60; % BTU /ft^2 /degR /h
    A = 13; % ft^2
    Tex = 220 + 459.67; % degR
    dHvap = 966; % BTU/lb
    % known
    R = 1.987; % BTU/lbmol
    % calculated
    nA0 = CA0*V;
    nB0 = CB0*V;

    % BSTR model function
    function [t, nA, nB, nD, nZ, nU, T] = bstr_model_variables(fB)
        % set the initial values
        ind_0 = 0;
        dep_0 = [nA0, nB0, 0, 0, 0, T0];

        % set the stopping criterion
        stopVar = 2;
        stopVal = nB0*(1-fB);

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
        
        % extract the bstr model variables to be returned
        nA = dep(:,1);
        nB = dep(:,2);
        nD = dep(:,3);
        nZ = dep(:,4);
        nU = dep(:,5);
        T = dep(:,6);
    end

    % BSTR derivatives function
    function ddt = bstr_derivatives(~, dep)
        % extract the dependent variables
        nA = dep(1);
        nB = dep(2);
        nD = dep(3);
        T = dep(6);

        % calculate the additional unknowns
        k1 = k01*exp(-E1/R/T);
        k2 = k02*exp(-E2/R/T);
        CA = nA/V;
        CB = nB/V;
        CD = nD/V;
        r1 = k1*CA*CB;
        r2 = k2*CD*CB;
        Q = U*A*(Tex-T);

        % evaluate the derivatives
        dnAdt = -r1*V;
        dnBdt = -(r1 + r2)*V;
        dnDdt = (r1 - r2)*V;
        dnZdt = (r1 + r2)*V;
        dnUdt = r2*V;
        dTdt = (Q - (r1*dH1 + r2*dH2)*V)/V/Cp;

        % return the derivatives as a column vector
        ddt = [dnAdt; dnBdt; dnDdt; dnZdt; dnUdt; dTdt];
    end

    % deliverables function
    function deliverables()
        % solve the BSTR design equations
        [~, ~, nB, nD, ~, ~, T] = bstr_model_variables(0.99);

        % calculate the quantities to be plotted
        fB = 100*(nB0 - nB)/nB0;
        Yield = nD/nA0;
        Qdot = U*A*(Tex - T);
        m_H2O = Qdot/dHvap;

        % generate, show, and save the graphs
        figure
        plot(fB,Yield,'k','LineWidth', 2)
        ylabel('Yield (mol D per mol A fed)','FontSize', 14)
        xlabel('Conversion of B (%)','FontSize', 14)
        set(gca, 'FontSize', 14)
        saveas(gcf, 'yield_vs_conversion.png')

        figure
        plot(fB,T - 459.67,'k','LineWidth', 2)
        ylabel('Temperature (°F)','FontSize', 14)
        xlabel('Conversion of B (%)','FontSize', 14)
        set(gca, 'FontSize', 14)
        saveas(gcf,'T_vs_conversion.png')

        figure
        plot(fB,m_H2O,'k','LineWidth', 2)
        ylabel('Water Flow Rate (lb_m /h)','FontSize', 14)
        xlabel('Conversion of B (%)','FontSize', 14)
        set(gca, 'FontSize', 14)
        saveas(gcf,'H2O_vs_conversion.png')
    end

    % perform the calculations
    deliverables()
end