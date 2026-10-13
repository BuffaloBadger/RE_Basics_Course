function activity_18()
% Calculations for the Class 18 Learning Activity from REB, The Course
    %global constants available to all functions
    %given
    dH1 = -16600; %cal /mol
    CA_0 = 2; %mol /L
    T0 = 20 + 273.15; %K
    k01 = 2.5E8; %L /mol /min
    E1 = 14300; %cal /mol
    V = 0.4; %L
    V_ex = 100; %cm^3
    A_ex = 32.5; %cm^2
    U_ex = 0.2; %cal /min /cm^2 /K
    Tex0 = 10 + 273.15; % K
    Tex_in = 10 + 273.15; %K
    rho_ex = 1; %g /cm^3
    Cp_ex = 1; %cal /g /°C
    mDot0 = 217; %g /min
    d_mDot = 212; %g /min
    tf = 60; %min
    Cp = 0.42; %cal /g /°C
    rho = 879; %g /L
    %known
    R = 1.987; %cal / mol / K
    %calculated
    nA0 = CA_0*V;
    
    % given equation
    mDot = @(t) d_mDot/tf * t + mDot0;

    %BSTR reactor function
    function [t, nA, nZ, T, Tex] = bstr_model_variables()
        %set the initial values
        ind_0 = 0;
        dep_0 = [nA0, 0, T0, Tex0];

        %set the stopping criterion
        stopVar = 0;
        stopVal = tf;

        %solve the design equations
        odes_are_stiff=true;
        [t, dep, success, message] = solve_ivodes(ind_0, dep_0, stopVar...
            , stopVal, @bstr_derivatives, odes_are_stiff);
        
        %check for solver issues
        if ~success
            disp('')
            disp(["BSTR model issue: ", message])
            disp('Press any key to continue or CTRL-C to exit')
            pause
        end
        
        % return the bstr model variables
        nA = dep(:,1);
        nZ = dep(:,2);
        T = dep(:,3);
        Tex = dep(:,4);
    end

    % BSTR derivatives function
    function ddt = bstr_derivatives(t, dep)
        %extract the dependent variables
        nA = dep(1);
        T = dep(3);
        Tex = dep(4);

        %calculate the additional unknowns
        k1 = k01 * exp(-E1/(R*T));
        CA = nA / V;
        r1 = k1 * CA;
        Q = U_ex * A_ex * (Tex - T);

        %evaluate the derivatives
        dnA_dt = -r1 * V;
        dnZ_dt = r1 * V;
        dT_dt = (Q - r1*V*dH1) / (rho*V*Cp);
        dTex_dt = (-Q - mDot(t)*Cp_ex*(Tex - Tex_in)) / (rho_ex*V_ex*Cp_ex);

        %return the derivatives as a column vector
        ddt = [dnA_dt; dnZ_dt; dT_dt; dTex_dt];
    end

    % items of interest function
    function items_of_interest()
        %solve the BSTR design equations
        [t, nA, ~, T, Tex] = bstr_model_variables();

        %tabulate, show, and save results
        item = ["Conversion"; "Maximum T"; "Minimum T"];
        value = [100*(nA0 - nA(end))/nA0; max(T) - 273.15; min(T) - 273.15];
        units = ["%"; "°C"; "°C"];
        resultsTable = table(item, value, units);
        disp('')
        disp(resultsTable)
        disp('')
        writetable(resultsTable, 'results.csv');

        figure
        plot(t, T-273.15, t, Tex-273.15, 'LineWidth', 2)
        xlabel('Time (min)', 'FontSize', 14)
        ylabel('Temperature (°C)','FontSize', 14)
        legend({"Reactor","Cooling Water"})
        set(gca, 'FontSize', 14)
        saveas(gcf,'T_profiles.png')

        figure
        plot(t, 100*(nA0-nA)/nA0, 'LineWidth', 2)
        xlabel('Time (min)','FontSize', 14)
        ylabel('Conversion of A (%)','FontSize', 14)
        saveas(gcf,'fA_profile.png')
    end

    % perform the calculations
    items_of_interest();
end
    