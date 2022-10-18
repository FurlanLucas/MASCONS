function out = gaussTable(alpha, beta, gamma, myPrint)
    %% Gauss table calculation function
    %
    % Function to gauss coeficients and wights calculation. Returns a table
    % with the values of the points x_m, y_m, z_m and the weights w_m. The
    % function are based in the roots of legendre plynomials, wich are
    % defined in [-1, 1]. The three firsts inputs are the number of points
    % in x, y and z domains, respectivel, having a maximum value of 10. The
    % fourth parameter is an optional boolean variable wich enables the
    % results printing in the command windown.
    %
    % See also legendreTable.
    
    %%
    if ~exist('myPrint', 'var')
        myPrint = false;
    end
    
    %% Main
    % Get legendre table
    table_x = legendreTable(alpha);
    table_y = legendreTable(beta);
    table_z = legendreTable(gamma);

    cm = zeros(alpha*beta*gamma, 1);
    xm = zeros(alpha*beta*gamma, 1);
    ym = zeros(alpha*beta*gamma, 1);
    zm = zeros(alpha*beta*gamma, 1);

    % Calculation
    m = 1;
    for i = 1:alpha
        for j = 1:beta
            for k = 1:gamma
                cm(m) = ((1-table_x(i,1))^2)*(1-table_y(j,1))*( ...
                    table_x(i,2)*table_y(j,2)*table_z(k,2))/64;
                xm(m) = (1+table_x(i,1))/2;
                ym(m) = (1-table_x(i,1))*(1+table_y(j,1))/4;
                zm(m) = (1-table_x(i,1))*(1-table_y(j,1))*...
                    (1+table_z(k,1))/8;
                if myPrint
                    fprintf("%0.15f \t\t %0.15f \t\t %0.15f \t\t" + ...
                        " %0.15f \n", xm(m), ym(m), zm(m), cm(m));
                end
                m = m + 1;
            end
        end
    end
    
    out = [xm, ym, zm, cm];
end

