% Alan Perez
% Project

function[ ] = MAT4800_Project( )
%   The project function will apply Newton's Method and the
% Broyden-Fletcher-Glodfarb-Shanmo (BFGS) Method to a data-fitting problem
% and to present the analytical and numerical results. 
% 
%   Given length-K vectors of observed data, together forming K data points
% of the form ( xi, yi ) with 1 <= i <= k, we wish to "fit" a sinusoidal
% funciton to these data points be determining a vector of parameters *a*
% which MINIMIZES the mean squared error between the given data points and
% the values produced by the function. 
%   The nonlinear optimization problem is to find vector *a* to minimize
% our objective function. 

    % Variables
    %    a         Vector of parameters
    %    K         Number of data points in the lists
    %    x         The x value in the pair of data points
    %    y         The y value in the pair of data points
    %    afXData   The x values of the data we will try fitting a sinusoid to
    %    afYData   The y values of the data we will try fitting a sinusoid to
    
    % Load the data into MATLAB
    afXData = load( 'xData.txt' );
    afYData = load( 'yDataNoisy.txt' ); 
    
    % Matrix containing figure values
    aaFigures = [ 1 2 3 4 ; 
                  5 6 7 8 ; 
                  9 10 11 12 ];
    
    % Matrix containing experiments 1, 2, and 3 as its columns
    aaExperiments = [ 1.23249720 1.90064870 1.82428450 ; 
                      0.14844971 2.79581380 2.36446250 ; 
                      3.01837270 2.80402290 2.86998480 ; 
                      0.51600794 2.92119350 0.93243968 ];
    
    % Prompt user for which experiment they want to run
    fprintf("Select which experiment you would like to run ")
    experimentinput = input("(1, 2, 3, 4, 5, or 6 (stop)): "); 
    
    % For Experiments 1, 2, and 3
    while experimentinput ~= 6
        if experimentinput == 1 || experimentinput == 2 || experimentinput == 3
            % Vector of figure values
            aFigures = aaFigures( experimentinput, : );

            % Make plots for Newton's method
            achMethod = 'Newtons Method'; 
            [ aaParams, aError, k ] = newtons_method( aaExperiments( :, experimentinput ), afXData, afYData );

            % Iterative fit functions
            fig = aFigures( 1 ); 
            color = '-g';
            plot_iterative_fit( afXData, afYData, aaParams, k, fig, achMethod, experimentinput, color );

            % Error Estimates versus iterations
            fig = aFigures( 2 ); 
            plot_error_estimates( aError, k, fig, achMethod, experimentinput )

            % PrintFinal a(k), corresponding k, and corresponding f(a, x, y) value
            fprintf("Experiment %d: %s Results\n", experimentinput, achMethod)
            print_results( aaParams, k-1, afXData, afYData )

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Make plots for BFGS Method
            achMethod = 'BFGS Method'; 
            [ aaParams, aError, k ] = BFGS_method( aaExperiments( :, experimentinput ), afXData, afYData );

            % Iterative fit plots
            fig = aFigures( 3 ); 
            color = '-g';
            plot_iterative_fit( afXData, afYData, aaParams, k, fig, achMethod, experimentinput, color );

            % Error Estimates plots
            fig = aFigures( 4 ); 
            plot_error_estimates( aError, k, fig, achMethod, experimentinput )

            % Print Final a(k), corresponding k, and corresponding f(a, x, y) value
            fprintf("Experiment %d: %s Results\n", experimentinput, achMethod)
            print_results( aaParams, k-1, afXData, afYData )


        elseif experimentinput == 4
            fprintf("Experiment 4")
            run_experiment_4( afXData, afYData );

        elseif experimentinput == 5
            fprintf("Experiment 5\n")
            run_experiment_5( afXData, afYData ); 
        end 
        
        fprintf("If you would like to run another experiment, select which one")
        experimentinput = input(" (1, 2, 3, 4, 5, or 6 (stop)): ");
    end
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function[ fFunctionAta ] = func( a, x, y )
% func will take in the initial guess for the parameter vector. It will
% then evaluate the function with these parameters
% Input
%
% Output
%

    % Local Variables
    %   K    length of the x and y vectors, how many data points there are
    
    % Get length of the data (how many data points there are)
    K = length( x ); 
    
    % Initialize sum 
    fFunctionAta = 0; 
    
    % anonymous function f
    f = @( a0, a1, a2, a3, yi, xi ) (yi - a0 - a1*sin(a2*xi+a3))^2; 
    
    % Evaluate the function with those parameters
    for i = 1 : K
        fFunctionAta = fFunctionAta + f( a(1), a(2), a(3), a(4), y(i), x(i));
    end
    
    fFunctionAta = (1/K)*fFunctionAta;  
    
end

function[ aGradf ] = gradient_of_f( a, x, y )

    % Variables
    %   K     The number of data points
    %   x     Vector of x data
    %   y     Vector of y data
    
    % Get the length of the data vectors
    K = length( x ); 
    
    % Initialize vector of gradient values
    aGradf = zeros( 4, 1 ); 
    
    % Start filling in those values 
    % Lol useless code
    for ai = 1 : 4
        fSum = 0; 
        
        if ai == 1
            for i = 1 : K
                fSum = fSum + (y(i) - a(1) - a(2)*sin(a(3) * x(i) + a(4)));
            end
            
        elseif ai == 2
            for i = 1 : K
                fSum = fSum + (y(i) - a(1) - a(2)*sin(a(3) * x(i) + a(4))) *...
                    (sin(a(3) * x(i) + a(4)));
            end
            
        elseif ai == 3
            for i = 1 : K
                fSum = fSum + (y(i) - a(1) - a(2)*sin(a(3) * x(i) + a(4))) *...
                    (a(2)*x(i) * cos(a(3)*x(i) + a(4)));
            end
            
        elseif ai == 4
            for i = 1 : K
                fSum = fSum + (y(i) - a(1) - a(2)*sin(a(3) * x(i) + a(4))) *...
                    (a(2)*cos(a(3) * x(i) + a(4)));
            end
                    

        end
        
        aGradf( ai ) = (-2/K) * fSum; 
    end

end

function[ aaHessian ] = hessian_of_f( a, x, y )

    % Variables
    %
    
    % Get the length of the data set
    K = length( x ); 
    
    % Initialize the hessian matrix
    aaHessian = zeros( 4, 4 ); 
    
    % Start filling it in
    for iRow = 1 : 4
        for iCol = 1 : 4
            fSum = 0;
            
            % ROW 1
            if iRow == 1 && iCol == 1
                for i = 1 : K
                    fSum = fSum - 1;
                end
                
            elseif iRow == 1 && iCol == 2 || iRow == 2 && iCol == 1
                for i = 1 : K
                    fSum = fSum - sin(a(3)*x(i) + a(4));
                end
                
            elseif iRow == 1 && iCol == 3 || iRow == 3 && iCol == 1
                for i = 1 : K
                    fSum = fSum - a(2)*x(i)*cos(a(3)*x(i) + a(4));
                end
                
            elseif iRow == 1 && iCol == 4 || iRow == 4 && iCol == 1
                for i = 1 : K
                    fSum = fSum - a(2)*cos(a(3)*x(i) + a(4));
                end
                
            % ROW 2
            elseif iRow == 2 && iCol == 2 
                for i = 1 : K
                    fSum = fSum - (sin(a(3)*x(i) + a(4)))^2;
                end
                
            elseif iRow == 2 && iCol == 3 || iRow == 3 && iCol == 2
                for i = 1 : K
                    fSum = fSum - a(2)*x(i)*sin(2*(a(3)*x(i) + a(4))) + ...
                        x(i)*y(i)*cos(a(3)*x(i) + a(4)) - a(1)*x(i)*cos(a(3)*x(i) + a(4)); 
                end
                
            elseif iRow == 2 && iCol == 4 || iRow == 4 && iCol == 2
                for i = 1 : K
                    fSum = fSum - a(2)*sin(2*(a(3)*x(i) + a(4))) + y(i)*cos(a(3)*x(i) + a(4))...
                        - a(1)*cos(a(3)*x(i) + a(4)); 
                end
                
            % ROW 3
            elseif iRow == 3 && iCol == 3
                for i = 1 : K
                    fSum = fSum + a(2)*x(i)*(-a(2)*x(i)*(cos( a(3)*x(i)+a(4)) )^2 ...
                        - x(i) * sin( a(3)*x(i)+a(4) ) * ...
                        (y(i) - a(1) - a(2)*sin( a(3)*x(i)+a(4) )));
                end
            
            elseif iRow == 3 && iCol == 4 || iRow == 4 && iCol == 3
                for i = 1 : K
                    fSum = fSum + a(2)*x(i) * (-a(2)*(cos( a(3)*x(i) + a(4) ))^2 ...
                        - sin( a(3)*x(i) + a(4) )*(y(i) - a(1) - a(2)*sin(a(3)*x(i) + a(4))));
                end
            
            % ROW 4
            elseif iRow == 4 && iCol == 4
                for i = 1 : K
                    fSum = fSum + a(2)*(-a(2)*(cos( a(3)*x(i) + a(4) ))^2 ...
                        - sin(a(3)*x(i) + a(4)) * (y(i) - a(1) - a(2)*sin(a(3)*x(i) + a(4)))); 
                end
            end

            aaHessian( iRow, iCol ) = (-2/K)*fSum;
            
        end
    end
    
end

function[ aaParams, aError, k] = newtons_method( a0, x, y )
% newtons_method will take in a vector of parameters, the x data points,
% the y data points. It will then run newton's method to find the vector
% *a* that will minimize our objective function. 
% Input
%    a         Vector of parameters; will be updated throughout
%    x         The x value in the pair of data points
%    y         The y value in the pair of data points
% Output
%    aaParams  Matrix containing the parameters at different iterations
%    aError    Vector containing the error at different iterations
%    k         The number of iterations it took to find optimized
%                 parameters

    % Local Variables
    %   epsilon     Threshol of convergence
    %   nMaxIts     Maximum number of iterations
    %   a0          The initial guess for *a*
    %   e0          The initial error, set to 999 as approximate for inifin
    %   k           Iteration number
    %   tk          our initial tkstar guess
    %   tkstar      Step size parameter 
    
    % Assign the threshold of Convergence
    epsilon = 1e-6; 
   
    % Assign maximum number of iterations and initial parameters guess
    nMaxIts = 100; 
    
    % Initialize matrix to hold parameters for each iteration
    aaParams = zeros( 4, nMaxIts ); 
    aaParams( :, 1 ) = a0;
    
    % Assign initial error 
    e0 = 999; 
    
    % Initialize error vector
    aError = zeros( 1, nMaxIts );
    aError( 1 ) = e0;
    
    % Initialize iteration to zero
    k = 0; 
    
    % Run Newton's method with the modifications
    while aError( k+1 ) > epsilon && k < nMaxIts
        
        % Calculate the gradient and the hessian for the current a vector
        aGradf = gradient_of_f( aaParams( :, k+1), x, y );
        aaHessian = hessian_of_f( aaParams( :, k+1), x, y );
        
        % Find r as defined in project
        r = aaHessian \ aGradf; 
        
        % Find a value for tk* for each iteration of NM
        fGamma = 0.95;
        tk = 1;     % 1 for NM, 5 for BFGS
        tkstar = tk; 
        b = 0.5*( aGradf )'*r;
        
        % Approximate the optimal tkstar 
        while func( aaParams( :, k+1 ), x, y ) - ...
                func( aaParams( :, k+1) - tk*r , x, y) < tk*b && tk > 1e-10
            tk = fGamma * tk;
            if func( aaParams( :, k+1) - tk*r , x, y) < ...
                    func( aaParams( :, k+1) - tkstar*r , x, y)
                tkstar = tk; 
            end
           
        end

        % Calculate the next set of parameter values
        aaParams( :, k+2 ) = aaParams( :, k+1) - tkstar * (aaHessian \ aGradf);
        
        % Calculate the error
        aError( 1, k+2 ) = sqrt( (aaParams( 1, k+2 ) - aaParams( 1, k+1 ))^2 +...
            (aaParams( 2, k+2 ) - aaParams( 2, k+1 ))^2 + ...
            (aaParams( 3, k+2 ) - aaParams( 3, k+1 ))^2 + ...
            (aaParams( 4, k+2 ) - aaParams( 4, k+1 ))^2);
        
        
        % Update the counter     
        k = k + 1; 
        
    end
end

function[ aaParams, aError, k ] = BFGS_method( a, x, y ) 
% BFGS_method will take in a vector containing the initial guess of the
% parameters and the data points as its x and y components. 
% Input
%
% Output
%  

    % Local Variables
    %
    
    % Set the threshold for convergence
    epsilon = 1e-6;
    
    % Assign the maximum number of iterations
    nMaxIts = 100; 
    
    % Initialize the matrix to hold parameters of each iteration
    aaParams = zeros( 4, nMaxIts ); 
    aaParams( :, 1) = a;
    
    % Initialize vector to hold the error approximation of each iteration
    aError = zeros( 1, nMaxIts );
    aError( 1 ) = 99; 
    
    % Initialize counter to be 0
    k = 0;
    
    % Initialize D0 to be the identity matrix
    aaDK = eye( 4 ); 
    
    % Start the BFGS method
    while aError( k+1 ) > epsilon && k < nMaxIts
        
        % Find the gradient at the current parameter values
        aGradf = gradient_of_f( aaParams( :, k+1), x, y );
        
        % calculate the p(k) as defined 
        pk = aaDK \ aGradf;
        
        % Find the tkstar that minimizes f(xk - tpk) 
        fGamma = 0.95;
        tk = 5; 
        tkstar = tk; 
        b = 0.5*( aGradf )' * pk; 
        
        
        while func( aaParams( :, k+1 ), x, y ) - ...
                func( aaParams( :, k+1 ) - tk*pk, x, y ) < tk*b && tk > 1e-10
            tk = fGamma * tk; 
            if func( aaParams( :, k+1 ) - tk*pk, x, y ) < ...
                    func( aaParams( :, k+1 ) - tkstar*pk, x, y )
                tkstar = tk;
            end
        end
        
        % Calculate the next set of parameters
        aaParams( :, k+2 ) = aaParams( :, k+1 ) - tkstar * pk;
        
        % Calculate the error 
        aError( k+2 ) = sqrt( ( aaParams( 1, k+2 ) - aaParams( 1, k+1 ) )^2 + ...
            ( aaParams( 2, k+2 ) - aaParams( 2, k+1 ) )^2 + ...
            ( aaParams( 3, k+2 ) - aaParams( 3, k+1 ) )^2 + ...
            ( aaParams( 4, k+2 ) - aaParams( 4, k+1 ) )^2 );
        
        % Calculate the dk value
        adk = aaParams( :, k+2 ) - aaParams( :, k+1 ); 
        
        % Calculate the yk value
        ayk = gradient_of_f( aaParams( :, k+2), x, y ) - ...
            gradient_of_f( aaParams( :, k+1), x, y );
        
        % Calculate the next Dk value
        aaDK = aaDK + ( ayk * ayk' )/( ayk' * adk ) - ...
            ( (aaDK * adk) * (aaDK * adk)' )/( (aaDK * adk)' * adk );
        
        % Update k value
        k = k + 1; 
        
    end

end

function[ aaParams, aError, k ] = DFP_method( a, x, y )
% DFP_method will run the DFP method on a set of parameters and return an
% optimized set of parameter values. 
% Input
%   a           Column vector of parameters that will be treated as an initial
%                   guess
%   x           Vector of x values
%   y           Vector of y values
% Outpu
%   aaParams    Matrix of sequence of parameters ( where each column is a
%                 set of parameters )
%   aError      Vector of error values 
%   k           Iteration that terminates the algorithm 

    % Local variables
    %    epsilon      Threshold for convergence
    %    
    
    % Assign threshold for convergence
    epsilon = 1e-8; 
    
    % Initialize matrix of sequential parameter values
    aaParams = zeros( 4, 100 ); 
    aaParams( :, 1 ) = a; 
    
    % Initialize vector of error values
    aError = zeros( 1, 100 ); 
    aError( 1 ) = 999; 
    
    % Initialize Ek matrix to be identity matrix
    aaEk = eye( 4 ); 
    
    % Initialize iteration value 
    k = 0; 
    
    % Run the DFP Method 
    while aError( k+1 ) > epsilon && k < 100
        
        % Find the gradient at the current parameter values
        aGradf = gradient_of_f( aaParams( :, k+1), x, y );
        
        % Calculate p(k)
        pk = aaEk * aGradf; 
        
        % Find the tkstar value via line search algorithm 
        fGamma = 0.95; 
        tk = 5;
        tkstar = tk; 
        b = 0.5*( aGradf )'*pk; 
        
        while func( aaParams( :, k+1 ), x, y ) - ...
                func( aaParams( :, k+1 ) - tk*pk, x, y ) < tk*b && tk > 1e-10
            tk = fGamma * tk; 
            
            if func( aaParams( :, k+1 ) - tk*pk, x, y ) < ...
                    func( aaParams( :, k+1) - tkstar*pk, x, y )
                tkstar = tk; 
            end
        end
        
        % Calculate the next set of parameters
        aaParams( :, k+2 ) = aaParams( :, k+1 ) - tkstar * pk; 
        
        % Calculate the next set of error values
        aError( k+2 ) = sqrt( (aaParams( 1, k+2 ) - aaParams( 1, k+1 ))^2 + ...
            (aaParams( 2, k+2 ) - aaParams( 2, k+1 ))^2 + ...
            (aaParams( 3, k+2 ) - aaParams( 3, k+1 ))^2 + ...
            (aaParams( 4, k+2 ) - aaParams( 4, k+1 ))^2 );
        
        % Calculate dk 
        dk = aaParams( :, k+2 ) - aaParams( :, k+1 ); 
        
        % Calculate yk 
        ayk = gradient_of_f( aaParams( :, k+2), x, y ) - ...
            gradient_of_f( aaParams( :, k+1), x, y ); 
        
            
        % Calculate Ek+1
        aaEk = aaEk + ( dk * dk' )/( ayk' * dk ) - ...
            ( (aaEk*ayk)*(aaEk*ayk)' )/( (aaEk*ayk)'*ayk );
        
        % Update k           
        k = k+1; 
        
    end
end

function[ ] = plot_iterative_fit( x, y, aaParams, k, fig, achMethod, experiment, color )
% plot_iterative_fit will take in the data points, the parameters, the
% error, and the number of iterations. It will then plot he y values versus
% the x value in open red circle markers. Then, for each iteration using a
% solid green line and no markers, plot the curve generated by evaluating
% h(a, x) using the ak calculated that iteration. 
% Input
%
% Output
%

    % Local Variables
    %

    % Initialize vector of x values
    xvals = [ -5 : 0.01 : 10 ];
    
    % Define anonymous function
    h = @( a, x ) a(1) + a(2) * sin( a(3)*x + a(4) );
    
    % Plot the data points provided
    figure( fig ); 
    clf;
    hold on; 
    plot( x, y, 'ro' )
    
    % For each iterations, plot the curve generated by h( *a*, *x* ) using
    % the a(k) calculated that iteration
    for iteration = 1 : k + 1
        plot( xvals, h( aaParams( :, iteration ), xvals ), color )
    end
    
    % Plot curve generated by evaluating h(a, x) using the final a(k) 
    plot( xvals, h( aaParams( :, k+1 ), xvals ), 'b-' )
    hold off
    
    % Add legend, title, and axes
    legend( 'Data', 'Functions at a^{(k)}')
    xlabel('x')
    ylabel('h(a^{(k)},x)')
    title( [achMethod, ' Experiment ', num2str(experiment),', Iterative Fit Functions'] );
    
end

function[ ] = plot_error_estimates( aError, k, fig, achMethod, experiment )
% plot_error_estimates will take in the vector error values, and the 
% final iteration. Plot the error estimates vs. iteration for k >= 1. 
% Input
%
% Output
%

    % Local Variables
    %

    % Create vector of x values
    aIterations = [ 1 : k ];
    
    % Plot the error estimate values vs iteration
    figure( fig )
    clf; 
    semilogy( aIterations, aError( 2 : k + 1), '-*k' )
    
    % Label axes, legends, and title
    title( [achMethod, ' Experiment ', num2str(experiment),', Error Values vs. Iteration'] );
    xlabel( 'Iteration' )
    ylabel( 'Error Value' )
    
    
end

function[ ] = print_results( aaParams, k, x, y )
% print_results will print the Final a(k), corresponding k, and 
% corresponding f(a, x, y) value. 
% Input
%
% Output
%

    % Local Variables
    %
    
    % Print the results
    fprintf("The final parameter values are: \n")
    for i = 1 : 4
        fprintf(" %0.4f ", aaParams( i, k+1 ) )
    end
    fprintf("\n")
    fprintf("The corresponding k value is: %d\n", k+1)
    fprintf("The final objective function f(a, x, y) value is: %0.5f\n", ...
        func( aaParams( :, k+1 ), x, y ))

    fprintf("\n")
    
    
    
end

function[ ] = run_experiment_4( x, y )
% run_experiment_4 will take in x and y data to be evaluated. 
% Input
%
% Output
%

    % Local Variables
    %  aaRandomParams         Matrix of random parameter values
    %  aOBjectiveFuncVals     Vector that holds the function value evaluated
    %                             at optimized parameter a
    % 

    % Create matrix of 100 randomized parameter values
    aaRandomParams = (pi - 0)*rand( 4, 100 ) + 0;
    
    % Initialize matrix to hold the resulting objective function values
    % Row 1 is for newton's values and Row 2 is for BFGS Method
    aObjectiveFuncVals = zeros( 2, 100 ); 
    
    % Initialize matrices to hold optimized parameter values for each trial
    aaOptimizedNMParams = zeros( 4, 100 ); 
    aaOptimizedBFGSParams = zeros( 4, 100 ); 
    
    % Evaluate the objective function at each random parameter
    for iTrial = 1 : 100
        % Run for Newton's Method
        [ aaParams, aError, k ] = newtons_method( aaRandomParams( :, iTrial), x, y  );
        aObjectiveFuncVals( 1, iTrial ) = func( aaParams( :, k+1 ), x, y );
        
        % Store optimized parameter in appropriate vector
        aaOptimizedNMParams( :, iTrial ) = aaParams( :, k+1 ); 
        
        % Run for BFGS Method
        [ aaParams, aError, k ] = BFGS_method( aaRandomParams( :, iTrial), x, y );
        aObjectiveFuncVals( 2, iTrial ) = func( aaParams( :, k+1 ), x, y ); 
        
        % Store optimized parameter in appropriate vector
        aaOptimizedBFGSParams( :, iTrial ) = aaParams( :, k+1 ); 
        
    end
     
    % Initialize counter for best values
    cGoodFitsNM = 0; 
    cGoodFitsBFGS = 0; 
        
    % Find the best fit in the list of function values (best function value)
    aBestFitNM = aObjectiveFuncVals( 1, 1 );
    aBestaParamNM = aaOptimizedNMParams( :, 1 ); 
    
    aBestFitBFGS = aObjectiveFuncVals( 2, 1 );
    aBestaParamBFGS = aaOptimizedBFGSParams( :, 1 );
    
    
    for i = 1 : 100
        % For newton's method
        if aObjectiveFuncVals( 1, i ) < aBestFitNM
            aBestFitNM = aObjectiveFuncVals( 1, i );
            aBestaParamNM = aaOptimizedNMParams( :, i );
        end
        
        if aObjectiveFuncVals( 1, i ) < 0.1
            cGoodFitsNM = cGoodFitsNM + 1; 
        end
        
        % For BFGS Method
        if aObjectiveFuncVals( 2, i ) < aBestFitBFGS
            aBestFitBFGS = aObjectiveFuncVals( 2, i );
            aBestaParamBFGS = aaOptimizedBFGSParams( :, i ); 
        end
        
        if aObjectiveFuncVals( 2, i ) < 0.1
            cGoodFitsBFGS = cGoodFitsBFGS + 1; 
        end
    end
    
    % Print out results
    fprintf("\nNewton's Method: \n")
    fprintf("Best objective function value over all trials : %d\n", aBestFitNM)
    fprintf("Associated parameters: \n")
    for i = 1 : 4
        fprintf("%3.4f ", aBestaParamNM( i ))
    end
    fprintf("\n")
    fprintf("Percentage of Trials yielding good fits: %d percent\n", ...
        cGoodFitsNM)
    fprintf("\n")
    
    fprintf("BFGS Method: \n")
    fprintf("Best objective function value over all trials: %d\n", aBestFitBFGS)
    fprintf("Associated parameters: \n")
    for i = 1 : 4
        fprintf("%3.4f ", aBestaParamBFGS( i ))
    end
    fprintf( "\n" )
    fprintf("BFGS Method: Percentage of Trials yielding good fits: %d percent\n", ...
        cGoodFitsBFGS)
    
    fprintf("\n")
    % plot_iterative_fit( x, y, aaParams, k, fig, achMethod, experiment )
    % Plot Best Fit Function for Newton's Method
    k = 0; 
    fig = 13; 
    achMethod = 'Newtons Method'; 
    experiment = 4;
    color = '-b';
    plot_iterative_fit( x, y, aBestaParamNM,  k, fig, achMethod, experiment, color ) 
    
    % Plot Objective funciton values for Newton's Method
    figure( 14 ); 
    clf; 
    plot( [ 1 : 100 ], aObjectiveFuncVals( 1, : ), '-ok' )
    xlabel('Trial Number');
    ylabel('Function Value');
    title("Newton's Method, Experiment 4:", "Objective Function Values vs Trial Number" )

    
    % Plot Best Fit Function for BFGS's Method
    k = 0; 
    fig = 15; 
    achMethod = 'BFGS Method';
    plot_iterative_fit( x, y, aBestaParamBFGS, k, fig, achMethod, experiment, color )
    
    % Plot Objective funciton values for BFGS's Method
    figure( 16 ); 
    clf; 
    plot( [ 1 : 100 ], aObjectiveFuncVals( 2, : ), '-ok' )
    xlabel('Trial Number');
    ylabel('Function Value');
    title("BFGSs Method, Experiment 4:", "Objective Function Values vs Trial Number" )
    
end

function[ ] = run_experiment_5( x, y )
% run_exeriment_5 will take in the x and y values. It will then run
% experiment 5 (my extension portion of the project), being the
% david-fletcher-powell (DFP) method for 100 randomized trials (100
% random initial guesses for the vector a). We will consider the objective
% function value f(a, x, y) resulting from each trial and determine the
% best objective function value over all trials and the associated a.
% Additionally, we will determine the percentage of trials which yield good
% fits (where we consider an objective function value less than or equal to
% 0.1 to be a "good fit." Lastly, to display our findings, we will plot the
% best fit function and the objective function values versus the trial
% number. 
% Input
%   x      Vector of x values
%   y      Vector of y values
% Output
%

    % Local Variables
    %  aaRandomParams        Matrix of 100 randomized parameter values
    %                           with each column being a set of
    %                           parameters 
    %  aObjectiveFuncVals    Vector to hold objective function values
    %  fBestObjFuncVal       Best objective function value 
    %  aBestParams           parameter value corresponding to the best
    %                           objective function value 
    
    % Create matrix of 100 randomized parameter values
    aaRandomParams = (pi - 0)*rand( 4, 100 ) + 0;
    
    % Initialize vector to hold objective function values after running
    % each parameter trial through the DFP method
    aObjectiveFuncVals = zeros( 2, 100 ); 
    
    % Initialize matrix to hold optimized parameter values for each trial
    aaOptimizedDFPParams = zeros( 4, 100 ); 
    aaOptimizedBFGSParams = zeros( 4, 100 ); 
    
    % Evaluate the objective function at each random parameter
    for iTrial = 1 : 100
        % Run for DFP Method
        [ aaParams, aError, k ] = DFP_method( aaRandomParams( :, iTrial), x, y  );
        aObjectiveFuncVals( 1, iTrial ) = func( aaParams( :, k+1 ), x, y );
        
        % Store optimized parameter in appropriate vector
        aaOptimizedDFPParams( :, iTrial ) = aaParams( :, k+1 ); 
        
        % Run for BFGS Method
        [ aaParams, aError, k ] = BFGS_method( aaRandomParams( :, iTrial), x, y ); 
        aObjectiveFuncVals( 2, iTrial ) = func( aaParams( :, k+1 ), x, y );
        
        % Store optimized parameter in appropriate vector
        aaOptimizedBFGSParams( :, iTrial ) = aaParams( :, k+1 ); 
    end
    
    % Initialize counter for trials when function value yields good fit
    cnGoodFitsDFP = 0; 
    cnGoodFitsBFGS = 0; 
    
    % Initialize the best objective function value and correspondin params
    fBestObjFuncValDFP = aObjectiveFuncVals( 1, 1 ); 
    aBestParamsDFP = aaOptimizedDFPParams( :, 1 ); 
    
    fBestObjFuncValBFGS = aObjectiveFuncVals( 2, 1 ); 
    aBestParamsBFGS = aaOptimizedBFGSParams( :, 1 ); 
    
    % Find the best objective function value and corresponding paramters
    for i = 1 : 100
        % For DFP method
        if aObjectiveFuncVals( 1, i ) < fBestObjFuncValDFP
            fBestObjFuncValDFP = aObjectiveFuncVals( 1, i );
            aBestParamsDFP = aaOptimizedDFPParams( :, i );
        end
        
        if aObjectiveFuncVals( 1, i ) < 0.1
            cnGoodFitsDFP = cnGoodFitsDFP + 1; 
        end
        
        % For BFGS Method
        if aObjectiveFuncVals( 2, i ) < fBestObjFuncValBFGS
            fBestObjFuncValBFGS = aObjectiveFuncVals( 2, i );
            aBestParamsBFGS = aaOptimizedBFGSParams( :, i ); 
        end
        
        if aObjectiveFuncVals( 2, i ) < 0.1
            cnGoodFitsBFGS = cnGoodFitsBFGS + 1; 
        end
    end
    
    % Print the results for experiment 5
    fprintf("DFP Method: \n")
    fprintf("Best objective funciton value over all trials : %d\n", fBestObjFuncValDFP)
    fprintf("Associated parameters: \n")
    for i = 1 : 4
        fprintf("%3.4f ", aBestParamsDFP( i ))
    end
    fprintf("\n")
    fprintf("Percentage of trials yielding good fits: %d percent\n", ...
        cnGoodFitsDFP)
    fprintf("\n")
    
    fprintf("BFGS Method: \n")
    fprintf("Best objective funciton value over all trials : %d\n", fBestObjFuncValBFGS)
    fprintf("Associated parameters: \n")
    for i = 1 : 4
        fprintf("%3.4f ", aBestParamsBFGS( i ))
    end
    fprintf("\n")
    fprintf("Percentage of trials yielding good fits: %d percent\n", ...
        cnGoodFitsBFGS)
    fprintf("\n")
    
    % plot_iterative_fit( x, y, aaParams, k, fig, achMethod, experiment )
    % Plot Best Fit Function for DFP Method
    k = 0; 
    fig = 17; 
    achMethod = 'DFP Method'; 
    experiment = 5;
    color = 'b'; 
    plot_iterative_fit( x, y, aBestParamsDFP,  k, fig, achMethod, experiment, color) 
    
    % Plot Objective funciton values for Newton's Method
    figure( 18 ); 
    clf; 
    plot( [ 1 : 100 ], aObjectiveFuncVals( 1, : ), '-ok' )
    title("DFP Method, Experiment 5:", "Objective Function Values vs Trial Number" )
    xlabel('Trial Number')
    ylabel('Function Values')
    
    
    % Plot best fit function for BFGS Method
    k = 0; 
    fig = 19; 
    achMethod = 'BFGS Method';
    experiment = 5; 
    color = '-b'; 
    plot_iterative_fit( x, y, aBestParamsBFGS, k, fig, achMethod, experiment, color)
    
    % Plot objective function values for BFGS Method
    figure( 20 );
    clf; 
    plot( [1 : 100], aObjectiveFuncVals( 2, : ), '-ok' )
    title("BFGS Method, Experiment 5:", "Objective Function Values vs Trial Number")
    xlabel('Trial Number')
    ylabel('Function Values')
    
end