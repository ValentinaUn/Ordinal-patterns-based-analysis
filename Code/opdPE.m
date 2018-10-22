function [ ePE, ordPatPlot, ordDistPlot ] = ...
  opdPE( x, delay, order, orderOPD, windowSize, timeScale, timeUnits )
% OPdistVis.m efficiently computes [UK13] and plots the following ordinal characteristics for time series x:
% sequence of ordinal patterns, ordinal distribution and permutation entropy [BKP02] 
% according to idea from [KM05]
% INPUT 
%   - x          - input time series
%   - delay      - delay between points in ordinal patterns (delay = 1 means successive points)
%   - order      - order of ordinal patterns used for computing ePE and ordinal paatterns distributions
%   - orderSeq   - order of ordinal patterns for visualizing sequence of
%         ordinal patterns (we use two different orders for better visualization)
%   - windowSize - window size
% OPTIONAL INPUT PARAMETERS
%   - timeScale  - optional time axis for xticks
%   - timeUnits  - optional time units for xlabel
% OUTPUT:
%   - ePE - values of permutation entropy
% OPTIONAL OUTPUT PARAMETERS:
%   - ordPatPlot  - sequence of ordinal patterns numbers with time
%   - ordDistPlot - distributions of ordinal patterns with time
%
% REFERENCES
% [BKP02] Bandt C., Pompe B., 2002. Permutation entropy: a natural complexity 
% measure for time series. Physical review letters, APS
% [KM05] Keller, K., and M. Sinn, 2005. Ordinal analysis of time series. 
% Physica A: Statistical Mechanics and its Applications 356.1: 114--120
% [UK13] Unakafova, V.A., Keller, K., 2013. Efficiently Measuring 
% Complexity on the Basis of Real-World Data. Entropy, 15(10), 4392-4415.
%
% @author Valentina Unakafova
% @date 18.07.17
% @email UnakafovaValentina(at)gmail.com
    
 load( ['table' num2str( order ) '.mat'] );
 nPoints   = max( size( x ) ); % length of time series x
 order1    = order+1;          % for fast computation
 dTau      = order*delay;
 nPatterns = factorial( order1 );   % amount of ordinal patterns of order order               
 ordDist   = zeros( 1, nPatterns ); % distribution of the ordinal patterns
 ePE       = zeros( 1, nPoints - windowSize - dTau ); % values of permutation entropy    
                    
 previousPatternAr = zeros( 1, delay );      % last ordinal numbers
 ordNumArray       = zeros( 1, windowSize ); % ordinal numbers for all tau
 ancNumArray       = nPatterns./factorial( 2:order1 ); % ancillary numbers for computation ordinal numbers    
 patternsTable     = eval( ['table' num2str( order )] );
    
 tablePE = zeros( 1, windowSize ); 
 tablePE( 1:windowSize ) = -( 1:windowSize ).*log( 1:windowSize ); % table of values $g(j)$  
 tablePE( 2:windowSize ) = ( tablePE( 2:windowSize ) - tablePE( 1:windowSize - 1 ) )./windowSize;       
    
 % variables for plots
 try
  ordDistPlot = zeros( nPatterns, nPoints - windowSize - dTau );
 catch
   error( 'For plotting reuslts we need factorial(order+1) x nPoints, currently MATLAB does not provide that much memory, try to decrease order' );
 end
 ordPatPlot  = zeros( 1, nPoints - windowSize - dTau );
             
 % for the first window
 for iTau = 1:delay  
     counter = iTau; 
     % computing the first ordinal number
     invNum = zeros( 1, order );    
     invNum( 1 ) = ( x( dTau + iTau - delay ) >= x( dTau + iTau ) );
     for j = 2:order
         invNum( j ) = sum( x( ( order - j )*delay + iTau ) >= ...
           x( ( order1 - j )*delay + iTau:delay:dTau + iTau ) );
     end        
     ordNumArray( counter ) = sum( invNum.*ancNumArray ); % the first ordinal pattern
     ordDist( ordNumArray( counter ) + 1 ) = ordDist( ordNumArray( counter ) + 1 ) + 1;
     % ordinal distribution in the first window
     for j = dTau + delay + iTau:delay:windowSize + dTau
        counter = counter + delay;                               
        pos = sum( x( j-dTau:delay:j ) >= x( j ) );
        ordNumArray( counter ) = patternsTable(ordNumArray( counter - delay )*order1 + pos);
        ordDist( ordNumArray( counter ) + 1 ) = ordDist(ordNumArray( counter ) + 1 ) + 1;            
     end  
     previousPatternAr( iTau ) = ordNumArray( counter );
 end
       
 ordPatPlot( 1:windowSize ) = ordNumArray;% for plot
 ordDistNorm = ordDist/windowSize;        % normalisation of the ordinal distribution
 ordDistPlot( :, 1 ) = ordDist;
 ePE( 1 ) = - nansum( ordDistNorm( 1:nPatterns ).*log( ordDistNorm( 1:nPatterns ) ) );
    
 clear ancNumArray invNum pointPos ordDistNorm pointPos;
    
 iTau = mod( windowSize, delay ) + 1; % current shift 1:TAU
 iPattern = 1;% position of the current pattern in the window
 % a cycle over all points in time series x
 x1 = x( windowSize:nPoints );
 for iPos = 2:nPoints - windowSize - dTau
     % computing the next ordinal number
     pos = 1;
     for j = iPos:delay:iPos + dTau - delay          
       pos = pos + ( x1( j ) >= x1( iPos + dTau ) );         
     end                  
        
     nNew = patternsTable( previousPatternAr( iTau )*order1 + pos );% the next ordinal pattern inside the window
     % permutaion entropy for current window computation
     nOut = ordNumArray( iPattern );          
     previousPatternAr( iTau ) = nNew;
     ordNumArray( iPattern )   = nNew; 

     nNew = nNew + 1;
     nOut = nOut + 1;       
     if nNew ~= nOut
         % update the ordinal distribution  
         ordDist( nNew ) = ordDist( nNew ) + 1; % new pattern
         ordDist( nOut ) = ordDist( nOut ) - 1; % old pattern
         ePE( iPos ) = ePE( iPos - 1 ) + ...
           tablePE( ordDist( nNew ) ) - tablePE( ordDist( nOut ) + 1 );
     else
         ePE( iPos ) = ePE( iPos - 1 );
     end
        
     iTau     = iTau + 1;
     iPattern = iPattern + 1;
     if ( iTau > delay ) 
       iTau = 1; 
     end
     
     if ( iPattern > windowSize ) 
       iPattern = 1;       
     end
        
     % for plots
     ordDistPlot( :, iPos )              = ordDist;
     ordPatPlot( iPos + windowSize - 1 ) = nNew;
 end
    
 % PLOTTING
 ePE = ePE/order;     
      
 figure;
 plotFontSize  = 8;
 plotLineWidth = 0.2;
 set( gca, 'FontSize', plotFontSize );
 forAxisPE  = length( x )  - length( ePE ) + 1:length( x );
 forAxisOP  = length( x )  - length( ordDistPlot( 1, : ) ) + 1:length( x );
 forAxisOP2 = order*delay:length( x );
 if ( exist( 'timeScale' ) )
    forAxisPE  = timeScale( length( timeScale ) - length( forAxisPE ) + 1:end );
    forAxisOP  = timeScale( length( timeScale ) - length( forAxisOP ) + 1:end );
    forAxisOP2 = timeScale( order*delay:end );
 end
 endRec = forAxisPE( end );
 ordDistPlot = ordDistPlot/windowSize;
 ordDistPlot = ( fliplr( ordDistPlot' ) )';
    
 ax1 = subplot( 4, 1, 1 );
 if ( exist( 'timeScale' ) )
   plot( timeScale, x, 'k', 'LineWidth', plotLineWidth ); 
 else
   plot( x, 'k', 'LineWidth', plotLineWidth );
 end
 axis( [ 0, endRec, min( x ), max( x )] );
 set( gca, 'ytick', [] );    
 title( '(a) ORIGINAL TIME SERIES', 'FontSize', plotFontSize );
 if ( exist( 'timeUnits' ) )
   xlabel( timeUnits, 'FontSize', plotFontSize );
 else
   xlabel( 'Time', 'FontSize', plotFontSize );
 end
    
 % preparing ordinal pattern plot                    
 ax2 = subplot( 4, 1, 3 ); 
 set( gca, 'FontSize', plotFontSize); 
 plot( forAxisOP, ordDistPlot( 1, : ), 'k', 'LineWidth', plotLineWidth ); 
 hold on; grid on;
 axis( [ 0, endRec, 0, 1] );
 set( gca, 'ytick', [] );    
 title( '(c) ORDINAL PATTERNS DISTRIBUTION', 'FontSize', plotFontSize );
 if ( exist( 'timeUnits' ) )
   xlabel( timeUnits, 'FontSize', plotFontSize );
 else
   xlabel( 'Time', 'FontSize', plotFontSize );
 end
        
 for i = 2:nPatterns
     ordDistPlot( i, : ) = ordDistPlot( i, : ) + ordDistPlot( i - 1, : );                
     plot( forAxisOP, ordDistPlot( i, 1:end ), 'k', 'LineWidth', plotLineWidth ); 
     hold on;
 end      
 grid on;
                
 ax3 = subplot( 4, 1, 4 ); 
 set( gca, 'FontSize', plotFontSize ); 
 plot( forAxisPE, ePE, 'k'); grid on;
 axis( [ 0, endRec, min( ePE ) - 0.01, max( ePE ) + 0.01 ] );
 if ( exist( 'timeUnits' ) )
   xlabel( timeUnits, 'FontSize', plotFontSize );
 else
   xlabel( 'Time', 'FontSize', plotFontSize );
 end
 title( '(d) PERMUTATION ENTROPY', 'FontSize', plotFontSize );   
    
 % plot of the sequence of ordinal patterns numbers for order order
 if ( exist( 'timeUnits' ) )
   ax4 = OPplot( x, orderOPD, delay, windowSize, forAxisOP2, timeUnits );
 else
   ax4 = OPplot( x, orderOPD, delay, windowSize, forAxisOP2 );
 end
 linkaxes( [ ax1, ax2, ax3, ax4 ], 'x' );
