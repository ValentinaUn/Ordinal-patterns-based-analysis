function ax4 = OPplot( x, order, delay, windowSize, forAxis, timeUnits )
% @brief plots the sequence of ordinal patterns for time series x
% @author Valentina Unakafova
% @date 18.07.2017
% @email UnakafovaValentina(at)gmail.com
    
nPoints = max( size( x ) ); % length of time series x
load( ['table' num2str( order ) '.mat'] );
patTable    = eval(['table' num2str(order)]);
dTau        = order*delay;
order1      = order + 1;
ordNumArray = zeros( 1, windowSize ); % ordinal numbers for all tau
nPatterns   = factorial( order1 );      % the number of ordinal patterns of order d    
ancNumArray = nPatterns./factorial( 2:order1 ); % the ancillary numbers for computation ordinal numbers
ordDist     = zeros( 1, nPatterns ); % the distribution of the ordinal patterns
previousPatternAr = zeros( 1, delay );
% for the first window
for iTau = 1:delay  
    counter = iTau; 
    % computation of the first ordinal number
    invNum      = zeros( 1, order );    
    invNum( 1 ) = ( x( dTau + iTau - delay ) >= x( dTau + iTau ) );
    for j = 2:order
        invNum( j ) = sum( x( ( order - j )*delay + iTau ) >= ...
          x( ( order1 - j )*delay + iTau:delay:dTau + iTau ) );
    end        
    ordNumArray( counter ) = sum( invNum.*ancNumArray ); % the first ordinal pattern
    ordDist( ordNumArray( counter ) + 1 ) = ordDist( ordNumArray( counter ) + 1 ) + 1;
    % the ordinal distribution for the first window
    for j = dTau+delay+iTau:delay:windowSize+dTau
        counter = counter+delay;                               
        pos = sum( x( j - dTau:delay:j ) >= x( j ) );
        ordNumArray( counter ) = patTable( ordNumArray( counter - delay )*order1 + pos );
        ordDist( ordNumArray( counter ) + 1 ) = ordDist( ordNumArray( counter ) + 1 ) + 1;            
    end  
    previousPatternAr( iTau ) = ordNumArray( counter );
 end
       
 ordPatPlot( 1:windowSize ) = ordNumArray;% for plot    
 clear ancNumArray invNum pointPos pointPos;
 
 iTau = delay;  % current shift 1:TAU
 iPattern = 1;% position of the current pattern in the window
 % a cycle over the all points in the time series
 x = x( windowSize:nPoints );
 for iPos = 2:nPoints-windowSize-dTau
   % next ordinal number computation                   
   pos = 1;
   for j = iPos:delay:iPos+dTau-delay
       pos = pos + ( x( j ) >= x( iPos + dTau ) );         
   end                  
        
   nNew = patTable( previousPatternAr( iTau )*order1 + pos );% the next ordinal pattern inside the window  
   % permutation entropy for current window computation
   nOut = ordNumArray( iPattern );          
   previousPatternAr( iTau ) = nNew;
   ordNumArray( iPattern )   = nNew; 
   
   nNew = nNew + 1;
   nOut = nOut + 1;       
   if nNew ~= nOut
            % update the ordinal distribution  
      ordDist( nNew ) = ordDist( nNew ) + 1; % new pattern
      ordDist( nOut ) = ordDist( nOut ) - 1; % old pattern
   end
        
   iTau     = iTau + 1;
   iPattern = iPattern + 1;
   if ( iTau > delay ) 
     iTau = 1; 
   end
   if ( iPattern > windowSize ) 
     iPattern = 1; 
   end
   ordPatPlot( iPos + windowSize - 1 ) = nNew; % for plots
 end
    
 % PLOTTING
 plotFontSize   = 8;
 % please, change MarkerSize to 0.1 when saving figure to EPS file
 plotMarkerSize = 0.5;
 ordPatPlot     = ordPatPlot/nPatterns;  
 ordPatPlot     = 1 - ordPatPlot;
 ax4            = subplot( 4, 1, 2 ); 
 set( gca, 'FontSize', plotFontSize ); 
    
 plotLength = length( ordPatPlot );
 plot( forAxis( 1:plotLength ), ordPatPlot( 1:plotLength ), 'ko', ...
                                          'MarkerSize', plotMarkerSize ); 
 grid on;
 axis( [ 0 forAxis( plotLength ) 0 1.0 ] );
 if ( exist( 'timeUnits' ) )
    xlabel( timeUnits, 'FontSize', plotFontSize );
 else
    xlabel( 'Time', 'FontSize', plotFontSize );
 end
 title( '(b) ORDINAL PATTERNS', 'FontSize', plotFontSize );   
