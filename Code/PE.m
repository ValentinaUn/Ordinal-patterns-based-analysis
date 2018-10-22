function ePE = PE( x, delay, order, windowSize )
% @brief PE efficiently [UK13] computes values of permutation entropy [BKP02] in 
% maximally overlapping sliding windows
%
% INPUT 
%   - x - considered time series 
%   - delay - delay between points in ordinal patterns (delay = 1 means successive points)
%   - order - order of the ordinal patterns (order+1 - number of points in ordinal patterns)
%   - windowSize - size of sliding window
% OUTPUT 
%   - ePE - values of permutation entropy
%
% REFERENCES
% [BKP02] Bandt C., Pompe B., 2002. Permutation entropy: a natural complexity 
% measure for time series. Physical review letters, APS
% [UK13] Unakafova, V.A., Keller, K., 2013. Efficiently measuring complexity 
% on the basis of real-world data. Entropy, 15(10), 4392-4415.
%
% @author Valentina Unakafova
% @date 11.07.2017
% @email UnakafovaValentina(at)gmail.com

 load( ['table' num2str( order ) '.mat'] ); % the precomputed table
 patternsTable = eval( ['table' num2str( order )] );    
 nPoints    = numel( x );            % length of the time series
 opOrder1   = order + 1; 
 orderDelay = order*delay; 
 nPatterns  = factorial( opOrder1 );   % amount of ordinal patterns of order d               
 patternsDistribution = zeros( 1, nPatterns ); % distribution of ordinal patterns
 ePE = zeros( 1, nPoints );            % permutation entropy    
 inversionNumbers = zeros( 1, order ); % inversion numbers of ordinal pattern (i1,i2,...,id)
 prevOP = zeros( 1, delay );           % previous ordinal patterns for 1:opDelay
 ordinalPatternsInWindow = zeros( 1, windowSize ); % ordinal patterns in the window
 ancNum = nPatterns./factorial( 2:opOrder1 );  % ancillary numbers       
 peTable( 1:windowSize ) = -( 1:windowSize ).*log( 1:windowSize  ); % table of precomputed values  
 peTable( 2:windowSize ) = ( peTable( 2:windowSize ) - peTable( 1:windowSize - 1 ) )./windowSize;       
 for iTau = 1:delay 
     cnt = iTau; 
     inversionNumbers( 1 ) = ( x( orderDelay + iTau - delay ) >= x( orderDelay + iTau ) );
     for j = 2:order
         inversionNumbers( j ) = sum( x( ( order - j )*delay + iTau ) >= ...
           x( ( opOrder1 - j )*delay + iTau:delay:orderDelay + iTau ) );
     end        
     ordinalPatternsInWindow( cnt ) = sum( inversionNumbers.*ancNum ); % first ordinal patterns
     patternsDistribution( ordinalPatternsInWindow( cnt )+ 1 ) = ...
       patternsDistribution( ordinalPatternsInWindow( cnt ) + 1 ) + 1;  
     for j = orderDelay + delay + iTau:delay:windowSize + orderDelay % loop for the first window
         cnt = cnt + delay;                               
         posL = 1; % the position of the next point
         for i = j - orderDelay:delay:j - delay
             if( x( i ) >= x( j ) ) 
                 posL = posL + 1; 
             end
         end  
         ordinalPatternsInWindow( cnt ) = ...
           patternsTable( ordinalPatternsInWindow( cnt - delay )*opOrder1 + posL );
         patternsDistribution( ordinalPatternsInWindow( cnt ) + 1 ) = ...
           patternsDistribution( ordinalPatternsInWindow( cnt ) + 1 ) + 1;            
     end  
     prevOP( iTau ) = ordinalPatternsInWindow( cnt );
 end    
 ordDistNorm = patternsDistribution/windowSize;
 ePE( windowSize + delay*order ) = ...
   -nansum( ordDistNorm( 1:nPatterns ).*log( ordDistNorm( 1:nPatterns ) ) );       

 iTau = mod( windowSize, delay ) + 1;   % current shift $1:\tau$
 patternPosition = 1;                   % position of the current pattern in the window
 for t = windowSize + delay*order + 1:nPoints   % loop over all points
     posL = 1;                 % the position of the next point
     for j = t-orderDelay:delay:t-delay
         if( x( j ) >= x( t ) ) 
             posL = posL + 1; 
         end
     end                          
     nNew = patternsTable( prevOP( iTau )*opOrder1 + posL );  % "incoming" ordinal pattern         
     nOut = ordinalPatternsInWindow( patternPosition ); % "outcoming" ordinal pattern 
     prevOP( iTau ) = nNew;
     ordinalPatternsInWindow( patternPosition ) = nNew; 
     nNew = nNew + 1;
     nOut = nOut + 1;       
     if ( nNew ~= nOut )           % update the distribution of ordinal patterns    
         patternsDistribution( nNew ) = patternsDistribution( nNew ) + 1; % "incoming" ordinal pattern
         patternsDistribution( nOut ) = patternsDistribution( nOut ) - 1; % "outcoming" ordinal pattern
         ePE( t ) = ePE( t - 1 ) + ( peTable( patternsDistribution( nNew ) ) - ...
                                        peTable( patternsDistribution( nOut ) + 1  ) );
     else
         ePE( t ) = ePE( t - 1 );
     end      
     iTau = iTau + 1; 
     patternPosition = patternPosition + 1;
     if ( iTau > delay ) 
       iTau = 1; 
     end
     if ( patternPosition > windowSize ) 
       patternPosition = 1; 
     end
 end 
 ePE = ePE( windowSize + delay*order:end )/order;