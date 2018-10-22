function robustPE = rePE( x, delay, order, windowSize, lowThreshold, highThreshold )
% @brief rePE efficiently computes values of robust permutation entropy [KUU14,U15] 
% in maximally overlapping sliding windows
% INPUT 
%   - x - considered time series 
%   - delay - delay between points in ordinal patterns (delay = 1 means successive points)
%   - order - order of ordinal patterns (order+1 = number of points in ordinal patterns)
%   - windowSize - size of a sliding window
%   - lowThreshold, upperThreshold - lower and upper thresholds, correspondingly (see [KUU14,U15] for more details)
% OUTPUT 
%   - robustPE - values of robust permutation entropy
%
% REFERENCES:
% [KUU14] Keller, K., Unakafov, A.M. and Unakafova, V.A., 2014.
% Ordinal patterns, entropy, and EEG. Entropy, 16(12), pp.6212-6239.
% [U15] Unakafova, V.A., 2015. Investigating measures of complexity
% for dynamical systems and for time series (Doctoral dissertation, University of Luebeck).
%
% @author Valentina Unakafova
% @date 18.07.2017
% @email unakafovavalentina(at)gmail.com

load( ['table' num2str(order) '.mat'] ); % the precomputed table
peTable    = eval( [ 'table' num2str( order ) ] );    
nPoints    = numel( x );                      % time series length
order1     = order + 1; 
orderDelay = order*delay; 
nPatterns  = factorial( order1 );             % amount of ordinal patterns of order d               
patternsDistribution = zeros( 1, nPatterns ); % distribution of ordinal patterns
inversionNumbers = zeros(1, order);           % inversion numbers
ancNum = nPatterns./factorial( 2:order1 );    % ancillary numbers       
MDthr  = ( order + 1 )*order/8;
prevOP = zeros( 1, delay );                % previous ordinal patterns for 1:delay
patternsInWindow = zeros( 1, windowSize ); % ordinal patterns in the window
robustPE = zeros( 1, nPoints );

MD = zeros( 1, nPoints );
for iTau = 1:delay
    MDar1 = zeros( 1, order );
    MDar2 = zeros( 1, order );
    for i  = 1:order
        MDar1( i ) = sum( abs( x( iTau + ( i - 1 )*delay ) - ...
          x( iTau + i*delay:delay:iTau + orderDelay ) ) < lowThreshold );
        MDar2( i ) = sum( abs( x( iTau + ( i - 1 )*delay ) - ...
          x( iTau + i*delay:delay:iTau + orderDelay ) ) > highThreshold );
    end
    MD( iTau ) = sum( MDar1 )+sum( MDar2 );
    MDar1( 1:order - 1 ) = MDar1( 2:order );
    MDar2( 1:order - 1 ) = MDar2( 2:order );
    MDar1( order ) = 0;
    MDar2( order ) = 0;
    for iPos = iTau+delay:delay:nPoints - orderDelay% - delay + 1
        for j = 0:order-1
            MDar1( j + 1 ) = MDar1( j + 1 ) + ...
              ( abs( x( iPos + j*delay ) - x( iPos + orderDelay ) ) < lowThreshold );
            MDar2( j + 1 ) = MDar2( j + 1 ) + ...
              ( abs( x( iPos + j*delay ) - x( iPos + orderDelay ) ) > highThreshold );
        end
        MD( iPos )           = sum( MDar1 ) + sum( MDar2 );
        MDar1( 1:order - 1 ) = MDar1( 2:order );
        MDar1( order )       = 0;
        MDar2( 1:order - 1 ) = MDar2( 2:order );
        MDar2( order )       = 0;
    end
end

for iTau = 1:delay                     % the first sliding window
  cnt = iTau; 
  inversionNumbers( 1 ) = ( x( orderDelay + iTau - delay ) >= x( orderDelay + iTau ) );
  for j = 2:order
    inversionNumbers( j ) = sum( x( ( order - j )*delay + iTau ) >= ...
      x( ( order1 - j )*delay + iTau:delay:orderDelay + iTau ) );
  end        
  patternsInWindow( cnt ) = sum( inversionNumbers .* ancNum );        % the first ordinal pattern
  OPnumber = patternsInWindow( cnt );              
  if ( MD( cnt ) < MDthr )
        patternsDistribution( OPnumber + 1 ) = patternsDistribution( OPnumber + 1 ) + 1;  
  end
  for j = orderDelay + delay + iTau:delay:windowSize + orderDelay  % loop for the first window
    cnt = cnt + delay;                               
    posL = 1;                        % the position $l$ of the next point
    for i = j - orderDelay:delay:j - delay
        if ( x( i ) >= x( j ) ) 
          posL = posL + 1; 
        end
    end  
    patternsInWindow( cnt ) = peTable( patternsInWindow( cnt - delay )*order1 + posL );
    OPnumber = patternsInWindow( cnt );
    if ( MD( cnt ) < MDthr )
        patternsDistribution( OPnumber + 1 ) = patternsDistribution( OPnumber + 1 ) + 1;            
    end
  end 
  prevOP( iTau ) = patternsInWindow( cnt );
end    
ordDistNorm = patternsDistribution/sum( patternsDistribution );
robustPE( windowSize+delay*order ) = ...
  - nansum( ordDistNorm( 1:nPatterns ).*log( ordDistNorm( 1:nPatterns ) ) );       

iTau = mod( windowSize, delay ) + 1;    % current shift 1:delay
iPat = 1;                               % position of the current pattern in the window
for tPos = windowSize + delay*order + 1:nPoints % loop over all points
  posL = 1;                             % position of the next point
  for j = tPos - orderDelay:delay:tPos - delay
    if( x( j ) >= x( tPos ) ) 
        posL = posL+1; 
    end
  end                          
  nNew = peTable( prevOP( iTau )*order1 + posL ); % "incoming" ordinal pattern         
  nOut = patternsInWindow( iPat );                % "outcoming" ordinal pattern 
  prevOP( iTau ) = nNew;
  patternsInWindow( iPat ) = nNew; 
  nNew = nNew + 1;
  nOut = nOut + 1;       
  % update the distribution of ordinal patterns 
  if ( MD( tPos - orderDelay ) < MDthr )   
     patternsDistribution( nNew ) = patternsDistribution( nNew ) + 1; % "incoming" ordinal pattern
  end
  if ( MD( tPos - windowSize - orderDelay ) < MDthr )
     patternsDistribution( nOut ) = patternsDistribution( nOut ) - 1; % "outcoming" ordinal pattern
  end
  if ( sum( patternsDistribution ) > 0 ) 
    ordDistNorm = patternsDistribution/sum( patternsDistribution );
  else
    ordDistNorm = zeros( 1, length( ordDistNorm ) );
  end
  robustPE( tPos ) = - nansum( ordDistNorm( 1:nPatterns ).*log( ordDistNorm( 1:nPatterns ) ) );
  
  iTau = iTau + 1; 
  iPat = iPat + 1;
  if ( iTau > delay ) 
    iTau = 1; 
  end
  if ( iPat > windowSize ) 
    iPat = 1; 
  end
end 
robustPE = robustPE( windowSize + delay*order:end )/order;
