function ePE = PEeq( x, delay, order, windowSize )
% @brief PEeq.m efficiently [UK13] computes permutation entropy [BKP02] 
% for the case of ordinal patterns with tied ranks [UK13,U15]
% in maximally overlapping sliding windows
% INPUT 
%   - x - considered time series
%   - delay - delay between points in ordinal patterns (delay = 1 means successive points)
%   - order - order of the ordinal patterns with tied ranks
%   - windowSize - size of a sliding window in time stamps
% OUTPUT 
%   - ePE - values of permutation entropy
%
% REFERENCES
% [BKP02] Bandt C., Pompe B., 2002. Permutation entropy: a natural complexity 
% measure for time series. Physical review letters, APS
% [UK13] Unakafova, V.A., Keller, K., 2013. Efficiently Measuring 
% Complexity on the Basis of Real-World Data. Entropy, 15(10), 4392-4415.
% [U15] Unakafova, V.A., 2015. Investigating measures of complexity 
% for dynamical systems and for time series (Doctoral dissertation, University of Luebeck).
%
% @author Valentina Unakafova
% @date 11.07.17
% @email UnakafovaValentina(at)gmail.com

load( ['tableEq' num2str( order ) '.mat'] ); % the precomputed table
patternsTable = eval(['tableEq' num2str( order ) ] );% of successive ordinal patterns    
nPoints = numel( x );                        % length of time series
dTau    = order*delay;
nPatterns = 1;                 
for i = 3:2:2*order+1
  nPatterns = nPatterns*i;           
end        
patternsDistribution = zeros( 1, nPatterns ); % distribution of the modified ordinal patterns
ePE = zeros( 1, nPoints );          % values of permutation entropy    
isEquality = zeros( delay, order ); % indicator of equality                    
prevOP = zeros( 1, delay );         % previous modified ordinal patterns for all delays
ordinalPatternsInWindow = zeros( 1, windowSize ); % modified ordinal patterns in the window
ancNum = ones( 1, order );          % ancillary numbers
for j = 2:order
  ancNum( j ) = ancNum( j - 1 )*( 2*j - 1 );
end    
peTable( 1:windowSize ) = -( 1:windowSize ).*log( 1:windowSize ); % table of precomputed values  
peTable( 2:windowSize ) = ( peTable( 2:windowSize ) - peTable( 1:windowSize - 1 ) )./windowSize;        
for iTau = 1:delay         % all delays     
  cnt = iTau;
  inversionNumbers = zeros( 1, order );
  t = dTau + iTau;         % current time of the last point in modified ordinal pattern  
  for j = 1:d              % determining modified ordinal patterns               
    for i = j-1:-1:0
      if ( i == 0 || isEquality( iTau, i ) == 0 )
        if ( x( t - j*delay ) > x( t - i*delay ) )
            inversionNumbers( j ) = inversionNumbers( j ) + 2;
        elseif ( x( t - j*delay ) == x( t - i*delay ) )
            isEquality( iTau, j ) = 1;
        end
      end
    end    
  end        
  inversionNumbers( 1:order ) = inversionNumbers( 1:order ) + isEquality( iTau, 1:order ); % add equality indicator  
  ordinalPatternsInWindow( cnt ) = sum( inversionNumbers.*ancNum );
  patternsDistribution(ordinalPatternsInWindow( cnt ) + 1 ) = ...
    patternsDistribution(ordinalPatternsInWindow( cnt ) + 1 ) + 1; 
  cnt = cnt + delay;
  for t = iTau + delay*( order + 1 ):delay:windowSize + delay*order % loop for the first window
    isEquality( iTau, 2:order ) = isEquality( iTau, 1:order-1 );    % renew equality keeping numbers
    isEquality( iTau, 1 ) = 0;                           
    posL = 1;                         % position of the next point
    eqFlag = 0;                       % indicator of equality
    for i = 1:order                   % determining the position
      if ( isEquality( iTau, i ) == 0 )
        if ( x( t - i*delay ) > x( t ) )
           posL = posL + 2;
        elseif ( x( t ) == x( t - i*delay ) )
           eqFlag = 1;
           isEquality( iTau, i ) = 1;
        end                                    
      end
    end
    posL = posL + eqFlag;  % position of the next point           
    ordinalPatternsInWindow( cnt ) = ...
            patternsTable( ordinalPatternsInWindow( cnt - delay )*( 2*order + 1 ) + posL );
    patternsDistribution( ordinalPatternsInWindow( cnt ) + 1 ) = ...
      patternsDistribution( ordinalPatternsInWindow( cnt ) + 1 ) + 1;            
    cnt = cnt + delay;
  end  
  prevOP( iTau ) = ordinalPatternsInWindow( t - dTau );
end    
OPDnorm = patternsDistribution/windowSize; % normalization of the ordinal distribution
ePE( windowSize + delay*order ) = -nansum( OPDnorm( 1:nPatterns ).*log( OPDnorm( 1:nPatterns ) ) );    

iTau = mod( windowSize, delay )+1;   % current delay
patternPosition = 1;                 % position of the current pattern in the window    
for t = windowSize+delay*order+1:nPoints  % loop for all points in a time series                       
  isEquality( iTau, 2:order ) = isEquality( iTau, 1:order - 1 );
  isEquality( iTau, 1 ) = 0;         
  posL = 1;
  eqFlag = 0;            % is x(j)==x(i)?
  for i = 1:order        % determining the position posL
    if ( isEquality( iTau, i ) == 0 )
      if ( x( t - i*delay ) > x( t ) )
         posL = posL + 2;
      elseif ( x( t ) == x( t - i*delay ) )
         eqFlag = 1;
         isEquality(iTau, i) = 1;
      end
     end
  end
  posL = posL + eqFlag;                                          % position of the next point                                         
  nNew = patternsTable( prevOP( iTau )*( 2*order + 1 ) + posL ); % "incoming" ordinal pattern        
  nOut = ordinalPatternsInWindow( patternPosition );             % "outcoming" ordinal pattern     
  prevOP( iTau ) = nNew;
  ordinalPatternsInWindow( patternPosition ) = nNew; 
  nNew = nNew+1;
  nOut = nOut+1;       
  if nNew ~= nOut                     % if nNew == nOut, ePE does not change
    patternsDistribution( nNew ) = patternsDistribution( nNew ) + 1; % "incoming" ordinal pattern
    patternsDistribution( nOut ) = patternsDistribution( nOut ) - 1; % "outcoming" ordinal pattern
    ePE( t ) = ePE( t - 1 ) + peTable( patternsDistribution( nNew ) ) - ...
                              peTable( patternsDistribution( nOut ) + 1 );
  else
    ePE( t ) = ePE( t - 1 );
  end        
  iTau = iTau + 1; 
  patternPosition  = patternPosition + 1;
  if ( iTau > delay ) 
    iTau = 1; 
  end
  if( patternPosition > windowSize) 
    patternPosition = 1; 
  end
end
ePE = ePE( windowSize + delay*order:end )/order;