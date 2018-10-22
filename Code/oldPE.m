function ePE = oldPE( x, delay, order, windowSize )
% @brief oldPE computes values of permutation entropy [BKP02] divided by order in
% sliding windows according to a method from [KES07]
% INPUT
%   - x     - time series
%   - delay - delay between points in ordinal patterns (delay = 1 means successive points)
%   - order - order of ordinal patterns (order+1 - number of points in ordinal patterns)
%   - windowSize - size of a sliding window
% OUTPUT
%   - ePE   - values of permutation entropy
% REFERENCES
% [BKP02] Bandt C., Pompe B., 2002. Permutation entropy: a natural complexity 
% measure for time series. Physical review letters, APS
% [KES07] Keller, K., Emonds, J., Sinn, M., 2007
% Time series from the ordinal viewpoint. Stoch. Dynam., 2, 247-272.
%
% @author Valentina Unakafova
% @date 18.07.17
% @email UnakafovaValentina(at)gmail.com

nPoints   = numel( x );                   % length of the time series
order1    = order + 1; 
dTau      = order*delay; 
nPatterns = factorial( order1 );          % amount of ordinal patterns of order d               
opd       = zeros( 1, nPatterns );        % distribution of ordinal patterns
ePE       = zeros( 1, nPoints );          % values of permutation entropy    
inversionNumbers = zeros( delay, order ); % ordinal pattern $(i_1,i_2,...,i_d)$
ordinalPatternsInWindow = zeros( 1, windowSize ); % ordinal patterns in the window
ancNum = nPatterns./factorial( 2:order1 );  % ancillary numbers           
for iDelay = 1:delay                        % loop for the first window
    cnt = iDelay;
    inversionNumbers( iDelay, 1 ) = ( x( dTau + iDelay - delay ) >= x( dTau + iDelay ) );
    for k = 2:order
        inversionNumbers( iDelay, k ) = sum( x( ( order - k )*delay + iDelay ) >= ...
          x( ( order1 - k )*delay + iDelay:delay:dTau + iDelay ) );
    end   
    % the first ordinal patterns
    ordinalPatternsInWindow( cnt )        = sum( inversionNumbers( iDelay, : ).*ancNum ) + 1;       
    opd( ordinalPatternsInWindow( cnt ) ) = opd(ordinalPatternsInWindow( cnt ) ) + 1;  
    for t = dTau+delay+iDelay:delay:windowSize+dTau   % loop for the next ord. patterns 
        inversionNumbers( iDelay, 2:order ) = inversionNumbers( iDelay, 1:order - 1 );
        inversionNumbers( iDelay, 1 ) = ( x( t - delay ) >= x( t ) );
        for j = 2:order
            if ( x( t - j*delay ) >= x( t ) )
              inversionNumbers( iDelay, j ) = inversionNumbers( iDelay, j ) + 1;
            end
        end        
        ordinalPatternNumber        = sum( inversionNumbers( iDelay, : ).*ancNum ) + 1;
        opd( ordinalPatternNumber ) = opd( ordinalPatternNumber ) + 1;
        cnt = cnt + delay;
        ordinalPatternsInWindow( cnt ) = ordinalPatternNumber;            % the next ordinal pattern
    end  
end  
ordDistNorm = opd/windowSize;
ePE( windowSize + delay*order ) = ...
  -nansum( ordDistNorm( 1:nPatterns ).*log( ordDistNorm( 1:nPatterns ) ) ); 

iDelay = mod( windowSize, delay ) + 1;  % current shift $1:\tau$
iPattern = 1;                           % current pattern in the window
for t = windowSize + dTau + 1:nPoints   % loop for whole time-series
    inversionNumbers( iDelay, 2:order) = inversionNumbers( iDelay, 1:order-1 );
    inversionNumbers( iDelay, 1 ) = ( x( t - delay ) >= x( t ) );
    for j = 2:order
        if ( x( t - j*delay ) >= x( t ) )
          inversionNumbers( iDelay, j ) = inversionNumbers( iDelay, j ) + 1;
        end
    end
    nNew = sum( inversionNumbers( iDelay, : ).*ancNum ) + 1;  % "incoming" ordinal pattern n  
    nOut = ordinalPatternsInWindow( iPattern );               % "outcoming" ordinal pattern          
    ordinalPatternsInWindow( iPattern ) = nNew;
    if ( nNew ~= nOut )                 % update the distribution  
       opd( nNew ) = opd( nNew ) + 1;   % "incoming" ordinal pattern
       opd( nOut ) = opd( nOut ) - 1;   % "outcoming" ordinal pattern
       ordDistNorm = opd/windowSize; 
      ePE( t )  = -nansum( ordDistNorm( 1:nPatterns ).*log( ordDistNorm( 1:nPatterns ) ) );
    else
       ePE( t ) = ePE( t - 1 );
    end  
    iDelay = iDelay + 1;  
    iPattern = iPattern + 1;
    if ( iDelay > delay ) 
      iDelay = 1; 
    end
    if ( iPattern > windowSize ) 
      iPattern = 1; 
    end
end      
ePE = ePE( windowSize + delay*order:end )/order;
