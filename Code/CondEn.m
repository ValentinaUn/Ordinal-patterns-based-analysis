function eCE = CondEn( x, delay, order, windowSize )
% CondEn.m efficiently [UK13] calculates values conditional entropy of ordinal patterns [UK14] 
% in sliding windows
% INPUT 
%   - x - time series 
%   - delay - delay between points in ordinal patterns (delay = 1 means successive points)
%   - order - order of ordinal patterns (order+1 - number of points in ordinal patterns)
%   - windowSize - size of a sliding window)
% OUTPUT
%   - eCE - values of conditional entropy of ordinal patterns
%
% REFERENCES
% [UK13] Unakafova, V.A., Keller, K., 2013. Efficiently measuring 
% complexity on the basis of real-world data. Entropy, 15(10), 4392-4415.
% [UK14] Unakafov, A.M. and Keller, K., 2014. Conditional entropy of 
% ordinal patterns. Physica D: Nonlinear Phenomena, 269, pp.94-102.
%
% @author Valentina Unakafova
% @date 18.07.17
% @email UnakafovaValentina(at)gmail.com
    
load( ['table' num2str(order) '.mat'] );% the precomputed table
nPoints      = max( size( x ) );        % the length of the time series
order1       = order + 1;               % for fast computation
orderDelay   = order*delay;
order1Delay  = order1*delay;    
nPatterns    = factorial( order1 );     % the number of ordinal patterns of order opOrder         
patternsDist = zeros( 1, nPatterns );   % the distribution of the ordinal patterns
eCE          = zeros( 1, nPoints );     % values of conditional entropy of ordinal patterns
wordsDist    = zeros( 1, nPatterns*order1 ); % the distribution of (2,d)-words
    
inversionNumbers = zeros( 1, order );       % ordinal pattern inversion numbers (i1,i2,...,id)
prevOP           = zeros( 1, delay );       % previous ordinal patterns for 1:delay
prevWord         = zeros( 1, delay );       % previous (2,d)-words for 1:Delay
patternsInWindow = zeros( 1, windowSize );  % ordinal patterns in the window   
wordWin          = zeros( 1, windowSize );  % (2,d)-words in the window   
peTable( 1:windowSize ) = -( 1:windowSize ).*log( 1:windowSize ); % table of precomputed PE values 
peTable( 2:windowSize ) = ( peTable( 2:windowSize ) - peTable( 1:windowSize-1 ) )./windowSize;
ancNum = nPatterns./factorial( 2:order1 );  % the ancillary numbers    
patternsTable = eval( [ 'table' num2str( order ) ] );
         
for iDelay = 1:delay 
  cnt = iDelay; 
  inversionNumbers( 1 ) = ...
          ( x( orderDelay + iDelay - delay ) >= x( orderDelay + iDelay ) );
  for j = 2:order
    inversionNumbers( j ) = sum( x( ( order - j )*delay + iDelay) >= ...
               x( ( order1 - j )*delay + iDelay:delay:orderDelay + iDelay ) );
  end       
  patternNumber = sum( inversionNumbers.*ancNum ); % first ordinal pattern
    
  % ordinal distribution for current window
  for j = order1Delay + iDelay:delay:windowSize + ( order + 1 )*delay
        word2 = patternNumber*order1; 
        for l = j - orderDelay:delay:j-1             
             if ( x( l ) >= x( j ) )
                  word2 = word2 + 1;
             end
        end      
        patternNumber = patternsTable( word2 + 1 );
        patternsInWindow( cnt ) = patternNumber;
        patternsDist( patternNumber + 1 ) = patternsDist( patternNumber + 1 ) + 1; 
        wordWin( cnt ) = word2;  
        wordsDist( word2 + 1 ) = wordsDist( word2 + 1 ) + 1; 
        cnt = cnt + delay;   
  end       
  prevOP( iDelay )   = patternsInWindow( cnt - delay );
  prevWord( iDelay ) = wordWin( cnt - delay );      
end    
ordDistNorm  = patternsDist/windowSize;
wordDistNorm = wordsDist/windowSize;
ePE = - nansum( ordDistNorm.*log( ordDistNorm ) );  
eCE( windowSize+delay*( order + 1 ) ) = - ePE - nansum( wordDistNorm.*log( wordDistNorm ) );
    
iDelay = mod( windowSize, delay ) + 1; % current shift 1:delay
patternPos = 1;                        % position of the current pattern in the window
% loop over all points
for timePosition = windowSize+delay*( order + 1 ) + 1:nPoints  
  nextPointPos = 0;                    % the position of the next point
  for j = timePosition - orderDelay:delay:timePosition- delay
    if( x( j ) >= x( timePosition ) ) 
        nextPointPos = nextPointPos + 1; 
    end
  end    
  newWord = prevOP( iDelay )*order1 + nextPointPos; % incoming (2,d)-word     
  oldWord = wordWin( patternPos );                  % outcoming (2,d)-word
  prevWord( iDelay ) = newWord;
  wordWin( patternPos ) = newWord;
  
  nNew = patternsTable( newWord + 1 );   % incoming ordinal pattern         
  nOut = patternsInWindow( patternPos ); % outcoming ordinal pattern 
  prevOP( iDelay ) = nNew;  
  patternsInWindow( patternPos ) = nNew;  
  
  newWord = newWord + 1; 
  oldWord = oldWord + 1;  
  nNew    = nNew    + 1; 
  nOut    = nOut    + 1;    
  % update the distribution of (2,d)-words  
  if newWord ~= oldWord                  
     if ( nNew ~= nOut )
        patternsDist( nNew ) = patternsDist( nNew ) + 1; % incoming ordinal pattern
        patternsDist( nOut ) = patternsDist( nOut ) - 1; % outcoming ordinal pattern
        wordsDist( newWord ) = wordsDist( newWord ) + 1; % incoming (2,d)-word
        wordsDist( oldWord ) = wordsDist( oldWord ) - 1; % outcoming (2,d)-word
        eCE( timePosition ) = eCE( timePosition - 1 ) - peTable( patternsDist( nNew ) ) + ...
          peTable( patternsDist( nOut ) + 1 ) + peTable( wordsDist( newWord ) ) - ...
                                      peTable( wordsDist( oldWord ) + 1 );
     else
        wordsDist( newWord ) = wordsDist( newWord ) + 1; % incoming (2,d)-word
        wordsDist( oldWord ) = wordsDist( oldWord ) - 1; % outcoming (2,d)-word
        eCE( timePosition ) = eCE( timePosition - 1 ) + ...
          peTable( wordsDist( newWord ) ) - peTable( wordsDist( oldWord ) + 1 );
     end
  else
     eCE( timePosition ) = eCE( timePosition - 1 );
  end
  iDelay = iDelay + 1; 
  patternPos = patternPos + 1;
  if( iDelay > delay ) 
    iDelay = 1; 
  end
  if( patternPos > windowSize ) 
    patternPos = 1; 
  end
end 
eCE = eCE( windowSize + delay*( order + 1 ):end );   
