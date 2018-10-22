function LE = LEofLogisticMap( rStart, rEnd, rStep )
% @brief calculates Lyapunov exponent of logistic map for r within the interval (3.5,4) 
% x(t+1) = r*x(t)*(1-x(t))
% using derivative for values of control parameter from rStart to rEnd with step rStep
% INPUT
%  - rStart - first value of control parameter r
%  - rEnd   - last values of control parameter r
%  - rStep  - step
% OUTPUT
%  - LE     - values of estimated Lyapunov exponent according to [S03]
% REFERENCES
% [S03] J.C. Sprott, 2003. Chaos and time-series analysis, volume 69. Oxford University Press Oxford.
%
% @author Valentina Unakafova
% @email UnakafovaValentina(at)gmail.com
% @date 11.07.2017

rValues = rStart:rStep:rEnd;  
nPoints = length( rValues );
nIterations = 1000; % number of iterations
LE      = zeros( 1, nPoints );
x       = zeros( 1, nIterations + 1 );
x( 1 ) = 0.1;                
for k = 1:nPoints
  sum = 0;
  for i = 1:nIterations
    x( i + 1 ) = rValues( k )*x( i )*( 1 - x( i ) );
    sum = sum + log( abs( rValues( k ) - 2*rValues( k )*x( i ) ) );   
  end
  LE( k ) = sum / nIterations;     
end