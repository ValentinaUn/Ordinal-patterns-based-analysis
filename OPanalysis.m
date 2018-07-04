function outdata = OPanalysis( cfg, indata )
% OPanalysis( cfg, indata ) efficiently computes from indata the
% following ordinal-patterns-based measures that are successfully used for non-linear 
% analysis of multivariate time series (see [AKU15,BKP02,KUU14] and ref. there for more details):
% - permutation entropy (cfg.method = 'PE') [BKP02]
% - permutation entropy for ordinal patterns with tied ranks (cfg.method = 'eqPE') [KUU14]
% - permutation entropy and ordinal patterns distributions (cfg.method = 'opdPE')
% - conditional entropy of ordinal patterns (cfg.method = cePE') [UK13]
% - robust permutation entropy (cfg.method = 'rePE') [KUU14,U15]
% - new ordinal-patterns-based measures are to be added!
%
% INPUT (see [KUU14, RMW13, U15] for discussion of parameters choice):
%   - indata    - data to be analyzed (nChannels x nPoints)
%   - cfg.delay - delay between points in ordinal patterns (often delay = 1 is recommended)
%   - cfg.order - order of ordinal patterns (usually order = 3..7 is recommended [AZS2008])
%   - cfg.windowSize - size of sliding in data windows in samples in which measures values
%       are computed (the assumption is that data are stationary in windows
%       [BKP02,KUU14,UK15]). This parameter is defined as the number of ordinal
%       patterns in the sliding window (see [KUU14] for details)
% OPTIONAL INPUT (these values may be omitted)
%   - cfg.lowerThreshold - lower threshold for robust permutation entropy,
%   applicable only for 'rePE' method, (see [KUU14,U15] for details)
%   - cfg.upperThreshold - upper threshold for robust permutation entropy,
%   applicable only for 'rePE' method, (see [KUU14,U15] for details)
%   - cfg.orderSeq - order of ordinal patterns required for 'opdPE' method 
%                             for visualizing sequence of ordinal patterns
%   - cfg.plot  - whether to plot the resulting values (0 - 'no', 1 - 'yes'), '1' by default
%   - cfg.time  - correct time axis corresponding to indata (used for plotting)
%   - cfg.units - units of time passed in cfg.time (string)
% OUTPUT:
%   - outdata - values of a complexity measure computed in successive (maximally overlapping)
%               sliding windows (nChannels x nPoints - windowSize)
%
% EXAMPLE of use (see more examples in examples.m):
%
% cfg            = [];
% cfg.method     = 'opdPE'; % try also 'PE', 'CE', 'rePE' and 'all' here
% cfg.order      = 3;       % for ordinal pattens of order 3 (4-points ordinal patterns)
% cfg.delay      = 1;       % for delay 1 between points in ordinal patterns (successive points)
% cfg.windowSize = 512;     % for window size 512 in time points
% indata         = rand( 1, 7777 );
% for i = 4000:7000         % change of data complexity
%    indata( i ) = 4*indata( i - 1 )*( 1 - indata( i - 1 ) );
% end 
% outdata        = OPanalysis( cfg, indata );
%
% REFERENCES (alphabetical order):
% [AKU15] Amigo, J.M., Keller, K. and Unakafova, V.A., 2015. On entropy, 
% entropy-like quantities, and applications. Discrete & Continuous Dynamical 
% Systems-Series B, 20(10).
% [AZS08] Amigo, J.M., Zambrano, S. and Sanjuan, M.A., 2008. 
% Combinatorial detection of determinism in noisy time series. 
% EPL (Europhysics Letters), 83(6), p.60005.
% [BKP02] Bandt C., Pompe B., 2002. Permutation entropy: a natural complexity 
% measure for time series. Physical review letters, APS
% [BQM2012] Bian, C., Qin, C., Ma, Q.D. and Shen, Q., 2012. Modified 
% permutation-entropy analysis of heartbeat dynamics. Physical Review E, 85(2), p.021906.
% [KM05] Keller, K., and M. Sinn, 2005. Ordinal analysis of time series. 
% Physica A: Statistical Mechanics and its Applications 356.1: 114--120
% [KUU14] Keller, K., Unakafov, A.M. and Unakafova, V.A., 2014. 
% Ordinal patterns, entropy, and EEG. Entropy, 16(12), pp.6212-6239.
% [RMW13] Riedl, M., Muller, A. and Wessel, N., 2013. Practical considerations 
% of permutation entropy. The European Physical Journal Special Topics, 222(2), pp.249-262.
% [UK13] Unakafova, V.A., Keller, K., 2013. Efficiently measuring 
% complexity on the basis of real-world Data. Entropy, 15(10), 4392-4415.
% [U15] Unakafova, V.A., 2015. Investigating measures of complexity 
% for dynamical systems and for time series (Doctoral dissertation, 
% University of Luebeck).
% 
% @author Valentina Unakafova
% @date 18.07.2017
% @email UnakafovaValentina(at)gmail.com

% add paths
addpath( 'Code' );
addpath( 'Data' );
addpath( 'Tables' );

% check dimensionality and type of data
try
  if ( ~isnumeric( indata ) )
    error( 'The data are not numeric' );
  end
  nChannels = size( indata, 1 );
  nPoints   = size( indata, 2 ); 
catch
  error( 'Something is wrong with either dimensions or type of indata' );
end

if ( nChannels > nPoints )
  error('Check dimensionality of indata, it is supposed to be in format nChannels x nPoints');
end

% check delay and order
if ( isfield( cfg, 'order' ) )
  if ( isnumeric( cfg.order ) && mod( cfg.order, 1 ) == 0 )
  else
    error( 'cfg.order must be integer and numeric' );
  end
  
  if ( cfg.order >= 1 && cfg.order < 9 )
    if ( nPoints < 5*factorial( cfg.order ) )
      warning( 'It is recommended to use length( indata ) > 5*factorial( order ). Otherwise the results may be biased (see [AZS08,KUU14] for more details)' );
    end
  else
    if ( cfg.order >= 9 )
      if ( isfield( cfg, 'method' ) && strcmp( cfg.method, 'opdPE' ) )
        warning( 'Orders d>9 are not supported for this method, cfg.order is set to 8' );
        cfg.order = 8;
      else
        warning( 'Orders d>9 take more time to be computed since not efficient version of algorithm is used for this case. Usually orders d = 3...7 are recommended' );
        if ( cfg.order < 11 )
          if( isfield( cfg, 'method' ) && strcmp( cfg.method, 'PE' ) )
            cfg.method = 'oldPE';
          else
            if ( isfield( cfg, 'method' ) )
              warning( ['Orders d>9 are not supported for the method ' cfg.method ', order is set to 8'] );
              cfg.order = 8;
            end
          end
        else
          error( 'This version does not support orders high than 10 due to memory limitations...' );
        end
      end
    else
      error( 'cfg.order is wrong' );
    end
  end
else
  warning( 'The cfg.order field was not specified, set to 3' );
  cfg.order = 3;
end % if isfield( cfg, 'order' )

% check orderSeq
if ( strcmp( cfg.method, 'opdPE' ) == 1 )
  if ( isfield( cfg, 'orderSeq' ) )
    if ( isnumeric( cfg.orderSeq ) && mod( cfg.orderSeq, 1 ) == 0 )
    else
      error( 'cfg.orderSeq must be integer and numeric' );
    end
  
    if ( cfg.orderSeq >= 1 && cfg.orderSeq < 9 )
      if ( nPoints < 5*factorial( cfg.orderSeq ) )
        warning( 'It is recommended to use length( indata ) > 5*factorial( orderSeq ). Otherwise the results may be biased (see [AZS08,KUU14] for more details)' );
      end
    else
      if ( cfg.orderSeq >= 9 )
          warning( 'Orders d>9 are not supported for opdPE method, cfg.orderSeq is set to 8' );
          cfg.orderSeq = 8;
      else
      end
    end
  else % if( isfield( cfg, 'orderSeq' ) )
    warning( 'The cfg.orderSeq field was not specified, set to 3' );
    cfg.orderSeq = 3;
  end  % if isfield( cfg, 'orderSeq' )
end

if ( isfield( cfg, 'delay' ) )
  if ( isnumeric( cfg.delay ) && mod( cfg.delay, 1 ) == 0 )
  else
    error( 'cfg.delay must be integer and numeric' );
  end
  if ( cfg.delay >= 1 && cfg.delay < 10 )
  else
  end
else
  warning( 'The cfg.delay field was not specified, set to 1' );
  cfg.delay = 1;
end

% check window size taking into account order, orederSeq and delay
if ( isfield( cfg, 'windowSize' ) )
  if ( isnumeric( cfg.windowSize ) && mod( cfg.windowSize, 1 ) == 0 )
  else
    error( 'cfg.windowSize must be integer and numeric' );
  end
  
  if ( strcmp( cfg.method, 'CE' ) == 1 || strcmp( cfg.method, 'all' ) == 1 )
    if ( cfg.windowSize + ( cfg.order + 1 )*cfg.delay > nPoints )
      cfg.windowSize = nPoints - ( cfg.order + 1 )*cfg.delay;
      outdata = zeros( nChannels, 1 );
      warning('The cfg.windowSize was too large, it is set now as one window for whole time series');
    else
      outdata = zeros( nChannels, nPoints - cfg.windowSize - ( cfg.order + 1 )*cfg.delay + 1 );
    end
  else
    if ( cfg.windowSize + cfg.order*cfg.delay > nPoints )
       cfg.windowSize = nPoints - cfg.order*cfg.delay;
       outdata = zeros( nChannels, 1 );
       warning('The cfg.windowSize was too large, it is set now as one window for whole time series');
    else
      if ( strcmp( cfg.method, 'opdPE' ) == 1 )
        % check correctness of windowSize for both order and orderSeq
        if ( isfield( cfg, 'orderSeq' ) )
        else
          warning( 'cfg.orderSeq was not specified, set to cfg.order' );
          cfg.orderSeq = cfg.order;
        end
        maxOrder = max( cfg.order, cfg.orderSeq );
        windowSize = nPoints - cfg.windowSize - ( maxOrder + 1 )*cfg.delay;
        if ( windowSize > 1 )
        else
          error( 'Too short time series for the set parameters (check cfg.order (cfg.orderSeq for opdPE method), cfg.delay and cfg.windowSize)' );
        end
    
        outdata = zeros( nChannels, nPoints - cfg.windowSize - ( cfg.order + 1 )*cfg.delay + 1 );
      else
        outdata = zeros( nChannels, nPoints - cfg.windowSize - cfg.order*cfg.delay + 1 );
      end
    end
  end
  % check if windowSize is correct
  if ( cfg.windowSize > 1 )
  else
    error( 'Too short time series for the set parameters (check cfg.order (cfg.orderSeq for opdPE method), cfg.delay and cfg.windowSize)' );
  end
  
else % if ( isfield( cfg, 'windowSize' ) )
  % size of the windows must be different for different methods
  warning( ['Size of the window was not set, ' cfg.method ' is computed for the whole time series'] );
  if ( strcmp( cfg.method, 'all' ) == 1 || strcmp( cfg.method, 'CE' ) == 1 )
    cfg.windowSize = nPoints - ( cfg.order + 1 )*cfg.delay;
  else
    if ( strcmp( cfg.method, 'opdPE' ) == 1 )
      % check correctness of windowSize for both order and orderSeq
      if ( isfield( cfg, 'orderSeq' ) )
      else
        warning( 'cfg.orderSeq was not specified, set to cfg.order' );
        cfg.orderSeq = cfg.order;
      end
      maxOrder = max( cfg.order, cfg.orderSeq );
      windowSize = nPoints - cfg.windowSize - ( maxOrder + 1 )*cfg.delay;
      if ( windowSize > 1 )
      else
        error( 'Too short time series for the set parameters (check cfg.order (cfg.orderSeq for opdPE method), cfg.delay and cfg.windowSize)' );
      end
    end
    cfg.windowSize = nPoints - cfg.order*cfg.delay;
  end
  % check whether window size is correct
  if ( cfg.windowSize > 1 )
  else
    error( 'Too short time series for the set parameters (check cfg.order (cfg.orderSeq for opdPE method), cfg.delay and cfg.windowSize)' );
  end
end

if ( isfield( cfg, 'time' ) )
  if ( isnumeric( cfg.time ) )
      if ( numel( cfg.time ) == length( indata ) )
          if ( all( diff( cfg.time ) > 0 ) )
          else
            warning( 'cfg.time values were not monotonically increasing and cfg.time was removed' );
            cfg = rmfield( cfg, 'time' );
          end
      else
        warning( 'cfg.time is not taken into account since its dimension is wrong' );
        cfg = rmfield( cfg, 'time' );
      end
  else
    warning( 'cfg.time values are not numeric and removed' );
    cfg = rmfield( cfg, 'time' );
  end
end

% check what method is to be used
if ( isfield( cfg, 'method' ) )
  if ( strcmp( cfg.method, 'PE' ) )
    figureTitle = 'Values of permutation entropy';
    fn = @PE;
  else
    if ( strcmp( cfg.method, 'PEeq' ) )
      fn = @PE;
      figureTitle = 'Values of permutation entropy for ordinal patterns with tied ranks';
    else
      if ( strcmp( cfg.method, 'CE' ) )
        fn = @CondEn;
        figureTitle = 'Values of conditional entropy of ordinal patterns';
      else
        if ( strcmp( cfg.method, 'oldPE' ) )
          fn = @oldPE;
          figureTitle = 'Values of permutation entropy';
        else
          if ( strcmp( cfg.method, 'rePE' ) )
            fn = @rePE;
            figureTitle = 'Values of robust permutation entropy';
          else
            if ( strcmp( cfg.method, 'opdPE' ) )
              fn = @opdPE;
            else
               if ( strcmp( cfg.method, 'all' ) )
               else
                  error( ['OPA toolbox does not support the method ' cfg.method ', please, see the supported methods in help'] );
               end
            end
          end
        end
      end
    end
  end
else
  warning( 'The cfg.method field was not specified, the cfg.method is set to permutation entropy' );
  cfg.method = 'PE';
  fn = @PE;
end

if ( isfield( cfg, 'plot' ) )
else
  cfg.plot = 1;
end

if ( strcmp( cfg.method, 'rePE' ) )
  if ( isfield( cfg, 'lowerThreshold' ) && isfield( cfg, 'upperThreshold' ) )
  else
    warning('You have not specified lower and upper thresholds for robust permutation entropy, they are set by default as min and max values of indata' );
    cfg.lowerThreshold = min( indata( 1, : ) );
    cfg.lowerThreshold = max( indata( 1, : ) );
  end
end

for iChannel = 1:nChannels
  if ( strcmp( cfg.method, 'all' ) == 1 && cfg.plot == 1 )
    clear outdata;
    outdataPE = PE( indata( iChannel, : ), cfg.delay, cfg.order, cfg.windowSize );
    % check window size
    outdataCE = CondEn( indata( iChannel, : ), cfg.delay, cfg.order, cfg.windowSize );
    if ( ~isfield( cfg, 'lowerThreshold' ) || ~isfield( cfg, 'upperThreshold' ) )
      cfg.lowerThreshold  = min( min( indata ) );
      cfg.upperThreshold = max( max( indata ) );
    end
    outdataRE = rePE( indata( iChannel, : ), cfg.delay, cfg.order, ...
          cfg.windowSize, cfg.lowerThreshold, cfg.upperThreshold );
        
    if ( isfield( cfg, 'time' ) )
      if ( isnumeric( cfg.time ) && numel( cfg.time ) == length( indata ) )
      else
        cfg = rmfield( cfg, 'time' );
        warning( 'cfg.time is incorrect and removed' );
      end
    end
    
    outdata.PE = outdataPE;
    outdata.RE = outdataRE;
    outdata.CE = outdataCE;
    
    if ( cfg.plot == 1 )
      figure;
      ax1 = subplot( 4, 1, 1 );
      if ( isfield( cfg, 'time' ) && numel( cfg.time ) == length( indata ) )
        plot( cfg.time( end - length( outdataPE ) + 1:end ), ...
              indata( iChannel, end - length( outdataPE ) + 1:end ), ...
              'k', 'LineWidth', 0.2 ); grid on;
      else
        forAxis = 1:length( indata( iChannel, : ) );
        plot( forAxis, indata( iChannel, : ), 'k', 'LineWidth', 0.2 ); grid on;
      end
    
      if ( isfield( cfg, 'units' ) )
        xlabel( cfg.units );
      end
      
      if ( nChannels > 1 )
        title( ['Original time series, channel ' int2str( iChannel ) ] );
      else
        title( 'Original time series' );
      end
    
      ax2 = subplot( 4, 1, 2 );
      if ( isfield( cfg, 'time' ) && numel( cfg.time ) == length( indata ) && isnumeric( cfg.time ) )
        plot( cfg.time( end - length( outdataPE ) + 1:end ), outdataPE, 'k', 'LineWidth', 0.2 ); grid on
      else
        plot( forAxis( end - length( outdataPE ) + 1:end ), outdataPE, 'k', 'LineWidth', 0.2 ); grid on;
      end
      if ( isfield( cfg, 'units' ) )
        xlabel( cfg.units );
      end
      if ( nChannels > 1 )
        title( ['Values of permutation entropy, channel ' int2str( iChannel ) ] );
      else
        title( 'Values of permutation entropy' );
      end
    
      ax3 = subplot( 4, 1, 3 );
      if ( isfield( cfg, 'time' ) && numel( cfg.time ) == length( indata ) && all( diff( cfg.time ) > 0 ) )
        plot( cfg.time( end - length( outdataCE ) + 1:end ), outdataCE, 'k', 'LineWidth', 0.2 ); grid on
      else
        plot( forAxis( end - length( outdataCE ) + 1:end ), outdataCE, 'k', 'LineWidth', 0.2 ); grid on;
      end
      
      if ( isfield( cfg, 'units' ) )
        xlabel( cfg.units );
      end
      
      if ( nChannels > 1 )
        title( ['Values of conditional entropy of ordinal patterns, channel ' int2str( iChannel ) ] );
      else
        title( 'Values of conditional entropy of ordinal patterns' );
      end
    
      ax4 = subplot( 4, 1, 4 );
      if ( isfield( cfg, 'time' ) && numel( cfg.time ) == length( indata ) )
        plot( cfg.time( end - length( outdataRE ) + 1:end ), outdataRE, 'k', 'LineWidth', 0.2 ); grid on
      else
        plot( forAxis( end - length( outdataRE ) + 1:end ), outdataRE, 'k', 'LineWidth', 0.2 ); grid on;
      end
      
      if ( isfield( cfg, 'units' ) )
        xlabel( cfg.units );
      end
      
      if ( nChannels > 1 )
        title( ['Values of robust permutation entropy, channel ' int2str( iChannel ) ] );
      else
        title( 'Values of robust permutation entropy' );
      end
      linkaxes( [ ax1, ax2, ax3, ax4 ], 'x' );
    end % if ( cfg.plot == 1 )
    
  else % if( strcmp( cfg.method, 'all' ) == 1 && cfg.plot == 1 )
    if( strcmp( cfg.method, 'rePE' ) == 1 )
      outdata( iChannel, : ) = fn( indata( iChannel, : ), cfg.delay, cfg.order, ...
          cfg.windowSize, cfg.lowerThreshold, cfg.upperThreshold );
    else
      if ( strcmp( cfg.method, 'opdPE' ) == 1 )
        if ( isfield( cfg, 'time' ) )
          if ( isfield( cfg, 'units' ) )
            outdata( iChannel, : ) = fn( indata( iChannel, : ), cfg.delay, cfg.order, cfg.orderSeq, cfg.windowSize, cfg.time, cfg.units );
          else
            outdata( iChannel, : ) = fn( indata( iChannel, : ), cfg.delay, cfg.order, cfg.orderSeq, cfg.windowSize, cfg.time );
          end % if ( isfield( cfg, 'units' ) )
        else
          if ( ~isfield( cfg, 'orderSeq' ) )
            cfg.orderSeq = cfg.order;
          end
          outdata( iChannel, : )   = fn( indata( iChannel, : ), cfg.delay, cfg.order, cfg.orderSeq, cfg.windowSize );
        end   % if ( isfield( cfg, 'time' ) )
      else
        outdata( iChannel, : )     = fn( indata( iChannel, : ), cfg.delay, cfg.order, cfg.windowSize );
      end     % if ( strcmp( cfg.method, 'opdPE' ) == 1 )
    end       % if( strcmp( cfg.method, 'rePE' ) == 1 )
  end         % if ( strcmp( cfg.method, 'all' ) == 1 && cfg.plot == 1 )
end

% now we plot the complexity values for all channels in separate figures
if ( strcmp( cfg.method, 'opdPE' ) ~= 1 && ...
    strcmp( cfg.method, 'all' ) ~= 1  && cfg.plot == 1 )
  lineWidth = 0.2;
  for iChannel = 1:nChannels
    figure;
    forAxis1 = 1:numel( indata( iChannel, : ) );
    a = numel( indata( iChannel, : ) ) - numel( outdata( iChannel, : ) ) + 1;
    forAxis2 = a:numel( indata( iChannel, : ) );
    ax1 = subplot( 2, 1, 1 );
    if ( isfield( cfg, 'time' ) )
      if ( numel( cfg.time ) == length( indata ) )
        plot( cfg.time, indata( iChannel, : ), 'k', 'LineWidth', lineWidth ); grid on;
        if ( isfield( cfg, 'units' ) )
          xlabel( cfg.units );
        end
        if ( nChannels > 1 )
          title( [ 'Original time series, channel ' num2str( iChannel ) ] );
        else
          title( 'Original time series' );
        end
        ax2 = subplot( 2, 1, 2 );
        plot( cfg.time( end - length( outdata( iChannel, : ) ) + 1:end ), ...
              outdata( iChannel, : ), 'k', 'LineWidth', lineWidth ); grid on;
        linkaxes( [ ax1, ax2 ], 'x' );
      else
        plot( forAxis1, indata( iChannel, : ), 'k', 'LineWidth', lineWidth );
        if ( isfield( cfg, 'units' ) )
          xlabel( cfg.units );
        end
        if ( nChannels > 1 )
          title( ['Original time series, channel ' num2str( iChannel ) ] );
        else
          title( 'Original time series' );
        end
        ax2 = subplot( 2, 1, 2 );
        plot( forAxis2, outdata( iChannel, : ), 'k', 'LineWidth', lineWidth ); grid on;
        linkaxes( [ ax1, ax2 ], 'x' );
      end
    else
      plot( forAxis1, indata( iChannel, : ), 'k', 'LineWidth', lineWidth ); grid on;
      if ( isfield( cfg, 'units' ) )
        xlabel( cfg.units );
      end
      if ( nChannels > 1 )
          title( ['Original time series, channel ' num2str( iChannel ) ] );
      else
          title( 'Original time series' );
      end
      ax2 = subplot( 2, 1, 2 );
      plot( forAxis2, outdata( iChannel, : ), 'k', 'LineWidth', lineWidth ); grid on;
      linkaxes( [ ax1, ax2 ], 'x' );
    end

    if ( isfield( cfg, 'units' ) )
      xlabel( cfg.units );
    end
    if ( nChannels > 1 )
        title( [figureTitle ', channel ' num2str( iChannel ) ] );
    else
        title( figureTitle );
    end
  end
end
