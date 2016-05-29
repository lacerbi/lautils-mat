% BINBUILD Build binned data.
% 
% Parameters:
% X             (1,nx) Vector of X data
% Y             (1,ny) Vector of Y data
% binCenters    (1,nbins) vector containing bin centers (optional)
% binSize       bin diameter (optional)
% outlierSD     remove outliers from each bin (default Inf, do not remove)
%
% Returns:
% binStruct     Structure containing binned data
%
% Typical usage for plotting mean biases:
% bCenters = -200:5:200; bSize = 60; [bins, bMeans, bStd] = binRawdata(rawData, bCenters, bSize); plot(bCenters, bMeans - bCenters);
%
function binStruct = binbuild(X, Y, binCenters, binSize, outliersd, robuststats)
    if ~exist('binCenters','var'); binCenters = []; end
    if ~exist('binSize','var'); binSize = []; end
    if ~exist('outliersd','var'); outliersd = []; end
    if ~exist('robuststats','var'); robuststats = []; end
    if isempty(robuststats); robuststats = 0; end

    % By default, go from min(X) to max(X), with a hundred bins
    if isempty(binCenters); binCenters = min(X):(max(X)-min(X))/99:max(X); end    
    
    % By default, bins are quite large (ten times the average bin step)
    if isempty(binSize); binSize = 0.5*[diff(binCenters(1:2)) diff(binCenters)] + 0.5*[diff(binCenters) diff(binCenters(end-1:end))]; end
    if isscalar(binSize); binSize = binSize*ones(1, length(binCenters)); end

    % By default, do not remove outliers
    if isempty(outliersd); outliersd = Inf; end
        
    % Number of bins
    nbins = size(binCenters, 2);
    
    statstring = {'ymeans', 'ystd', 'yskew', 'ykurt'};
    robuststatstring = {'ymedian', 'ykmad', 'yinterq', 'yinterq1', 'yinterq2'};
    
    % Initialize stats
    for k = 1:length(statstring)
        binstruct.(statstring{k}) = NaN(1, nbins);
    end
    if robuststats
        for k = 1:length(robuststatstring)
            binstruct.(robuststatstring{k}) = NaN(1, nbins);
        end
    end    
    binStruct.n = zeros(1, nbins);
        
    for i = 1:nbins        
       binStruct.y{i} = [];        
       binc = binCenters(i);
       
       f = ((X - binc) > -0.5*binSize(i)) & ((X - binc) <= 0.5*binSize(i));
       binStruct.y{i} = [binStruct.y{i} Y(f)];              
       
       % Compute bin statistics
       if ~isempty(binStruct.y{i}); fillBinStats(); end
       binStruct.n(i) = length(binStruct.y{i});
       
       % Remove outliers and recompute statistics
       if ~isinf(outliersd)
            if ~isempty(binStruct.y{i})
                binStruct.y{i}(abs(binStruct.y{i} - binStruct.ymeans(i)) > outliersd*binStruct.ystd(i)) = [];
            end
           if isempty(binStruct.y{i})
                for k = 1:length(statstring)
                    binstruct.(statstring{k})(i) = NaN;
                end
                if robuststats
                    for k = 1:length(robuststatstring)
                        binstruct.(robuststatstring{k})(i) = NaN;
                    end
                end    
           else
               fillBinStats();
           end
       end
       binStruct.discarded(i) = binStruct.n(i) - length(binStruct.y{i});
       binStruct.n(i) = length(binStruct.y{i});
       
    end
    binStruct.x = binCenters;
    binStruct.size = binSize;
    
    return;
    
    % Compute bin statistics
    function fillBinStats()
       binStruct.ymeans(i) = nanmean(binStruct.y{i});
       binStruct.ystd(i) = nanstd(binStruct.y{i});
       binStruct.yskew(i) = skewness(binStruct.y{i});
       binStruct.ykurt(i) = kurtosis(binStruct.y{i})-3;
       if robuststats
           binStruct.ymedian(i) = median(binStruct.y{i});
           binStruct.yinterq1(i) = quantile(binStruct.y{i}, 0.25);
           binStruct.yinterq2(i) = quantile(binStruct.y{i}, 0.75);
           binStruct.yinterq(i) = binStruct.yinterq2(i) - binStruct.yinterq1(i);
           binStruct.ykmad(i) = 1.4826*median(abs(binStruct.y{i} - binStruct.ymedian(i)));
       end        
    end
end