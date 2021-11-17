% This function rarefies a relative abundance otu table
% relabs is a matrix with relative abundances
% reads is the total number of reads in each sample, samples are in the
% same order than in relabs
% n_input is the number of reads to rarefied to

function out = AAD_rarefaction(relabs,reads,n_input)

% find dimension that has samples

dimsamples = find(size(relabs)==length(reads));
if isempty(dimsamples)
    error('dimension of reads vector does not match rows or columns in relabs')
elseif length(dimsamples)==2
    warning('relabs is square, assuming that rows are samples columns are taxa')
    dimsamples = 1;
end

% if samples are columns, transpose
if dimsamples==2
    relabs = relabs';
end

% make matrix with zeros
out = zeros(size(relabs));

for i = 1:size(relabs,1)
    % if size of sample is smaller than input, use all data
    n = min([n_input,reads(i)]);
    x=  relabs(i,:);
    if isnan(reads(i))
        out(i,:) = nan;
    elseif ~any(isnan(x))
        idx0 = find(x~=0);
        x0 = x(idx0);
        r = randperm(reads(i),n)/reads(i);
        rarefiedtemp = histcounts(r,[0,cumsum(x0)])/n;
        rarefiedx = zeros(size(x));
        rarefiedx(idx0) = rarefiedtemp;
        out(i,:) = rarefiedx/sum(rarefiedx);
    elseif all(isnan(x))
        %warning(['row ',num2str(i),' is all NaNs'])
    else
        error(['at least one element (and not all) of a sample (row ',num2str(i),' is NaNs'])
    end
end
end

