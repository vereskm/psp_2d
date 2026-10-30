function [localId,regions] = groupParameterSamples(samples,parameters,localParameters)
    %GROUPPARAMETERSAMPLES Assign rows to fixed Cartesian local-ROM regions.
    % Feature columns follow PARAMETERS. Match local features by device/field.
    % Region numbering follows ndgrid: the first local parameter varies fastest.
    % All regions are retained, including regions with no training samples.

    validateattributes(samples,{'numeric'},{'2d','real'},mfilename,'samples');
    assert(size(samples,2) == numel(parameters), ...
        'offline:groupParameterSamples:FeatureCount', ...
        'Sample columns must match the number of training parameters.');

    nLocal = numel(localParameters);
    featureIndices = zeros(1,nLocal);
    intervalCounts = zeros(1,nLocal);
    edgeSets = cell(1,nLocal);
    intervalId = zeros(size(samples,1),nLocal);

    for j = 1:nLocal
        parameter = localParameters(j);
        matches = arrayfun(@(p) isequal(p.device,parameter.device) && ...
            string(p.field) == string(parameter.field),parameters);
        featureIndex = find(matches);
        assert(isscalar(featureIndex), ...
            'offline:groupParameterSamples:FeatureMatch', ...
            ['Local parameter device %g, field %s must match exactly one ' ...
            'training parameter. Add it to training.parameters and ' ...
            'regenerate snapshots if it was not sampled.'], ...
            parameter.device,char(string(parameter.field)));
        assert(~ismember(featureIndex,featureIndices(1:j-1)), ...
            'offline:groupParameterSamples:DuplicateParameter', ...
            'A training feature cannot define two local partition axes.');

        edges = parameter.edges;
        validateattributes(edges,{'numeric'}, ...
            {'vector','real','finite'}, ...
            mfilename,'edges');
        assert(numel(edges) >= 2 && all(diff(edges(:)) > 0), ...
            'offline:groupParameterSamples:Edges', ...
            'Local parameter edges must contain at least two increasing values.');
        edges = edges(:).';
        values = samples(:,featureIndex);
        assert(all(isfinite(values)), ...
            'offline:groupParameterSamples:NonfiniteSamples', ...
            'Local parameter samples must be finite.');
        ids = discretize(values,edges);
        assert(all(~isnan(ids)), ...
            'offline:groupParameterSamples:OutsideEdges', ...
            'Some samples for device %g, field %s lie outside its local edges.', ...
            parameter.device,char(string(parameter.field)));

        featureIndices(j) = featureIndex;
        intervalCounts(j) = numel(edges)-1;
        edgeSets{j} = edges;
        intervalId(:,j) = ids;
    end

    nRegions = prod(intervalCounts); % prod([]) = 1: one global model
    localId = ones(size(samples,1),1);
    regionIntervals = zeros(nRegions,nLocal);
    lowerBounds = zeros(nRegions,nLocal);
    upperBounds = zeros(nRegions,nLocal);
    stride = 1;
    for j = 1:nLocal
        localId = localId + (intervalId(:,j)-1)*stride;
        indices = mod(floor((0:nRegions-1)'/stride),intervalCounts(j))+1;
        regionIntervals(:,j) = indices;
        edges = edgeSets{j};
        lowerBounds(:,j) = reshape(edges(indices),[],1);
        upperBounds(:,j) = reshape(edges(indices+1),[],1);
        stride = stride*intervalCounts(j);
    end

    regions.parameters = localParameters;
    regions.featureIndices = featureIndices;
    regions.intervalCounts = intervalCounts;
    regions.intervalIds = regionIntervals;
    regions.lowerBounds = lowerBounds;
    regions.upperBounds = upperBounds;
    regions.numRegions = nRegions;
end
