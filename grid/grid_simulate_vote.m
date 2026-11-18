function ballot = grid_simulate_vote(settings, pointIndex, seed)
% Simulate one initial condition and apply the classifier used in Fig. 1(a).
parameters = settings.p;
parameters.alpha = settings.alpha(pointIndex);
parameters.gmax = settings.gmax(pointIndex);
parameters.weightSeed = seed;
parameters.phaseSeed = 10000 + seed;
inputs = lif_make_inputs(parameters,settings.coordinates);
result = lif_simulate(parameters,inputs);
assert(result.finalState.step * parameters.dtMs == parameters.endMs);
spikeTime = double(result.timesMs(:));
neuron = double(result.nodeIds(:));
time = 15001:20000; % samples in ms
phase = zeros(parameters.N,numel(time));
for node = 1:parameters.N
    events = spikeTime(neuron == node);
    phase(node,:) = interp1(events,0:numel(events)-1,time,'linear',NaN);
end
assert(all(isfinite(phase),'all'),'Spikes do not bracket the observation window.');

digest = java.security.MessageDigest.getInstance('SHA-256');
for values = {spikeTime,neuron,time(:)}
    digest.update(typecast(values{1},'uint8'));
end
phaseHash = lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2)',1,[]));
bandwidths = [.015 .02 .025 .03 .04];
initial = lif_phase_clusters_20261012(phase);
groups = cell(size(bandwidths));
for band = 1:numel(bandwidths)
    scale = lif_phase_clusters_20261012(phase,struct('bandwidth',bandwidths(band),'reference',initial.reference));
    groups{band} = struct('K',scale.K,'cluster_id',scale.cluster_id);
end
candidate = struct('groups',{groups},'reference',initial.reference,'D',scale.D,'inputHash',phaseHash);
state = lif_n100_tree_state_20261014(phase,time,candidate,settings.cfg);
[code,partitions,passing] = classify_boundaries(state,phase,parameters,pointIndex);
ballot = struct('point_index',pointIndex,'seed',seed,'original_code',state.code, ...
    'code',code,'partitions',partitions,'passing_partitions',passing,'phase_hash',phaseHash, ...
    'spike_steps',uint32(round(spikeTime / parameters.dtMs)),'neuron_ids',uint16(neuron));
end

function [code,partitions,passing] = classify_boundaries(state,phase,parameters,pointIndex)
code = state.code;
partitions = 0;
passing = 0;
if code ~= 0 || ~strcmp(state.reason,'member_dynamics_unresolved')
    return
end
assert(any(state.support_mask) && ~any(state.D) && ~any(state.escape_mask));
seen = {};
for index = 1:numel(state.partition_diagnostics)
    labels = state.partition_diagnostics{index}.initial_membership(:);
    assert(all(labels>0) && numel(labels)==size(phase,1));
    pairs = labels == labels';
    if any(cellfun(@(previous) isequal(previous,pairs),seen))
        continue
    end
    seen{end+1} = pairs; %#ok<AGROW>
    partitions = partitions + 1;
    passes = true;
    for group = unique(labels)'
        if ~stable_membership(phase(labels==group,:),state.config)
            passes = false;
            break
        end
    end
    passing = passing + passes;
end
if partitions>0 && passing==partitions
    code = 2;
    return
end
% The published map applied the positive dispersion recheck only to the
% reviewed weak-coupling rectangle and to the separately reviewed cell 2414.
reviewed = (parameters.alpha<=.48 && parameters.gmax>=.00014 && ...
    parameters.gmax<=.00053) || pointIndex==2414;
if reviewed && partitions>0 && passing==0 && min(state.dispersion_fraction)>=state.config.persistence
    code = 1;
end
end

function yes = stable_membership(phase,cfg)
yes = false;
n = size(phase,1);
for node = 2:n
    delta = phase(1:node-1,:) - phase(node,:);
    if any(max(delta,[],2)-min(delta,[],2) > cfg.coreSpan+1e-12)
        return
    end
end
ordered = sort(mod(phase,1),1);
k = ceil(.9*n);
ordered = [ordered;ordered+1];
width = min(ordered(k:k+n-1,:)-ordered(1:n,:),[],1);
inside = width<=cfg.coreWidth+1e-12;
half = floor(numel(inside)/2);
yes = min(mean(inside(1:half)),mean(inside(half+1:end))) >= cfg.persistence;
end
