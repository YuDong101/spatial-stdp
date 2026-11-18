function phases = reproduce_reps(mode, outputDir)
%REPRODUCE_REPS Recompute the six Fig. 1(b-g) rasters, in manuscript order.
%   archived uses all recorded spikes. resimulate starts each LIF network
%   from its spatial coordinates and random seeds, then verifies its spikes.
%   Both modes reconstruct phases, detect groups, sort neurons, and compare
%   each final matrix to the published figure. Returned phase units are cycles.
arguments
    mode (1,1) string {mustBeMember(mode,["archived","resimulate"])} = "archived"
    outputDir (1,1) string = ""
end
here = fileparts(mfilename('fullpath'));
originalPath = path;
cleanup = onCleanup(@() path(originalPath)); %#ok<NASGU>
addpath(fullfile(fileparts(here),'grid','lib'),'-begin');
gridMs = 15001:20000;
sample = unique([1:5:numel(gridMs),numel(gridMs)]);
phases = cell(1,6);
report = table('Size',[6 9],'VariableTypes', ...
    {'string','double','double','double','double','double','double','logical','string'}, ...
    'VariableNames',{'panel','alpha','gmax','weightSeed','phaseSeed','reference','stateCode','phaseMatch','phaseSHA256'});
for j = 1:6
    panel = char('a'+j);
    s = load(fullfile(here,'data',['panel_' panel '.mat']));
    fprintf('Figure 1(%c), %s: alpha=%g, gmax=%g, weightSeed=%d\n', ...
        panel, mode, s.p.alpha, s.p.gmax, s.p.weightSeed);
    if mode=="resimulate"
        inputs = lif_make_inputs(s.p,s.coordinates);
        result = lif_simulate(s.p,inputs);
        spikeMs = result.timesMs(:);
        node = result.nodeIds(:);
        assert(isequal(spikeMs,s.timesMs) && isequal(node,double(s.nodeIds)), ...
            'Figure 1(%c): simulated spikes differ from the archive.', panel);
        clear result
    else
        spikeMs = s.timesMs;
        node = double(s.nodeIds);
    end
    q = nan(s.p.N,numel(gridMs));
    for i = 1:s.p.N
        times = spikeMs(node==i);
        q(i,:) = interp1(times,0:numel(times)-1,gridMs,'linear',nan);
    end
    assert(all(isfinite(q),'all'),'Incomplete spikes within the observation window.');
    state = classify(q,gridMs,s.cfg);
    assert(state.code==s.classCode && state.reference==s.reference, ...
        'Figure 1(%c): state or reference differs from its published value.',panel);
    rel = mod(q-q(state.reference,:),1);
    % The reference rotates the phase origin. Rank is fixed throughout the window.
    meanPhase = mod(angle(mean(exp(2i*pi*rel),2))/(2*pi),1);
    rankGroup = ranking_groups(state,panel,s.p.N);
    [~,order] = sortrows([rankGroup meanPhase],[1 2]);
    phases{j} = rel(order,sample);
    md = java.security.MessageDigest.getInstance('SHA-256');
    md.update(typecast(phases{j}(:),'uint8'));
    digest = lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
    assert(strcmp(digest,s.phaseSHA256),'Figure 1(%c): phase raster differs from the manuscript.',panel);
    report(j,:) = {string(panel),s.p.alpha,s.p.gmax,s.p.weightSeed,s.p.phaseSeed, ...
        state.reference,state.code,true,string(digest)};
end
if strlength(outputDir)>0
    if ~isfolder(outputDir),mkdir(outputDir);end
    writetable(report,fullfile(outputDir,'Figure1_representatives.csv'));
end
end

function state = classify(q,gridMs,cfg)
initial = lif_phase_clusters_20261012(q);
bands = [.015 .02 .025 .03 .04];
groups = cell(size(bands));
for j = 1:numel(bands)
    z = lif_phase_clusters_20261012(q,struct('bandwidth',bands(j),'reference',initial.reference));
    groups{j} = struct('K',z.K,'cluster_id',z.cluster_id);
end
legacy = struct('groups',{groups},'reference',initial.reference,'D',z.D,'inputHash',"");
state = lif_n100_tree_state_20261014(q,gridMs,legacy,cfg);
end

function group = ranking_groups(state,panel,n)
if panel=='f'
    cores = state.core_diagnostics;
    cores = cores([cores.band]==4 & [cores.accepted]);
    assert(numel(cores)==2,'The coexistence example must have two verified cores at bandwidth 0.03.');
    group = 3*ones(n,1);
    for k = 1:numel(cores)
        group(cores(k).members) = k;
    end
    group(state.D) = 4;
else
    group = state.macro_id;
    if ~any(group),group = double(state.support_mask);end
    group(group==0) = max(group)+1;
    group(state.D) = max(group)+1;
end
end
