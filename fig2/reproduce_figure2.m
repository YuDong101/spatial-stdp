function [metrics, validation] = reproduce_figure2(mode, outputDir)
%REPRODUCE_FIGURE2 Recompute both members of the paired formation experiment.
%   reproduce_figure2('archived', out) uses recorded spikes and weight frames.
%   reproduce_figure2('resimulate', out) regenerates inputs and integrates
%   both arms before performing exactly the same measurements and rendering.
if nargin < 1, mode = 'archived'; end
here = fileparts(mfilename('fullpath'));
if nargin < 2, outputDir = fullfile(here, '..', 'output'); end
mode = validatestring(mode, {'archived','resimulate'});
previousPath = path;
restorePath = onCleanup(@()path(previousPath)); %#ok<NASGU>
addpath(fullfile(here,'..','grid','lib'),'-begin');
if ~isfolder(outputDir), mkdir(outputDir); end
data = load(fullfile(here,'data','paired_records.mat'));
arms = {'S','N'};
for i = 1:2
    record = data.raw(i);
    p = record.parameters;
    assert(strcmp(p.arm,arms{i}));
    frames = record.weights;
    if isempty(frames)
        % Static weights in the no-learning arm are stored once, losslessly.
        frames = repmat(reshape(data.W0,1,p.N,p.N),numel(data.checkpointsMs),1,1);
    end
    if strcmp(mode,'resimulate')
        x = lif_make_inputs(p,data.coordinates);
        assert(isequal(x.W,data.W0),'Initial weight regeneration failed.');
        simulated = lif_simulate(p,x);
        assert(isequal(simulated.timesMs,record.timesMs) && ...
               isequal(simulated.nodeIds,double(record.nodeIds)), ...
               'Regenerated %s spikes differ from the archived run.',p.arm);
        assert(isequal(simulated.frameTimesMs,data.checkpointsMs) && ...
               max(abs(simulated.weights-frames),[],'all') < 1e-13, ...
               'Regenerated %s weights differ from the archived run.',p.arm);
        record.timesMs = simulated.timesMs;
        record.nodeIds = simulated.nodeIds;
        frames = simulated.weights;
    end
    metrics.(p.arm) = measure_formation(record,frames,data.checkpointsMs);
end
metrics.initialRelativeWeight = data.W0/data.raw(1).parameters.gmax;
metrics.snapshotTimesMs = [100 3320 6900 20000];
% Stable ordering: descending module size, then original neuron number.
[~,metrics.nodeOrder] = sortrows([metrics.S.modules,(1:data.raw(1).parameters.N)']);
fig = render_formation(metrics,fullfile(here,'data','layout.fig'));
cleanup = onCleanup(@()close(fig)); %#ok<NASGU>

% Only the audit reads the approved panel values. It accesses the actual
% plotted graphics objects, including the rank, point colors and boundaries.
oracle = load(fullfile(here,'data','reference_panels.mat'),'reference');
validation = compare_panels(fig,oracle.reference);
save(fullfile(outputDir,'Figure2_metrics.mat'),'metrics');
writetable(validation,fullfile(outputDir,'Figure2_validation.csv'));
exportgraphics(fig,fullfile(outputDir,'Figure2.pdf'),'ContentType','vector','BackgroundColor','white');
exportgraphics(fig,fullfile(outputDir,'Figure2.png'),'Resolution',300,'BackgroundColor','white');
fprintf('Figure 2: %s, %d panel objects verified; largest absolute difference %.3g.\n', ...
    mode,height(validation),max(validation.maxAbsError));
end

function m = measure_formation(record,frames,frameTimes)
% Linear interpolation of each neuron's cumulative spike count.
grid = 100:10:20000;
n = record.parameters.N;
q = zeros(n,numel(grid));
for node = 1:n
    spikes = double(record.timesMs(record.nodeIds==node));
    assert(numel(spikes)>=2 && all(diff(spikes)>0));
    q(node,:) = interp1(spikes,0:numel(spikes)-1,grid,'linear',NaN);
end
assert(all(isfinite(q),'all'),'Spikes must bracket the formation window.');
phaseVector = exp(2*pi*1i*q);
m.timeMs = grid;
m.phase = q;
m.R1 = abs(mean(phaseVector,1));
m.R2 = abs(mean(phaseVector.^2,1));
m.frameTimesMs = frameTimes;
m.Q = zeros(size(frameTimes));
endFrame = find(frameTimes==20000,1);
assert(~isempty(endFrame));
for k = 1:numel(frameTimes)
    W = squeeze(frames(k,:,:));
    % BCT's bundled spectral algorithm is deterministic (first optimal
    % swap). Identical adjacent frames therefore share one calculation.
    if k==1 || ~isequal(W,previousW)
        [labels,Q] = modularity_dir(W,1);
    end
    m.Q(k) = Q;
    previousW = W;
    if k==endFrame
        m.finalRelativeWeight = W/record.parameters.gmax;
        original = labels(:);
        originalIds = unique(original);
        sizes = arrayfun(@(v)sum(original==v),originalIds);
        [~,idx] = sortrows([-sizes(:),originalIds(:)]);
        m.modules = zeros(n,1);
        for j = 1:numel(idx)
            m.modules(original==originalIds(idx(j))) = j;
        end
    end
end
m.concentration = zeros(max(m.modules),numel(grid));
for j = 1:size(m.concentration,1)
    members = sort(mod(q(m.modules==j,:),1),1);
    count = size(members,1);
    wrapped = [members;members+1];
    largest = ones(1,numel(grid));
    for offset = 1:count-1
        span = wrapped(1+offset:count+offset,:)-members;
        largest(any(span<=0.05+1e-12,1)) = offset+1;
    end
    m.concentration(j,:) = largest/count;
end
end

function fig = render_formation(m,layoutFile)
fig = openfig(layoutFile,'new','invisible');
arms = {'S','N'};
top = {'a','b'};
snapshots = {'f','g'};
bottom = {'i','j'};
nodes = m.nodeOrder;
labels = m.S.modules(nodes);
colors = [0 0.4470 0.7410;0.8500 0.3250 0.0980]; % published module palette
assert(all(labels<=size(colors,1)));
for i = 1:2
    arm = m.(arms{i});
    time = arm.timeMs/1000;
    set_line(fig,top{i},1,time,arm.R1);
    set_line(fig,top{i},2,time,arm.R2);
    set_line(fig,'h',i,arm.frameTimesMs/1000,arm.Q);
    for k = 1:numel(m.snapshotTimesMs)
        col = arm.timeMs==m.snapshotTimesMs(k);
        assert(nnz(col)==1);
        points = tagged(fig,sprintf('%s%d_points',snapshots{i},k));
        set(points,'XData',1:numel(nodes),'YData',mod(arm.phase(nodes,col),1), ...
            'CData',colors(labels,:));
    end
    sigma = arm.concentration;
    for j = 1:size(sigma,1)
        set_line(fig,bottom{i},j,time,sigma(j,:));
    end
    set_line(fig,bottom{i},size(sigma,1)+1,time,mean(sigma,1));
end
matrices = {m.initialRelativeWeight,m.S.finalRelativeWeight,m.N.finalRelativeWeight};
panels = 'cde';
assert(max(labels)==2);
edge = nnz(labels==1)+0.5;
for j = 1:3
    panel = panels(j);
    imageHandle = tagged(fig,[panel '_image']);
    imageHandle.CData = matrices{j}(nodes,nodes);
    set(tagged(fig,[panel '_boundary1']),'Value',edge);
    set(tagged(fig,[panel '_boundary2']),'Value',edge);
end
drawnow;
end

function set_line(fig,panel,series,x,y)
lineHandle = tagged(fig,sprintf('%s_line%d',panel,series));
set(lineHandle,'XData',x,'YData',y);
end

function obj = tagged(fig,name)
obj = findobj(fig,'Tag',name);
assert(isscalar(obj),'Missing or duplicate figure object: %s.',name);
end

function checks = compare_panels(fig,reference)
fields = fieldnames(reference);
errors = zeros(numel(fields),1);
for k = 1:numel(fields)
    tag = fields{k};
    h = tagged(fig,tag);
    switch h.Type
        case 'line', value = [h.XData(:),h.YData(:)];
        case 'image', value = h.CData;
        case 'scatter', value = [h.XData(:),h.YData(:),h.CData];
        case 'constantline', value = h.Value;
        otherwise, error('Unexpected graphic type %s in %s.',h.Type,tag);
    end
    ref = reference.(tag);
    assert(isequal(size(ref),size(value)),'Size mismatch in %s.',tag);
    errors(k) = max(abs(double(value(:))-double(ref(:))));
    assert(errors(k)<1e-10,'Mismatch in %s: max error %.15g.',tag,errors(k));
end
checks = table(string(fields),errors,'VariableNames',{'graphic','maxAbsError'});
end
