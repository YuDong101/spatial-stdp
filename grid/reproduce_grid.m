function [alpha,gmax,winner,fraction] = reproduce_grid(mode,outputDir,workerCount,pointIndices)
% Recompute Figure 1(a) from archived ballots or from independent LIF runs.
% [a,g,c,p]=reproduce_grid('archived',out) returns 3000 column-vector rows.
% reproduce_grid('full',out,24) resumes a complete simulation sweep.
% reproduce_grid('full',out,1,2414) restricts the sweep to a named cell.
% 'pilot' replays one specified trajectory; its voting fraction is undefined.
if nargin<1, mode='archived'; end
here = fileparts(mfilename('fullpath'));
if nargin<2, outputDir=fullfile(here,'..','output'); end
if ~isfolder(outputDir), mkdir(outputDir); end
dataDir = fullfile(here,'data');
settings = load(fullfile(dataDir,'scan_inputs.mat'));
reference = readtable(fullfile(dataDir,'map_reference.csv'));
assert(isequal(reference.point_index,(1:3000)') && ...
    max(abs(reference.alpha-settings.alpha))<1e-12 && max(abs(reference.gmax-settings.gmax))<1e-12);
if nargin<4, pointIndices=(1:height(reference))'; end
pointIndices = pointIndices(:);
validateattributes(pointIndices,{'double'},{'integer','positive','nonempty','<=',height(reference)});
assert(numel(pointIndices)==numel(unique(pointIndices)));
mode = validatestring(mode,{'archived','pilot','full'});
if strcmp(mode,'archived')
    votes = readtable(fullfile(dataDir,'trajectory_votes.csv'));
    trialTotals = accumarray(votes.point_index,1,[height(reference) 1]);
    boundaries = [0; cumsum(trialTotals)];
    counts = zeros(numel(pointIndices),6);
    for row = 1:numel(pointIndices)
        point = pointIndices(row);
        slice = boundaries(point)+1:boundaries(point+1);
        assert(all(votes.point_index(slice)==point) && isequal(votes.seed(slice),(1:numel(slice))'));
        counts(row,:) = replay(votes.state_code(slice));
    end
else
    if isempty(which('linkage')) || ~license('test','Statistics_Toolbox')
        error('Classification from spikes requires Statistics and Machine Learning Toolbox.');
    end
    previousPath = path;
    cleanupPath = onCleanup(@()path(previousPath));
    addpath(fullfile(here,'lib'),'-begin');
    checkpointDir = fullfile(outputDir,'grid_checkpoints');
    if ~isfolder(checkpointDir), mkdir(checkpointDir); end
    signature = code_signature(here);
    if strcmp(mode,'pilot')
        point = 2414;
        ballot = get_ballot(settings,point,4,checkpointDir,signature);
        expected = 'c88b617163f3423603d6c57cda66d8126da4f5a90a106c31832f5590901dc2b9';
        assert(strcmp(ballot.phase_hash,expected) && ballot.original_code==0 && ...
            ballot.code==1 && ballot.partitions==2 && ballot.passing_partitions==0, ...
            'The pilot did not reproduce the archived spike record and boundary classification.');
        alpha = settings.alpha(point); gmax = settings.gmax(point);
        winner = ballot.code; fraction = NaN;
        result = struct2table(rmfield(ballot,{'spike_steps','neuron_ids'}),'AsArray',true);
        writetable(result,fullfile(outputDir,'figure1_grid_pilot.csv'));
        return
    end
    counts = zeros(numel(pointIndices),6);
    if nargin<3 || isempty(workerCount)
        workerCount = 1;
        if license('test','Distrib_Computing_Toolbox')
            cluster = parcluster('Processes');
            workerCount = min([24,cluster.NumWorkers,numel(pointIndices)]);
        end
    end
    validateattributes(workerCount,{'numeric'},{'scalar','integer','positive'});
    if workerCount==1
        for row = 1:numel(pointIndices)
            counts(row,:) = simulate_cell(settings,pointIndices(row),checkpointDir,signature);
        end
    else
        pool = gcp('nocreate');
        if isempty(pool)
            pool = parpool('Processes',workerCount);
            cleanupPool = onCleanup(@()delete(pool));
        end
        assert(isa(pool,'parallel.ProcessPool'),'Use a process pool for independent simulations.');
        parfor (row=1:numel(pointIndices),workerCount)
            maxNumCompThreads(1);
            counts(row,:) = simulate_cell(settings,pointIndices(row),checkpointDir,signature);
        end
    end
end

alpha = settings.alpha(pointIndices);
gmax = settings.gmax(pointIndices);
n = sum(counts,2);
[largest,index] = max(counts,[],2);
winner = index-1;
fraction = largest./n;
assert(all(sum(counts==largest,2)==1));
expectedCounts = reference{pointIndices,compose('count_%d',0:5)};
matches = all(counts==expectedCounts,2);
if strcmp(mode,'archived')
    assert(all(matches) && all(winner==reference.modal_state(pointIndices)) && ...
        max(abs(fraction-reference.modal_fraction(pointIndices)))<1e-12, ...
        'The archived votes and manuscript map disagree.');
end
result = [table(pointIndices,alpha,gmax,winner,fraction,n,'VariableNames', ...
    {'point_index','alpha','gmax','modal_state','modal_fraction','trials'}), ...
    array2table(counts,'VariableNames',compose('count_%d',0:5)),table(matches)];
if numel(pointIndices)==height(reference)
    filename = ['figure1_grid_' mode '.csv'];
else
    % Retain completed checks from other subsets in this output directory.
    filename = ['figure1_grid_' mode '_subset.csv'];
    if isfile(fullfile(outputDir,filename))
        previous = readtable(fullfile(outputDir,filename));
        retain = ~ismember(previous.point_index,pointIndices);
        result = sortrows([previous(retain,:);result],'point_index');
    end
end
temporaryTable = [tempname(outputDir) '.csv'];
writetable(result,temporaryTable);
movefile(temporaryTable,fullfile(outputDir,filename),'f');
if any(~matches)
    warning('Fresh simulated ballots differ from the archived state counts at %d grid points; see the matches column.',nnz(~matches));
end
end

function counts = replay(codes)
assert(numel(codes)>=10 && all(ismember(codes,0:5)));
counts = zeros(1,6);
for seed = 1:numel(codes)
    counts(codes(seed)+1) = counts(codes(seed)+1)+1;
    if seed>=10 && nnz(counts==max(counts))==1
        assert(seed==numel(codes),'Ballots present beyond the specified stopping point.');
        return
    end
end
error('The archived sequence ended before the tie was resolved.');
end

function counts = simulate_cell(settings,point,folder,signature)
counts = zeros(1,6);
seed = 0;
while seed<10 || nnz(counts==max(counts))~=1
    seed = seed+1;
    b = get_ballot(settings,point,seed,folder,signature);
    counts(b.code+1) = counts(b.code+1)+1;
end
fprintf('Figure 1(a): cell %04d, %d trials.\n',point,seed);
end

function ballot = get_ballot(settings,point,seed,folder,signature)
sub = fullfile(folder,sprintf('p%04d',point));
if ~isfolder(sub), mkdir(sub); end
filename = fullfile(sub,sprintf('s%03d.mat',seed));
if isfile(filename)
    z = load(filename,'signature','ballot');
    assert(strcmp(signature,z.signature) && z.ballot.point_index==point && ...
        z.ballot.seed==seed && ismember(z.ballot.code,0:5) && ~isempty(z.ballot.phase_hash), ...
        'Checkpoint mismatch. Select a new output directory for a different configuration.');
    ballot = z.ballot;
    return
end
ballot = grid_simulate_vote(settings,point,seed);
tmp = [tempname(sub) '.mat'];
save(tmp,'signature','ballot','-v7');
movefile(tmp,filename,'f');
end

function signature = code_signature(here)
files = dir(fullfile(here,'lib','*.m'));
filenames = [fullfile({files.folder},{files.name}), ...
    {fullfile(here,'grid_simulate_vote.m'),fullfile(here,'data','scan_inputs.mat')}];
digest = java.security.MessageDigest.getInstance('SHA-256');
for k = 1:numel(filenames)
    fid = fopen(filenames{k},'rb');
    assert(fid~=-1); handle = onCleanup(@()fclose(fid));
    digest.update(fread(fid,Inf,'*uint8'));
    clear handle
end
signature = lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2)',1,[]));
end
