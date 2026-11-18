function outputDir = run_figures12(mode, outputDir)
%RUN_FIGURES12 Run this file in MATLAB to reconstruct manuscript Figures 1 and 2.
%   run_figures12                   Recalculate from archived ballots/spikes/weights.
%   run_figures12('resimulate')     Reintegrate the examples and S/N pair; use archived map votes.
%   run_figures12('pilot')          Reintegrate a grid trajectory and audit its classification.
%   run_figures12('full')           Reintegrate the full map, examples, and S/N pair.
%   run_figures12(mode, outputDir)  Choose an output folder.
arguments
    mode (1,1) string {mustBeMember(mode,["archived","resimulate","pilot","full"])} = "archived"
    outputDir (1,1) string = ""
end
root = fileparts(mfilename('fullpath'));
if strlength(outputDir) == 0
    outputDir = fullfile(root,'output',mode);
end
if ~isfolder(outputDir), mkdir(outputDir); end
previousPath = path;
previousRng = rng;
restorePath = onCleanup(@() path(previousPath)); %#ok<NASGU>
restoreRng = onCleanup(@() rng(previousRng)); %#ok<NASGU>
addpath(fullfile(root,'function'),fullfile(root,'grid'),fullfile(root,'reps'),fullfile(root,'fig2'));
fprintf('Figures 1 and 2 | MATLAB %s | %s\n',version('-release'),mode);
if mode == "pilot"
    reproduce_grid('pilot',outputDir);
    fprintf('Grid pilot complete. Output: %s\n',outputDir);
    return
end
gridMode = "archived";
exampleMode = "archived";
if mode == "full", gridMode = "full"; end
if mode ~= "archived", exampleMode = "resimulate"; end

[alpha,gmax,winner,fraction] = reproduce_grid(gridMode,outputDir);
phases = reproduce_reps(exampleMode,outputDir);
plot_figure1(alpha,gmax,winner,fraction,phases,outputDir);
reproduce_figure2(exampleMode,outputDir);
fprintf('Complete. Open Figure1.pdf and Figure2.pdf in %s\n',outputDir);
end
