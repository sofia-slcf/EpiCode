try %en local
    scriptpath = matlab.desktop.editor.getActiveFilename;
catch %cluster
    scriptpath = mfilename('fullpath');
end

epicodepath = [fileparts(fileparts(fileparts(scriptpath))), filesep];

addpath (genpath([epicodepath,'development']))
addpath (genpath([epicodepath,'shared']))
addpath (genpath([epicodepath,'external']))
addpath (genpath([epicodepath,'templates']))
addpath (genpath([epicodepath,'projects', filesep, 'wod']))
addpath (genpath([epicodepath,'projects', filesep, 'dtx']))
addpath (genpath([epicodepath,'projects', filesep, 'wod',filesep,'wod_functions']))
addpath (genpath([epicodepath,'projects', filesep, 'dtx',filesep,'dtx_functions']))

if ispc
    addpath \\lexport\iss01.charpier\analyses\wod\fieldtrip-20200607
elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/wod/fieldtrip-20200607
end

ft_defaults


figure_depths=figure;
% Plot 1000
irow=0;
for irat=1:size(stats_all,2)
    if isempty(stats_all{irat})
        continue
    end
    for itrial=1:size(stats_all{irat}.oridepthclass,2)
       if stats_all{irat}.oridepthclass(itrial)==1000
           irow=irow+1;
Value_1000(irow) = stats_all{irat}.wod_origin_depth(itrial);
       end
    end
end

x = ones(1, length(Value_1000));
scatter(x, Value_1000, 'ko',  'LineWidth', 1,'MarkerFaceColor','w','jitter', 'on', 'jitterAmount', 0.02);
hold on;
% Plot 14000
irow=0;
for irat=1:size(stats_all,2)
    if isempty(stats_all{irat})
        continue
    end
    for itrial=1:size(stats_all{irat}.oridepthclass,2)
       if stats_all{irat}.oridepthclass(itrial)==1400
           irow=irow+1;
Value_1400(irow) = stats_all{irat}.wod_origin_depth(itrial);
       end
    end
end
x = 1.2 * ones(1, length(Value_1400));
scatter(x, Value_1400, 'ko',  'LineWidth', 1,'MarkerFaceColor','w','jitter', 'on', 'jitterAmount', 0.02);
% Set up axes.
xlim([0.5, 1.5]);
ylim([700, 2000]);
ax = gca;
ax.XTick = [1, 1.2];
ax.XTickLabels = {'1000','1400'};

fname=



scatter(x, Value_1000, 'ko',  'LineWidth', 1,'MarkerFaceColor','w','jitter', 'on', 'jitterAmount', 0.05);
