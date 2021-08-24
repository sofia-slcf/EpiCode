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

if ispc
    addpath \\lexport\iss01.charpier\analyses\wod\fieldtrip-20200607
elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/wod/fieldtrip-20200607
end

ft_defaults


config16 = wod_setparams;
config32 = wod_setparams_32chan;

cfg=[config16 config32];

for irat=1:size(cfg,2)
    
    if isempty(cfg{irat})
        continue
    end
    for iwave= ["WoD" "WoR"]
        fig_scatter=figure;hold
        
        Color='rx';
        
        
        for itrial=1:size(stats_all{irat}.(iwave).peak_time,2)
            
            if itrial>1
                Color='bx';
            end
            
            data_plot.timing(:,1)=stats_all{irat}.(iwave).peak_time(:,itrial);
            data_plot.depth(:,1)=stats_all{irat}.Depth(:,itrial);
            
            scatter(data_plot.timing,data_plot.depth,Color);
            
            
            end %itrial
            
            scatterpath=fullfile(cfg{irat}.imagesavedir,'delays','scatter',sprintf('%s',cfg{irat}.prefix));
            
            
            
            if ~isfolder(scatterpath)
                mkdir(scatterpath);
            end
            
           
            fname_scatter=fullfile(scatterpath,sprintf('%s_scatter',(iwave)));
            
            dtx_savefigure(fig_scatter,fname_scatter,'pdf','png','close');
            clear data_plot
        
        end %iwave
    end %irat