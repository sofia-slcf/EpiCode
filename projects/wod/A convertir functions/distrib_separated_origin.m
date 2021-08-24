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


config16 = wod_setparams;
config32 = wod_setparams_32chan;
config = config32;
cfgorig = config;
cfg= [config16 config32];

ipart= 1;

detectsavedir=fullfile(cfg{4}.imagesavedir,'detection');
detectiondatapath= fullfile(cfg{4}.datasavedir,'Detection');

temp= load(fullfile(detectiondatapath,'wod_wavedetection_allrat'));
stats_all=temp.stats_all;
clear temp


if ~isfolder(fullfile(detectsavedir,'histo'))
    mkdir(fullfile(detectsavedir,'histo'));
end



for irat= 1:size(cfg,2)
    
    if isempty(stats_all{irat})
        continue
    end
    
    
    irat_name=sprintf('Rat_%i',irat);
            for itrial= 1:size(stats_all{irat}.WoD.peak_time,2)
                %% WOD Find origin and store origin depth and timings
                
                A=stats_all{irat}.WoD.peak_time(:,itrial);
                B= stats_all{irat}.Depth(:,itrial);
                %find minimum timing
                origin_time=min(A);
               
                %find index of origin and get origin depth
                idx_origin=find(A==origin_time);
                origin_depth=stats_all{irat}.Depth(idx_origin,itrial);
                %store origin depth
                
                %security for 2 origin depth
                if size(origin_depth,1)>1
                    origin_depth=origin_depth(2);
                    A(idx_origin(1))=nan;
                    idx_origin=idx_origin(2);
                end
                
                if origin_depth > 1300 && origin_depth < 1700
                    origin_1400(irat,itrial)=origin_depth;
                elseif origin_depth < 700 || origin_depth > 1700
                    origin_others(irat,itrial)=origin_depth;
                elseif origin_depth>=800 && origin_depth<=1300
                    origin_1000(irat,itrial)=origin_depth;
                end
                
            end %itrial
end %irat

%% Delete 0 values and replace by NaNs
idx_zeros_1000=find(origin_1000==0);
origin_1000(idx_zeros_1000)=nan;
dist1000=reshape(origin_1000,[],1);

idx_zeros_1400=find(origin_1400==0);
origin_1400(idx_zeros_1400)=nan;
dist1400=reshape(origin_1400,[],1);

idx_zeros_others=find(origin_others==0);
origin_others(idx_zeros_others)=nan;
distothers=reshape(origin_others,[],1);

%% Histogram for origin around 1000 µm
fig_1000=figure;
histo_1000=histfit(dist1000,11,'kernel');
xlim([0 2000]);

fname_1000=fullfile(detectsavedir,'histo','histo_1000');
dtx_savefigure(fig_1000,fname_1000,'pdf','png','close');

%% Histogram for origin around 1400 µm
fig_1400=figure;
histo_1400=histfit(dist1400,11,'kernel');
xlim([0 2000]);


fname_1400=fullfile(detectsavedir,'histo','histo_1400');
dtx_savefigure(fig_1400,fname_1400,'pdf','png','close');

%% Histogram for outliers
fig_others=figure;
histo_others=histfit(distothers,11,'kernel');
xlim([0 2000]);

fname_others=fullfile(detectsavedir,'histo','histo_others');
dtx_savefigure(fig_others,fname_others,'pdf','png','close');

