try %en local
    scriptpath = matlab.desktop.editor.getActiveFilename;
catch %cluster
    scriptpath = mfilename('fullpath');
end

epicodepath = [fileparts(fileparts(fileparts(fileparts(scriptpath)))), filesep];

addpath (genpath([epicodepath,'development']))
addpath (genpath([epicodepath,'shared']))
addpath (genpath([epicodepath,'external']))
addpath (genpath([epicodepath,'templates']))
addpath (genpath([epicodepath,'projects',filesep, 'wod']))
addpath (genpath([epicodepath,'projects',filesep, 'dtx']))
addpath (genpath([epicodepath,'projects',filesep, 'dtx',filesep,'dtx_functions']))
addpath (genpath([epicodepath,'projects',filesep, 'wod',filesep,'wod_functions']))
addpath (genpath([epicodepath,'projects',filesep, 'wod',filesep,'Intra']))

if ispc
    addpath \\lexport\iss01.charpier\analyses\wod\fieldtrip-20200607
    addpath \\lexport\iss01.charpier\echanges\scripts-paul\Spike2_vers_MATLAB
elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/wod/fieldtrip-20200607
    addpath /network/lustre/iss01/charpier/echanges/scripts-paul/Spike2_vers_MATLAB
    
end

ft_defaults

configDC= DC_setparams;
configIntra= Intra_setparams;


DC_filt.time= DC_filt.time{1}-DC_events.markers.VentOff.synctime;
DC_raw.time= DC_raw.time{1}-DC_events.markers.VentOff.synctime;

Intra.time=Intra.time{1}-Intra_events.markers.VentOff.synctime;


fig_intra=figure;hold;

plot(Intra.time,Intra.trial{1},'Color',[0.8500 0.3250 0.0980]);
plot(DC_raw.time,DC_raw.trial{1}+10,'Color',[0.8500 0.3250 0.0980]);
plot(DC_filt.time,DC_filt.trial{1}+10,'Color','r');

ylim([-100 30]);
xlim([-10 200]);

figpath=fullfile(configDC{14}.imagesavedir,'DC_intra_overdraw');

if ~isfolder(figpath)
    mkdir(figpath);
end
    
fname=fullfile(configDC{14}.imagesavedir,'DC_intra_overdraw',sprintf('Intra_DC_L23'));
dtx_savefigure(fig_intra,fname,'pdf','png','close');





