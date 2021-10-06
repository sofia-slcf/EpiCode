%scripts d'analyse
addpath '\\lexport\iss01.charpier\analyses\tms\raw\';

%fieldtrip
addpath '\\lexport\iss01.charpier\analyses\tms\scripts\fieldtrip';
ft_defaults

%toolbox ced
addpath '\\lexport\iss01.charpier\echanges\scripts-paul\Spike2_vers_MATLAB\';
addpath '\\lexport\iss01.charpier\echanges\scripts-paul\Spike2_vers_MATLAB\CED_library\CEDS64ML';
CEDS64LoadLib('\\lexport\iss01.charpier\echanges\scripts-paul\Spike2_vers_MATLAB\CED_library\CEDS64ML');

%% paramètres :
datapath{ineuron}     = ('\\lexport\iss01.charpier\analyses\wod\Antoine\data\Cross_cor_from_matlab\test_crosscorr.smrx');
channame{1}    = 'Vm';
new_chan_number{1} = 107;
new_chan_name{1}   = 'Vm';

for ineuron = 1:size(datapath,1)
    
    % lire le canal d'intérêt
    data                = readCEDwaveforms(datapath{ineuron}, channame{ineuron});
    %ouvrir le fichier Spike2 pour pouvoir le modifier
    fid         = CEDS64Open(datapath{ineuron},0);
    
    
    %% faire les modifs voulues sur les données
    
    %read LF events
    Events{ineuron}=readCEDevents(fullfile(config{ineuron}.rawdir,config{ineuron}.directorylist{1}{1}));
    
    %define period to do xcorr on
    t_1= Events{ineuron}.markers.LF.synctime(1);
    t_2= Events{ineuron}.markers.LF.synctime(2);
    t_sel=[t_1 t_2];
    
    cfgtemp=[];
    cfgtemp.latency=t_sel;
    LF_Vm{ineuron}=ft_selectdata(cfgtemp,Vm_APless{ineuron});
    LF_EEG{ineuron}=ft_selectdata(cfgtemp,EEG{ineuron});
 
    data_modif                    = LF_Vm{ineuron};
    
    %align to Voff
    data_modif.time{1} = data_modif.time{1}(1) - temps_mrk;
    
    %vérifier le résultat
    % plot(data_detrend.time{1}, data_detrend.trial{1});
    
    %% écrire le nouveau canal dans le fichier Spike2
    
    %convert time from seconds to ticks
    start_time = CEDS64SecsToTicks(fid,data_modif.time{1}(1));
    waveform   = int16(data_modif.trial{1}*100); %*100 pour garder une précision des integer.
    
    %écrire le nouveau canal dans le fichier Spike2, et fermer le fichier
    i64Div = 1.0/(data_modif.fsample*CEDS64TimeBase(fid)); %trouver le nombre de divisions entre 2 points
    check = CEDS64SetWaveChan(fid,new_chan_number{1},i64Div,1,data_modif.fsample); %créer un nouveau canal
    check = CEDS64ChanTitle(fid, new_chan_number{1}, new_chan_name{1}); %renommer le nouveau canal
    check = CEDS64WriteWave(fid,new_chan_number{1},waveform,start_time); %écrire les données corrigées dans ce canal
    CEDS64Close(fid); %ferme le fichier spike2
    
end

unloadlibrary ceds64int;
