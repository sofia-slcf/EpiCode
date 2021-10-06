function wod_create_Muse_Marker_file(cfg)

% cfg.rawdir
% cfg.muse.templatemarker
% cfg.directorylist{ipart}{idir}

if ~isfolder(cfg.rawdir)
    error('cfg.rawdir does not exist');
end

ipart = 1;

for idir = 1 : size(cfg.directorylist{ipart}, 2)
    
    %find the .nev file
    temp = dir(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},'*.nev'));
    if size(temp,1) > 1
        error('there should be only one nev file')
    end
    if size(temp,1) == 0
        nev_data{idir} = [];
        continue
    end
    
    % Verify that there is no Muse marker file (to avoid overwriting)
    [a, filename, c] = fileparts(temp.name);
    markerfile = fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, sprintf('%s.mrk', filename));
    if exist(markerfile, 'file')
        warning('A ''.mrk'' file is already present, nothing is done. \n%s', markerfile)
        continue
    end
    
    % Copy the muse template marker
    copyfile(cfg.muse.templatemarker, markerfile);
    
    % add TTL triggers (recorded during the experiment) to the Muse marker file
    cfgtemp = cfg;
    cfgtemp.directorylist = {};
    cfgtemp.directorylist{1}{1}  = cfg.directorylist{ipart}{idir};
    MuseStruct = readMuseMarkers(cfgtemp, true);
    MuseStruct = read_nev_Muse(cfgtemp, MuseStruct);
    writeMuseMarkerfile(MuseStruct{1}{1}, markerfile);
    
end