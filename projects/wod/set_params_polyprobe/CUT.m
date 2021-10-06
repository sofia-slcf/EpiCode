%% read data
clear all
%import les données 'experience' et formes-en une matrice
Table_experience = readtable('\\lexport\iss01.charpier\raw\rat-wod\2probes\Data_table\Tableau_Experience_per_rat.xlsx')
Table_exp = table2cell(Table_experience)
Table_experience=[Table_exp(:,2:22)] 
Table_experience= cell2mat(Table_experience)
folder = Table_exp(:,1) 


protocol = dir('\\lexport\iss01.charpier\raw\rat-wod\2probes\Data_table\grid');
protocol = protocol(3:end);    
tableau.each=[]


for irat=1:size(Table_experience,1) %pour chaque les rats
    
    %selection le protocole 1 parmis les protocoles disponobles pour 1 rat
    rat_name = cell2mat(folder(irat))
    grid = [rat_name '1']
    
    %importe la gride de données spike2 du rat WOD1
    opts = detectImportOptions([protocol(irat).folder,'/',grid, '_grid.xlsx']);
    opts = setvartype(opts,'single');
    tableau(irat).each = readtable([protocol(irat).folder,'/',grid, '_grid.xlsx'], opts);
    tableau(irat).each=cell2mat(table2cell(tableau(irat).each))
    
    %prend les donnée d'experience
       maxisize = size(tableau(irat).each,1)
       Nrat(1:maxisize) = Table_experience(irat,1);
       
       if size(Nrat,1)==1
          Nrat=Nrat.'
       else
          Nrat=Nrat
       end
       
       Trial(1:maxisize)  = Table_experience(irat,2);
       if size(Trial,1)==1
          Trial=Trial.'
       else
          Trial=Trial
       end
       
       Recov = zeros(maxisize,2)  %= Table_experience(irat,3);
       Recov(:,1)= Table_experience(irat,3)
       Recov(:,2)= Table_experience(irat,4)
       
     
       ProbA=[]
       ProbB=[]
           for i=1:size(tableau(irat).each,1)
           if tableau(irat).each(i,1) >=17
               ProbA(i)=1
               ProbB(i)=0
           else
               ProbA(i)=0
               ProbB(i)=1
           end
           end
       len_ProbA = length(find(ProbA == 1))
       len_ProbB = length(find(ProbB == 1))
       ProbA=ProbA.'
       ProbB=ProbB.'
       
       ZE01probA = ones(len_ProbA, 1) * Table_experience(irat,13);
       ZE01probB = ones(len_ProbB, 1) * Table_experience(irat,20);
       ZE01=[ZE01probA;ZE01probB ]
      
       % insert les données experience dans un tableau maxi 
       tableau(irat).each = [ Nrat tableau(irat).each(:,1:end)]
       tableau(irat).each = [ tableau(irat).each(:,1) Trial tableau(irat).each(:,2:end)]
       tableau(irat).each = [ tableau(irat).each(:,1:2) Recov tableau(irat).each(:,3:end)]
       tableau(irat).each = [ tableau(irat).each(:,1:4) ProbA tableau(irat).each(:,5:end)]
       tableau(irat).each = [ tableau(irat).each(:,1:5) ProbB tableau(irat).each(:,6:end)]
       tableau(irat).each = [ tableau(irat).each(:,1:7) ZE01 tableau(irat).each(:,8:end)]
       % take the right stereotaxie
       
       
       
        X_bl = Table_experience(irat,5)
        Y_bl = Table_experience(irat,6)
        X_PA = Table_experience(irat,11)
        Y_PA = Table_experience(irat,12)
        X_PB = Table_experience(irat,18)
        Y_PB = Table_experience(irat,19)

        name(irat).x_pa  = ones(len_ProbA,1) * stereo(X_bl,Y_bl,X_PA,Y_PA,X_PB,Y_PB).x_pa;
        name(irat).y_pa = ones(len_ProbA,1) * stereo(X_bl,Y_bl,X_PA,Y_PA,X_PB,Y_PB).y_pa;
        name(irat).x_pb = ones(len_ProbB,1) * stereo(X_bl,Y_bl,X_PA,Y_PA,X_PB,Y_PB).x_pb;
        name(irat).y_pb   = ones(len_ProbB,1) * stereo(X_bl,Y_bl,X_PA,Y_PA,X_PB,Y_PB).y_pb;
       
        X_P = [name(irat).x_pa ;name(irat).x_pb]
        Y_P = [name(irat).y_pa ;name(irat).y_pb]
       
        tableau(irat).each = [ tableau(irat).each(:,1:7) X_P tableau(irat).each(:,8:end)]
        tableau(irat).each = [ tableau(irat).each(:,1:8) Y_P tableau(irat).each(:,9:end)]
       
       %%%% to change
       X_ch = zeros(1,len_ProbA+len_ProbB).';
       Y_ch = zeros(1,len_ProbA+len_ProbB).';
       tableau(irat).each = [ tableau(irat).each(:,1:10) X_ch tableau(irat).each(:,11:end)]
       tableau(irat).each = [ tableau(irat).each(:,1:11) Y_ch tableau(irat).each(:,12:end)]
       %%%%
        
        
        
        %positon electrode pour les probe A(1)=>pta  et B(2)=>s1
       % probeA
        X_bl         = Table_experience(irat,5); %commun pour probA et B 
        Z_bl         = Table_experience(irat,7); %commun pour probA et B 
        
        zE01_A       = Table_experience(irat,13);
        angle_A      = Table_experience(irat,14);
        nb_chan_A    = [1:Table_experience(irat,9)]
        inter_chan_A = Table_experience(irat,10)

        name(irat).name_channel_A   =  str2double(Table_A(angle_A, zE01_A, nb_chan_A, inter_chan_A, X_bl, Z_bl).prob_A.name_channel);
        name(irat).rename_channel_A = str2double(Table_A(angle_A, zE01_A, nb_chan_A, inter_chan_A, X_bl, Z_bl).prob_A.rename_channel);
        name(irat).absolute_depth_A = cell2mat(Table_A(angle_A, zE01_A, nb_chan_A, inter_chan_A, X_bl, Z_bl).prob_A.absolute_depth);
        name(irat).structure_A    = cell2mat(Table_A(angle_A, zE01_A, nb_chan_A, inter_chan_A, X_bl, Z_bl).prob_A.structure) ;

        %probeB
        zE01_B = Table_experience(irat,20)
        angle_B = Table_experience(irat,21)
        nb_chan_B = [1:Table_experience(irat,16)]
        inter_chan_B = Table_experience(irat,17)
        
        name(irat).name_channel_B   =  str2double(Table_B(angle_B, zE01_B, nb_chan_B, inter_chan_B, X_bl, Z_bl).prob_B.name_channel);
        name(irat).rename_channel_B = str2double(Table_B(angle_B, zE01_B, nb_chan_B, inter_chan_B, X_bl, Z_bl).prob_B.rename_channel);
        name(irat).absolute_depth_B = cell2mat(Table_B(angle_B, zE01_B, nb_chan_B, inter_chan_B, X_bl, Z_bl).prob_B.absolute_depth);
        name(irat).structure_B    = cell2mat(Table_B(angle_B, zE01_B, nb_chan_B, inter_chan_B, X_bl, Z_bl).prob_B.structure) ;
 
        name(irat).absolute_depth= [name(irat).absolute_depth_A;name(irat).absolute_depth_B];
        name(irat).rename_channel = [name(irat).rename_channel_A; name(irat).rename_channel_B];
        name(irat).structure = [name(irat).structure_A ;name(irat).structure_B];
        name(irat).sub_structure =  zeros(1,maxisize).'% à remplir à la main
        
        tableau(irat).each = [ tableau(irat).each(:,1:12) name(irat).sub_structure tableau(irat).each(:,13:end)];
        
        % NaN resahpe MAXI
        
        tableau(irat).eachshaped= NaN(48,23)
       nbNaN=[]
       for row=1:size(tableau(irat).eachshaped,1)
           ligne= row-length(nbNaN)
           array=49-row
           if tableau(irat).each(ligne,7).' == array
              tableau(irat).eachshaped(row,:) = [ tableau(irat).each(ligne,:)]
           else
              tableau(irat).eachshaped(row,:) = NaN
           end
           nbNaN= find(isnan(tableau(irat).eachshaped(1:row,6)))
       end

       
        
        tableau(irat).eachshaped = [ tableau(irat).eachshaped(:,1:12) name(irat).absolute_depth tableau(irat).eachshaped(:,13:end)];
        tableau(irat).eachshaped = [ tableau(irat).eachshaped(:,1:13) name(irat).rename_channel tableau(irat).eachshaped(:,14:end)];
        tableau(irat).eachshaped = [ tableau(irat).eachshaped(:,1:14) name(irat).structure tableau(irat).eachshaped(:,15:end)];
       


savedir= '\\lexport\iss01.charpier\raw\rat-wod\2probes\Data_table\maxi_cut\'
maxicut = tableau(irat).eachshaped(:,1:16)
save(fullfile(savedir, [rat_name]), 'maxicut')


clearvars -except  protocol folder Table_experience tableau.each


end













