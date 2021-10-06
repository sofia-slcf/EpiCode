%% read data
clear all
%import les données 'experience' et formes-en une matrice
Table_experience = readtable('C:\Users\s.carrionfalgarona\Desktop\Tableau_Experience.xlsx')
Table_exp = table2cell(Table_experience)
Table_experience=[Table_exp(:,2:21)] 
Table_experience= cell2mat(Table_experience)
%definir un nombre de rat comme le nombre de grid diponibles
rat = dir('C:\Users\s.carrionfalgarona\Desktop\grid');
rat=rat(3:end);
     
tableau.each=[]

for irat=1:size(rat,1) %pour chaque les rats
    
    %importe la gride de données spike 2 du rat
    opts = detectImportOptions([rat(irat).folder,'/',rat(irat).name]);
    opts = setvartype(opts,'single');
    tableau(irat).each = readtable([rat(irat).folder,'/',rat(irat).name], opts);
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

       Recov(1:maxisize)  = Table_experience(irat,3);
       
       if size(Recov,1)==1
          Recov=Recov.'
       else
          Recov = Recov
       end
       
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
       
       ZE01probA = ones(len_ProbA, 1) * Table_experience(irat,12);
       ZE01probB = ones(len_ProbB, 1) * Table_experience(irat,19);
       ZE01=[ZE01probA;ZE01probB ]
      
       % insert les données experience dans un tableau maxi 
       tableau(irat).each = [ Nrat tableau(irat).each(:,1:end)]
       tableau(irat).each = [ tableau(irat).each(:,1) Trial tableau(irat).each(:,2:end)]
       tableau(irat).each = [ tableau(irat).each(:,1:2) Recov tableau(irat).each(:,3:end)]
       tableau(irat).each = [ tableau(irat).each(:,1:3) ProbA tableau(irat).each(:,4:end)]
       tableau(irat).each = [ tableau(irat).each(:,1:4) ProbB tableau(irat).each(:,5:end)]
       tableau(irat).each = [ tableau(irat).each(:,1:6) ZE01 tableau(irat).each(:,7:end)]
       % take the right stereotaxie
       
       
       
        X_bl = Table_experience(irat,4)
        Y_bl = Table_experience(irat,5)
        X_PA = Table_experience(irat,10)
        Y_PA = Table_experience(irat,11)
        X_PB = Table_experience(irat,17)
        Y_PB = Table_experience(irat,18)

        name(irat).x_pa  = ones(len_ProbA,1) * stereo(X_bl,Y_bl,X_PA,Y_PA,X_PB,Y_PB).x_pa;
        name(irat).y_pa = ones(len_ProbA,1) * stereo(X_bl,Y_bl,X_PA,Y_PA,X_PB,Y_PB).y_pa;
        name(irat).x_pb = ones(len_ProbB,1) * stereo(X_bl,Y_bl,X_PA,Y_PA,X_PB,Y_PB).x_pb;
        name(irat).y_pb   = ones(len_ProbB,1) * stereo(X_bl,Y_bl,X_PA,Y_PA,X_PB,Y_PB).y_pb;
       
        X_P = [name(irat).x_pa ;name(irat).x_pb]
        Y_P = [name(irat).y_pa ;name(irat).y_pb]
       
        tableau(irat).each = [ tableau(irat).each(:,1:6) X_P tableau(irat).each(:,7:end)]
        tableau(irat).each = [ tableau(irat).each(:,1:7) Y_P tableau(irat).each(:,8:end)]
       
       %%%% to change
       X_ch = zeros(1,len_ProbA+len_ProbB).';
       Y_ch = zeros(1,len_ProbA+len_ProbB).';
       tableau(irat).each = [ tableau(irat).each(:,1:9) X_ch tableau(irat).each(:,10:end)]
       tableau(irat).each = [ tableau(irat).each(:,1:10) Y_ch tableau(irat).each(:,11:end)]
       %%%%
        
        
        
        %positon electrode pour les probe A(1)=>pta  et B(2)=>s1
       % probeA
        X_bl         = Table_experience(irat,4); %commun pour probA et B 
        Z_bl         = Table_experience(irat,6); %commun pour probA et B 
        
        zE01_A       = Table_experience(irat,12);
        angle_A      = Table_experience(irat,13);
        nb_chan_A    = [1:Table_experience(irat,8)]
        inter_chan_A = Table_experience(irat,9)

        name(irat).name_channel_A   =  str2double(Table_A(angle_A, zE01_A, nb_chan_A, inter_chan_A, X_bl, Z_bl).prob_A.name_channel);
        name(irat).rename_channel_A = str2double(Table_A(angle_A, zE01_A, nb_chan_A, inter_chan_A, X_bl, Z_bl).prob_A.rename_channel);
        name(irat).absolute_depth_A = cell2mat(Table_A(angle_A, zE01_A, nb_chan_A, inter_chan_A, X_bl, Z_bl).prob_A.absolute_depth);
        name(irat).structure_A    = cell2mat(Table_A(angle_A, zE01_A, nb_chan_A, inter_chan_A, X_bl, Z_bl).prob_A.structure) ;

        %probeB
        zE01_B = Table_experience(irat,19)
        angle_B = Table_experience(irat,20)
        nb_chan_B = [1:Table_experience(irat,15)]
        inter_chan_B = Table_experience(irat,16)
        
        name(irat).name_channel_B   =  str2double(Table_B(angle_B, zE01_B, nb_chan_B, inter_chan_B, X_bl, Z_bl).prob_B.name_channel);
        name(irat).rename_channel_B = str2double(Table_B(angle_B, zE01_B, nb_chan_B, inter_chan_B, X_bl, Z_bl).prob_B.rename_channel);
        name(irat).absolute_depth_B = cell2mat(Table_B(angle_B, zE01_B, nb_chan_B, inter_chan_B, X_bl, Z_bl).prob_B.absolute_depth);
        name(irat).structure_B    = cell2mat(Table_B(angle_B, zE01_B, nb_chan_B, inter_chan_B, X_bl, Z_bl).prob_B.structure) ;
 
        name(irat).absolute_depth= [name(irat).absolute_depth_A;name(irat).absolute_depth_B];
        name(irat).rename_channel = [name(irat).rename_channel_A; name(irat).rename_channel_B];
        name(irat).structure = [name(irat).structure_A ;name(irat).structure_B];
        name(irat).sub_structure =  zeros(1,maxisize).'% à remplir à la main
        
        tableau(irat).each = [ tableau(irat).each(:,1:11) name(irat).sub_structure tableau(irat).each(:,12:end)];
        
        % NaN resahpe MAXI
        
        tableau(irat).eachshaped= NaN(48,22)
       nbNaN=[]
       for row=1:size(tableau(irat).eachshaped,1)
           ligne= row-length(nbNaN)
           array=49-row
           if tableau(irat).each(ligne,6).' == array
              tableau(irat).eachshaped(row,:) = [ tableau(irat).each(ligne,:)]
           end
           nbNaN= find(isnan(tableau(irat).eachshaped(1:row,6)))
       end
%         if max < 48
%             tableau(irat).eachshaped= NaN(48,22)
%             if tableau(irat).each(1,6) < 48
%                 nbNaN =1
%             else 
%                 nbNaN=[]
%             end
%        for row=1:size(tableau(irat).eachshaped,1)
%            ligne= row-length(nbNaN)
%            array=49-row
%            if tableau(irat).each(ligne,6) == array
%               tableau(irat).eachshaped(row,:) = [ tableau(irat).each(ligne,:)]
%            end
%            nbNaN= find(isnan(tableau(irat).each(1:row,6)))
%        end
%         else
%             tableau(irat).eachshaped=tableau(irat).each
%         end
        
        tableau(irat).eachshaped = [ tableau(irat).eachshaped(:,1:11) name(irat).absolute_depth tableau(irat).eachshaped(:,12:end)];
        tableau(irat).eachshaped = [ tableau(irat).eachshaped(:,1:12) name(irat).rename_channel tableau(irat).eachshaped(:,13:end)];
        tableau(irat).eachshaped = [ tableau(irat).eachshaped(:,1:13) name(irat).structure tableau(irat).eachshaped(:,14:end)];
       



maxidata = tableau(irat).eachshaped
savedir= 'C:\Users\s.carrionfalgarona\Desktop\maxi_tab\'
if ~isnan(tableau(irat).eachshaped(irat,2))
save(fullfile(savedir, [Table_exp{irat,1}, mat2str(tableau(irat).eachshaped(irat,2))]), 'maxidata');
else
save(fullfile(savedir, [Table_exp{irat,1}, mat2str(tableau(irat).eachshaped(irat+1,2))]), 'maxidata');
end



end













