 function [eachtable] = Table_B( angle_B, zE01_B, nb_chan_B, inter_chan_B, X_bl, Z_bl)
 
    L = length(nb_chan_B)*inter_chan_B;

   %% -------- creation d'un tableau de parametre par electrode ------%%
    
   %% electrode de l'experience

        Table_B = {};
        chan_name = [];
        for i=1:length(nb_chan_B)
             varname = [ num2str(i)];
             Table_B{1,i} = varname;
             if length(Table_B{1,i})==1;
                Table_B{1,i} = ['0' Table_B{1,i}];
             end   
        end

    %% associer intervalle à chaque élément de chan_name
    
             depth_mini_B=[]; % les profondeur minimal pour chaque channel
             depth_maxi_B= []; %les profondeur maximal pour chaque channel
        
             for i = 1:length(nb_chan_B)
             Table_B{2,i} = L - (i-1)*inter_chan_B - 100; %depth_min
             Table_B{3,i} = L - (i-1)*inter_chan_B + 99; %depth_max
             depth_mini_B(i)= Table_B{2,i};
             depth_maxi_B(i)= Table_B{3,i};
             end

          depth_mini_B= min(depth_mini_B);% la profondeur la plus superficiel qu'il est possible d'admettre pour la probe
          depth_maxi_B= max(depth_maxi_B);% la profondeur la plus profonde(^^)qu'il est possible d'admettre pour la probe

     

    %% profondeur reel de chaque electrode
       
            for i = 1:length(nb_chan_B)
                Table_B{4,i} = (zE01_B) - (i-1) * inter_chan_B * cosd(angle_B)* cosd(atand(Z_bl/X_bl));% deduire la position de toutes les autres electrodes

            end
            
     %% electrode theorique assih=gné à chaque electrode experimentale
     
            for q = 1:length(nb_chan_B) %assigner la position de references à chaque channels expériementales
                for k = 1:length(nb_chan_B)

                        if (Table_B{4,q}>=Table_B{2,k})&(Table_B{4,q}<=Table_B{3,k})
                            Table_B{5,q} = Table_B{1,k}

                        end

      %% assigner les structures en fonction des profondeurs
      % CC, HPC, NC, Pta, S1, Th
      
                        if  Table_B{4,q} <= 2300.99
                             Table_B{6,q} = 5;
                        end
                        if  (Table_B{4,q} <= 2499.99 & Table_B{4,q} >= 2301) 
                             Table_B{6,q} = 1;
                        end
                        if   Table_B{4,q} >= 2500
                             Table_B{6,q} = 3;
                        end
       %% transformer les profondeurs aberrantes assignées par default en cell vide

                        if  (Table_B{4,q} < depth_mini_B | Table_B{4,q} > depth_maxi_B)
                           Table_B{4,q} = NaN;
                        end

                end
            end
            
%% mettre le tableau dans le bon sens et couper les données aberrantes
        for q = 1:length(nb_chan_B)
                if isnan(Table_B{4,q})
                        Table_B{1,q} = NaN
                        Table_B{2,q} = NaN
                        Table_B{3,q} = NaN
                        Table_B{5,q} = NaN
                        Table_B{6,q} = NaN
                end
            end
        
        Table_B=fliplr(Table_B);% retourne le tableau de droite à gauche



%% utilser les données dans setparam
        eachtable        = struct 
        prob_B           = struct 
        name_channel_B   = struct  %old
        rename_channel_B = struct   %new
        absolute_depth_B = struct   %new
        structure_B      = struct  
           
        name_channel_B   = {Table_B{1,:}} %old
        rename_channel_B = {Table_B{5,:}} %new
        absolute_depth_B = {Table_B{4,:}} %new
        structure_B      = {Table_B{6,:}}
       
        eachtable.prob_B.name_channel   = name_channel_B.';
        eachtable.prob_B.rename_channel = rename_channel_B.';
        eachtable.prob_B.absolute_depth = absolute_depth_B.';
        eachtable.prob_B.structure = structure_B.' ;
   
    end

       
