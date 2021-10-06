  %% codé par Pierre et Sofia

    %ceci est un programme pour redefinir la position des channels experiementaux en fonction de la position d'une probe de reference 
    %les variables de références constituent les parametre de la probe utilisée
    %les variables de l'expérences relatent de la position de la probe durant
    %l'expérience

    function [eachtable] = Table_A( angle_A, zE01_A, nb_chan_A, inter_chan_A, X_bl, Z_bl)
    
    L = length(nb_chan_A)*inter_chan_A;

   %% -------- creation d'un tableau de parametre par electrode ------%%
    
   %% electrode de l'experience

        Table_A = {};
        chan_name = [];
        for i=1:length(nb_chan_A)
             varname = [ num2str(i)];
             Table_A{1,i} = [varname];
             if length(Table_A{1,i})==1;
                Table_A{1,i} = ['0' Table_A{1,i}];
             end   
        end

    %% associer intervalle à chaque élément de chan_name
    
             depth_mini_A=[]; % les profondeur minimal pour chaque channel
             depth_maxi_A= []; %les profondeur maximal pour chaque channel
        
             for i = 1:length(nb_chan_A)
             Table_A{2,i} = L - (i-1)*inter_chan_A - 100; %depth_min
             Table_A{3,i} = L - (i-1)*inter_chan_A + 99; %depth_max
             depth_mini_A(i)= Table_A{2,i};
             depth_maxi_A(i)= Table_A{3,i};
             end

          depth_mini_A= min(depth_mini_A);% la profondeur la plus superficiel qu'il est possible d'admettre pour la probe
          depth_maxi_A= max(depth_maxi_A);% la profondeur la plus profonde(^^)qu'il est possible d'admettre pour la probe

     

    %% profondeur reel de chaque electrode
       
            for i = 1:length(nb_chan_A)
                Table_A{4,i} = (zE01_A) - (i-1) * inter_chan_A * cosd(angle_A)* cosd(atand(Z_bl/X_bl));% deduire la position de toutes les autres electrodes

            end
            
     %% electrode theorique assih=gné à chaque electrode experimentale
     
            for q = 1:length(nb_chan_A) %assigner la position de references à chaque channels expériementales
                for k = 1:length(nb_chan_A)

                        if (Table_A{4,q}>=Table_A{2,k})&&(Table_A{4,q}<=Table_A{3,k})
                            Table_A{5,q} = Table_A{1,k}

                        end
               
          

      %% assigner les structures en fonction des profondeurs
      % CC, HPC, NC, Pta, S1, Th
      
                        if  Table_A{4,q} <= 1500.99
                             Table_A{6,q} = 4;
                        end
                        if  (Table_A{4,q} <= 1649.99 & Table_A{4,q} >= 1501) 
                             Table_A{6,q} = 1;
                        end
                        if  (Table_A{4,q} <= 3600.99 & Table_A{4,q} >= 1650) 
                             Table_A{6,q} = 2;
                        end
                        if   Table_A{4,q} >= 3601
                             Table_A{6,q} = 6;
                        end

       %% transformer les profondeurs aberrantes assignées par default en cell vide

                        if  (Table_A{4,q} < depth_mini_A | Table_A{4,q} > depth_maxi_A)
                           Table_A{4,q} = NaN;
                        end
                end
            end
            
%% mettre le tableau dans le bon sens et couper les données aberrantes

            for q = 1:length(nb_chan_A)
                if isnan(Table_A{4,q})
                        Table_A{1,q} = NaN
                        Table_A{2,q} = NaN
                        Table_A{3,q} = NaN
                        Table_A{5,q} = NaN
                        Table_A{6,q} = NaN
                end
            end
        

        Table_A=fliplr(Table_A)% retourne le tableau de droite à gauche
        
 

%% utilser les données dans setparam
        eachtable        = struct 
        prob_A           = struct 
        name_channel_A   = struct  %old
        rename_channel_A = struct   %new
        absolute_depth_A = struct   %new
        structure_A      = struct  
           
        name_channel_A   = {Table_A{1,:}} %old
        rename_channel_A = {Table_A{5,:}} %new
        absolute_depth_A = {Table_A{4,:}} %new
        structure_A      = {Table_A{6,:}}
       
        eachtable.prob_A.name_channel   = name_channel_A.';
        eachtable.prob_A.rename_channel = rename_channel_A.';
        eachtable.prob_A.absolute_depth = absolute_depth_A.';
        eachtable.prob_A.structure = structure_A.' ;
   
    end

       
