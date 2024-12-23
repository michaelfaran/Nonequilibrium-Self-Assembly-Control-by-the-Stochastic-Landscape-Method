function[Built_mat_discreteP,Median_mat2,Median_mat_discrete]= Plot_Js_mu_map_many_results_MANY_targets_CALL(RESULTS_PATH,Hallmark_cd,energies,Healing_Driving_value)
% save(horzcat('Built_mat_discreteP',last_part),"Built_mat_discreteP");
% save(horzcat('Median_mat_continuous',last_part),"Median_mat2");
% save(horzcat('Median_mat_discrete',last_part),"Median_mat_discrete");

% clear;
create_me=1;
simplest_continious_time_limit_flag=1;
% Hallmark_cd='C:\Users\admin\Documents\Run_Matlab_Fast_Folder\FASTKMC';

% Rafiq_of_the_many='C:\Users\admin\Documents\RMFO\LoopaUP_1_time_0_85_amp2';
% Rafiq_of_the_many='C:\Users\admin\Documents\RMFO\LoopaUP_N_Control_C_1_T_25_P';
% Rafiq_of_the_many='C:\Users\admin\Documents\RMFO\LoopaUP_Ama_1_t_0_85_A_4_T_25_P';
Rafiq_of_the_many=RESULTS_PATH;
% Rafiq_of_the_many='C:\Users\admin\Documents\GitHub\GradProject\Results\KMC_many2';
% Rafiq_of_the_many='C:\Users\admin\Documents\RMFO\LoopaUP_N_Control_man';
% Rafiq_of_the_many='C:\Users\admin\Documents\RMFO\LoopaUP_N_Control_New_';
% Split the string by backslashes
parts = strsplit(Rafiq_of_the_many, '\');

% Get the last part
last_part = parts{end};
Main_dir_hazirim=horzcat('C:\Users\admin\Documents\RMFO2\',last_part);
if create_me == 1
    % Create the new directory
    mkdir(Main_dir_hazirim);
    disp(['Directory created at: ', Main_dir_hazirim]);
else
    disp('Directory creation skipped.');
end

cd(Main_dir_hazirim);

cd(Rafiq_of_the_many);
%Js_3_5_Mu_0_No_Beast_No_control_W_B_folder_num_0- name of the folder

% energies_true=[2 2.2 2.4 2.5 2.6 2.8 3 3.1 3.2 3.3 3.4 3.5];% 3.2 3.4
% % energies_true=[2 2.2 2.4 2.5 2.6 2.8 3 3.2 3.4 3.5];% 3.2 3.4
% energies_true=[ 3.4 3.5 3.6];% 3.2 3.4
energies_true=energies;

crirtical_num=19999;
min_max_N = find_min_of_max_N_in_subfolders(Rafiq_of_the_many);
turncoat=(min_max_N)*12;

[Targets_num, particles_num] = extract_X_Y(Rafiq_of_the_many);

non_eq_str = extract_non_eq_string(Rafiq_of_the_many);
if isempty(non_eq_str)
non_eq_str='eq';
end

% drive_vec_num=0:1:14;
drive_vec_num=Healing_Driving_value;
     % drive_vec=0:0.2:2.8;
 drive_vec=drive_vec_num;    
turncoat_mat=zeros(length(energies_true),length(drive_vec_num));
topLevelFolder = pwd;
files = dir(topLevelFolder);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags); % A structure with extra info.
% Get only the folder names into a cell array.
subFolderNames = {subFolders(3:end).name};

if create_me==1
    % cd 'C:\Users\admin\Documents\GitHub\GradProject\Results_New\Results_Drives_NEW\'
    % CallME='Tfas_drives';
    topLevelFolder = pwd; % or whatever, such as 'C:\Users\John\Documents\MATLAB\work'
    % Get a list of all files and folders in this folder.
    files = dir(topLevelFolder);
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags); % A structure with extra info.
    % Get only the folder names into a cell array.
    subFolderNames = {subFolders(3:end).name};
    dist_threshold_bottom=0;
    dist_threshold_up=200;
    
    energy_vec=0.0;
    % energies_true=[2 2.2 2.4 2.5 2.6 2.8 3 3.2 3.4 3.5];
    mimi_SS_god=zeros(length(energies_true),length(drive_vec_num),turncoat);
    mimi_SS_god_discrete=zeros(length(energies_true),length(drive_vec_num),turncoat);
     use_time_flag=1;
    % drive_vec=3;
    maxim=10000;%a big but reasonable number
    for kok=1:1:length(energies_true)
    mimi_SS={};
    mimi_SS_discrete={};
    Ap={};
    for ii=1:1:size(subFolderNames,2)
    pp=subFolderNames(1,ii);
    A=horzcat('mainname',num2str(ii),'=','''',pp{1},'''');
    eval ( A );
    Ap=[Ap horzcat('mainname',num2str(ii))];
    end
    Mergi={};
    for ii=1:1:size(subFolderNames,2)
    Mergi=[Mergi  eval(eval('Ap{1,ii}'))];
    end
    
    AD1=horzcat(Rafiq_of_the_many,'\',regexprep(num2str(energies_true(kok), '%5.1f'),'\.','_'));
    AD1 = strrep(AD1, '2_0', '2');
    AD1 = strrep(AD1, '3_0', '3');
    d = dir(AD1);
    % remove all files (isdir property is 0)
    dfolders = d([d(:).isdir]) ;
    % remove '.' and '..' 
    dfolders = dfolders(~ismember({dfolders(:).name},{'.','..'}));
    AD2=horzcat(AD1,'\', dfolders.name);
    
    for joj=1:1:length(Mergi) 
    
    TXT_filepathsgg={};
      %Meet the butcher
       mainname=Mergi{joj}; 
     mainFolder = horzcat(AD1);
     oldpath= addpath([AD1 ...
         '']);
             [~,message,~] = fileattrib([mainFolder,'\*']);
             fprintf('\n There are %i total files & folders.\n',numel(message));
             allExts = cellfun(@(s) s(end-2:end),{message.Name},'uni',0);% Get file ext
             TXTidx = ismember(allExts,'mat');% Search extensions for "CSV" at the end 
             TXT_filepaths = {message(TXTidx).Name};  % Use idx of TXTs to list paths.
             fprintf('There are %i files with *.mat file ext.\n',numel(TXT_filepaths));
      
            TXT_filepathsgg=[TXT_filepathsgg TXT_filepaths];
    end
    TXT_filepaths=TXT_filepathsgg;
    
    addpath ('C:\Users\admin\Documents\GitHub\GradProject\mi');
    addpath('C:\Users\admin\Documents\GitHub\GradProject\github_repo');
    addpath('C:\Users\admin\Documents\GitHub\GradProject\rp-master\rp-master');
     addpath ('C:\Users\admin\Documents\GitHub\GradProject');
     addpath ('C:\Users\admin\Documents\GitHub\GradProject\knee_pt');
    addpath 'C:\Users\admin\Documents\GitHub\GradProject\Beast\Matlab';
    
    %check b4 go, if the date of the simulation results is less then the 18_5
    %flag or not, and all others for resample, smothing etc
    inital_flag=1;
    resample_flag=1;
    %if smooth then change those
    data_smooth_flag=0;
    window_samples=1000000;
    helper=1;
    totalic_length=2*10^7;
    %nevertheless those as well about resample
    downsample_index=1000;
    post_mortem_flag_18_5_21=1;
    
     mimi_SS={};
     mimi_SS_discrete={};
    
    for idi=1:1:length(energy_vec)
    for jj=1:1:length(drive_vec)
    median_vec=zeros(length(energy_vec),length(drive_vec));
    
    
            
            %initate our dreams
            mu = drive_vec(jj);
            energy= energy_vec(idi);
            %num_of_targets_might_change
            %check check
    %3/10/22 here consider to change the mu format in num2str depends
            %on what is needed  num2str(mu) vs. num2str(mu, '%5.1f')
            %on what is needed
            %the figures of merit,"distance","energy","entropy
    if mod(mu,0.5)~=0
            if resample_flag==0
            energy_str =horzcat('energy_UP_vec_mu_', num2str(mu, '%5.1f'),'_energy_',char(sprintfc('%0.1f',energy)));
            entropy_str = horzcat('entropy_UP_vec_mu_', num2str(mu, '%5.1f'),'_energy_',char(sprintfc('%0.1f',energy)));
            distance_str =horzcat('distance_UP_vec_mu_', num2str(mu, '%5.1f'),'_energy_',char(sprintfc('%0.1f',energy)));
            time_str =horzcat('times_vec_mu_', num2str(mu, '%5.1f'),'_energy_',char(sprintfc('%0.1f',energy)));
    
            elseif resample_flag==1
    %         energy_str =horzcat('energy_vec_mu_', num2str(mu, '%5.1f'),'_energy_',char(sprintfc('%0.1f',energy)));
    %         entropy_str = horzcat('entropy_vec_mu_', num2str(mu, '%5.1f'),'_energy_',char(sprintfc('%0.1f',energy)));
    %         distance_str =horzcat('distance_vec_mu_', num2str(mu, '%5.1f'),'_energy_',char(sprintfc('%0.1f',energy)));
                    energy_str =horzcat('energy_UP_vec_mu_', num2str(mu, '%5.1f'),'_energy_',char(sprintfc('%0.1f',energy)));
            entropy_str = horzcat('entropy_UP_vec_mu_', num2str(mu, '%5.1f'),'_energy_',char(sprintfc('%0.1f',energy)));
            distance_str =horzcat('distance_UP_vec_mu_', num2str(mu, '%5.1f'),'_energy_',char(sprintfc('%0.1f',energy)));
            time_str =horzcat('times_vec_mu_', num2str(mu, '%5.1f'),'_energy_',char(sprintfc('%0.1f',energy)));
    
            end
    else
    
            if resample_flag==0
            energy_str =horzcat('energy_UP_vec_mu_', num2str(mu, '%5.1f'),'_energy_',char(sprintfc('%0.1f',energy)));
            entropy_str = horzcat('entropy_UP_vec_mu_', num2str(mu, '%5.1f'),'_energy_',char(sprintfc('%0.1f',energy)));
            distance_str =horzcat('distance_UP_vec_mu_', num2str(mu, '%5.1f'),'_energy_',char(sprintfc('%0.1f',energy)));
            time_str =horzcat('times_vec_mu_', num2str(mu, '%5.1f'),'_energy_',char(sprintfc('%0.1f',energy)));
    
            elseif resample_flag==1
    %         energy_str =horzcat('energy_vec_mu_', num2str(mu, '%5.1f'),'_energy_',char(sprintfc('%0.1f',energy)));
    %         entropy_str = horzcat('entropy_vec_mu_', num2str(mu, '%5.1f'),'_energy_',char(sprintfc('%0.1f',energy)));
    %         distance_str =horzcat('distance_vec_mu_', num2str(mu, '%5.1f'),'_energy_',char(sprintfc('%0.1f',energy)));
                    energy_str =horzcat('energy_UP_vec_mu_', num2str(mu, '%5.1f'),'_energy_',char(sprintfc('%0.1f',energy)));
            entropy_str = horzcat('entropy_UP_vec_mu_', num2str(mu, '%5.1f'),'_energy_',char(sprintfc('%0.1f',energy)));
            distance_str =horzcat('distance_UP_vec_mu_', num2str(mu, '%5.1f'),'_energy_',char(sprintfc('%0.1f',energy)));
            time_str =horzcat('times_vec_mu_', num2str(mu, '%5.1f'),'_energy_',char(sprintfc('%0.1f',energy)));
    
            end
    
    end       
    
    
    
            %distancething
            listt= ~cellfun('isempty',strfind(TXT_filepaths,distance_str));
            address=find (listt==1);
            %energyzing
            listt2= ~cellfun('isempty',strfind(TXT_filepaths,energy_str));
            address2=find (listt2==1);
            mutual_min_vec=zeros(length(address2),1);
            listt3= ~cellfun('isempty',strfind(TXT_filepaths,entropy_str));
            address3=find (listt3==1);   
            listt4= ~cellfun('isempty',strfind(TXT_filepaths,time_str));
            address4=find (listt4==1);  
            if isempty(address4)
                time_str =horzcat('denial_vec_mu_', num2str(mu, '%5.1f'),'_energy_',char(sprintfc('%0.1f',energy)));
                listt4= ~cellfun('isempty',strfind(TXT_filepaths,time_str));
                address4=find (listt4==1);              
            end
            
    
    
    
    A={};
    A_reduced={};
    
    mimi=zeros(1,turncoat);
    mimi_discrete=zeros(1,turncoat);    
        for ii=1:1:turncoat
            pivot2=load(TXT_filepaths{address2(ii)});
            pivot1=load(TXT_filepaths{address(ii)});
            pivot4=load(TXT_filepaths{address4(ii)});
    
            if resample_flag==1
               if post_mortem_flag_18_5_21==1
                Total_energy2=pivot2.foo; 
                else
    
                Total_energy2= sum(cumsum(pivot2.foo),2);
               end
    
                if length(Total_energy2)~=(totalic_length/downsample_index)
                    IL=floor(linspace(1,length(Total_energy2),length(Total_energy2)));
                else
                    IL=floor(linspace(1,length(Total_energy2),length(Total_energy2)));
                end
                Total_energy =Total_energy2(IL);
%             if use_time_flag==1
                time_o=pivot4.foo; 
                time= time_o(:,1);
                blak=sort(time);
                Mining=min(blak(blak~=0));
                time(find (time==0))=Mining;
                Aop=cumsum(time);
                CCp=floor(Aop(end));
                Distance=pivot1.foo;
                qqqqq=pivot1.foo;
                Distance=Distance(IL,:);
                New_time=linspace(0,CCp,length(IL));
        
                Total_energy=interp1(Aop, Total_energy,New_time);
        
                % Determine the number of columns (L) in Distance
                L = size(Distance, 2);
                
                % Preallocate the result matrix SSSS_array for interpolated values
                SSSS_array = zeros(length(New_time), L, 'single');
                
                % Loop through each column of Distance and perform interpolation
                for l = 1:L
                    % Perform interpolation for each column of Distance
                    SSSS_array(:, l) = interp1(Aop, single(Distance(:, l)), New_time, '');
                    
                    % Handle boundary conditions as described
                    SSSS_array(1, l) = SSSS_array(3, l);
                    SSSS_array(2, l) = SSSS_array(3, l);
                end
                
                % Convert interpolated results to int32 and reshape the Distance array
                Distance = int32(SSSS_array);
                
                % Initialize Total_energy
                Total_energy(1) = 0;
                
                % Assign values from qqqqq to the first row of Distance
                Distance(1, :) = qqqqq(1, :);

%                 else
%                 New_time=1:1:length(Total_energy);
%                         qqqqq=pivot1.foo;
%                         Distance=pivot1.foo;
%                 Distance=Distance(IL,:);
%             end
    
                elseif resample_flag==0
                Total_energy=ceil(pivot2.foo);
                end
              %dist_threshold_bottom
              %dist_threshold_top
    

              
              % mimi1_up=New_time(find(Distance(:,1)==dist_threshold_bottom,1));
              % mimi2_up=New_time(find(Distance(:,2)==dist_threshold_bottom,1));          
              % 
              % mimi1_bottom=New_time(find(Distance(:,1)<dist_threshold_up,1));
              % mimi2_bottom=New_time(find(Distance(:,2)<dist_threshold_up,1));
              % 
              % 
              % mimi1= mimi1_up-mimi1_bottom;
              % mimi2= mimi2_up-mimi2_bottom; 

MMM = size(Distance, 1);  % Trajectory length (MMM)
L = size(Distance, 2);    % Number of targets (L)

% Initialize arrays to store results for each target
mimi_upX = zeros(1, L);
mimi_bottomX = zeros(1, L);
mimiX = zeros(1, L);

% Generalized processing for all L targets
for l = 1:L
    % Find mimi_up for each target l
    idx_up = find(Distance(:, l) == dist_threshold_bottom, 1);
    if isempty(idx_up)
        mimi_upX(l) = New_time(end);  % Use fallback if not found
    else
        mimi_upX(l) = New_time(idx_up);
    end

    % Find mimi_bottom for each target l
    idx_bottom = find(Distance(:, l) < dist_threshold_up, 1);
    if isempty(idx_bottom)
        mimi_bottomX(l) = New_time(end);  % Use fallback if not found
    else
        mimi_bottomX(l) = New_time(idx_bottom);
    end

    % Calculate mimi for each target
    mimiX(l) = mimi_upX(l) - mimi_bottomX(l);
end



%               if (simplest_continious_time_limit_flag==1)
% 
%                   if isempty(mimi1)
%                   mimi1=New_time(length(Total_energy));  
%                   end
%                   if isempty(mimi2)
%                   mimi2=New_time(length(Total_energy));
%                   end
%                   mimi(ii)=min(mimi1,mimi2) ;
%                   if mimi(ii)==0
%                      mimi(ii)=eps;
%                   end
% %here is rewriting for discrete time
%                  qqqqq=pivot1.foo;
%                  New_time=1:1:size(qqqqq,1);
%                  Distance=pivot1.foo;
%                  Distance=Distance(IL,:);
%                  mimi1_up=New_time(find(Distance(:,1)==dist_threshold_bottom,1));
%                  mimi2_up=New_time(find(Distance(:,2)==dist_threshold_bottom,1));          
% 
%                  mimi1_bottom=New_time(find(Distance(:,1)<dist_threshold_up,1));
%                  mimi2_bottom=New_time(find(Distance(:,2)<dist_threshold_up,1));
% 
% 
%                  mimi1= mimi1_up-mimi1_bottom;
%                  mimi2= mimi2_up-mimi2_bottom;          
% 
%                  if isempty(mimi1)
%                  mimi1=New_time(size(qqqqq,1));  
%                  end
%                   if isempty(mimi2)
%                   mimi2=New_time(size(qqqqq,1));
%                   end
%                   mimi_discrete(ii)=min(mimi1,mimi2) ;
%                   if mimi_discrete(ii)==0
%                      mimi_discrete(ii)=eps;
%                   end
%               end

% Handling simplest_continious_time_limit_flag
if simplest_continious_time_limit_flag == 1
    % Find the minimum mimi across all targets
    mimi(ii) = min(mimiX);
    if mimi(ii) == 0
        mimi(ii) = eps;  % Avoid having zero
    end

    % Check for empty mimi values
    % for l = 1:L
    %     if isempty(mimi(l))
    %         mimi(l) = New_time(length(Total_energy));  
    %     end
    % end
end

% Discrete time processing
if simplest_continious_time_limit_flag == 1
    qqqqq = pivot1.foo;
    New_time = 1:1:size(qqqqq, 1);
    Distance = pivot1.foo;
    Distance = Distance(IL, :);  % Adjust based on IL
    
    % Re-initialize for discrete time processing
    mimi_up = zeros(1, L);
    mimi_bottom = zeros(1, L);
    mimi_discretex = zeros(1, L);

    % Generalized processing for discrete case for all targets
    for l = 1:L
        idx_up = find(Distance(:, l) == dist_threshold_bottom, 1);
        if isempty(idx_up)
            mimi_up(l) = New_time(end);  % Use fallback if not found
        else
            mimi_up(l) = New_time(idx_up);
        end

        idx_bottom = find(Distance(:, l) < dist_threshold_up, 1);
        if isempty(idx_bottom)
            mimi_bottom(l) = New_time(end);  % Use fallback if not found
        else
            mimi_bottom(l) = New_time(idx_bottom);
        end

        % Calculate mimi_discrete for each target
        mimi_discretex(l) = mimi_up(l) - mimi_bottom(l);
    end

    % Find the minimum across all targets for discrete time
    mimi_discrete(ii) = min(mimi_discretex);
    if mimi_discrete(ii)==crirtical_num
        mimi_discrete(ii)=mimi_discrete(ii)+1;
    end
    % if mimi_discrete(ii) == 0
    %     mimi_discrete(ii) = eps;  % Avoid having zero
    % end
end

%               else
% %26.2.24 here is without the 0.01 to tfas, but 0.01 out of the total time
%                   if isempty(mimi1)
%                   mimi1=New_time(length(Total_energy));  
%                   end
%                   if isempty(mimi2)
%                   mimi2=New_time(length(Total_energy));
%                   end
%                   mimi(ii)=min(mimi1,mimi2) ;
%                   if mimi(ii)==0
%                      mimi(ii)=eps;
%                   end
% 
%                  qqqqq=pivot1.foo;
%                  %New_time=1:1:size(qqqqq,1);
%                  mimi_discrete(ii)=size(qqqqq,1);

                  
%              end
        % end
    

        end
        mimi=mimi*downsample_index*helper;
        mimi_discrete=mimi_discrete*downsample_index*helper;
        mimi_SS=[mimi_SS; mimi];
        mimi_SS_discrete=[mimi_SS_discrete; mimi_discrete];
    end
    mimi_SS_god(kok,:,:)=cell2mat(mimi_SS);
    mimi_SS_god_discrete(kok,:,:)=cell2mat(mimi_SS_discrete);
    end
    end
% Assuming `Targets_num` and `particles_num` are already extracted

% Go to the specified directory
cd(Hallmark_cd);

% Construct the dynamic file names for saving
filename_god = sprintf(horzcat('mimi_SS_god_%d_T_%d_N_',non_eq_str,'.mat'), Targets_num, particles_num);
filename_god_discrete = sprintf(horzcat('mimi_SS_god_discrete_%d_T_%d_N_',non_eq_str,'.mat'), Targets_num, particles_num);
if length(energies)==13
filename_god_discrete = sprintf(horzcat('mimi_SS_god_discrete_%d_T_%d_N13_',non_eq_str,'.mat'), Targets_num, particles_num);
save(filename_god_discrete, 'mimi_SS_god_discrete');
end
% Save the .mat files with the dynamic file names
save(filename_god, 'mimi_SS_god');
save(filename_god_discrete, 'mimi_SS_god_discrete');
end
% Return to the specified directory
cd(Hallmark_cd);

% Load the dynamically named .mat files
load(filename_god);
load(filename_god_discrete);

% You can now proceed with the rest of your code
% create_Tfas_distrubtions(mimi_SS_god, mimi_SS_god_discrete, energies_true);
% In the above, I think there is a bug when all are not built
% continuous_time


Median_mat2=zeros(size(mimi_SS_god,1),length(drive_vec_num));
%Built_mat=zeros(10,15);
Median_mat_discrete=zeros(size(mimi_SS_god,1),length(drive_vec_num));
Built_mat_discrete=zeros(size(mimi_SS_god,1),length(drive_vec_num));
turncoat_vec=zeros(size(mimi_SS_god,1),length(drive_vec_num));
for ii=1:1:size(mimi_SS_god,1)
    for jj=1:1:size(mimi_SS_god,2)
        Median_mat_discrete(ii,jj)=median(mimi_SS_god_discrete(ii,jj,:));
        Built_mat_discrete(ii,jj)=length(find(mimi_SS_god_discrete(ii,jj,:)<2*10^7));
        Median_mat2(ii,jj)=median(mimi_SS_god(ii,jj,:));
        turncoat_vec(ii,jj)=length(find(mimi_SS_god_discrete(ii,jj,:)<30*10^7));
        %Built_mat_discrete(ii,jj)=length(find(mimi_SS_god(ii,jj,:)<2*10^7));
    end
end
addpath 'C:\Users\admin\Documents\GitHub\GradProject\Key Figures\KMC'%'HeatMap_KMC_J_mu.mat'

% cd('C:\Users\admin\Documents\GitHub\GradProject\Results_New');
% load('No_Control_Built_Mat.mat');
cd(Rafiq_of_the_many);
% Heated_debate_map_CALL(Main_dir_hazirim,energies_true,drive_vec,Built_mat_discrete,Median_mat_discrete,Median_mat);
Built_mat_discreteP=Built_mat_discrete./turncoat;
path_str =cd;
cd('C:\Users\admin\Documents\RMFO2');
% Split the string by backslashes
parts = strsplit(path_str, '\');

% Get the last part
last_part = parts{end};
cd(Hallmark_cd);
save(horzcat('Built_mat_discreteP',last_part),"Built_mat_discreteP");
save(horzcat('Median_mat_continuous',last_part),"Median_mat2");
save(horzcat('Median_mat_discrete',last_part),"Median_mat_discrete");
cd(Rafiq_of_the_many);
