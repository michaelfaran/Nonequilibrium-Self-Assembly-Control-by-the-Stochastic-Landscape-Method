%activate in main folder only
%Michael add 10/10/2022
% Adding results instead of the orginial Place where the simulation was
% stuck
clear;

time_Factor=4*10^-7;
create_me=1;
simplest_continious_time_limit_flag=1;
Hallmark_cd=cd;
addpath(Hallmark_cd);
% Rafiq_of_the_many='C:\Users\admin\Documents\RMFO\LoopaUP_1_time_0_85_amp2';
Rafiq_of_the_many=horzcat(Hallmark_cd,'\Data\SI_Fig3\L');
% Rafiq_of_the_many='C:\Users\admin\Documents\GitHub\GradProject\Results\KMC_many2';
% Rafiq_of_the_many='C:\Users\admin\Documents\RMFO\LoopaUP_N_Control_man';
% Rafiq_of_the_many='C:\Users\admin\Documents\RMFO\LoopaUP_N_Control_New_';
% Split the string by backslashes
parts = strsplit(Rafiq_of_the_many, '\');

% Get the last part
last_part = parts{end};
Main_dir_hazirim=horzcat(Hallmark_cd,'\',last_part);
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

energies_true=[2 2.2 2.4 2.5 2.6 2.8 3 3.1 3.2 3.3 3.4 3.5];% 3.2 3.4
energies_true=3.5;
% energies_true=[2 2.2 2.4 2.5 2.6 2.8 3 3.2 3.4 3.5];% 3.2 3.4
crirtical_num=19999;
min_max_N = find_min_of_max_N_in_subfolders(Rafiq_of_the_many);
turncoat=(min_max_N)*12;

[Targets_num, particles_num] = extract_X_Y(Rafiq_of_the_many);

drive_vec_num=0:1:14;
drive_vec_num=0;
     % drive_vec=0:0.2:2.8;
 drive_vec=0;    
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
    
    % addpath ('C:\Users\admin\Documents\GitHub\GradProject\mi');
    % addpath('C:\Users\admin\Documents\GitHub\GradProject\github_repo');
    % addpath('C:\Users\admin\Documents\GitHub\GradProject\rp-master\rp-master');
    %  addpath ('C:\Users\admin\Documents\GitHub\GradProject');
    %  addpath ('C:\Users\admin\Documents\GitHub\GradProject\knee_pt');
    % addpath 'C:\Users\admin\Documents\GitHub\GradProject\Beast\Matlab';
    
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
end

    A={};
    A_reduced={};
    




            figgg=figure; 
  bbb=get(figgg,'Position');
  h_factor=bbb(3)/bbb(4);
%here is achem instead of pnas where is 17.8
  new_width=17.5;
  set(figgg, 'Units', 'centimeters', 'Position',[2 2 new_width new_width]);

%     t = tiledlayout('flow');
    t=tiledlayout(4,3);
    mimi=zeros(1,length(address2));
    %popi=length(address2);
    popi=12;
    for ii=1:1:popi
        pivot1=load(TXT_filepaths{address(ii)});   
        pivot2=load(TXT_filepaths{address2(ii)});
        % Assuming TXT_filepaths{address(ii)} is a valid path string
    original_filepath = TXT_filepaths{address(ii)};
    
    % Split the string to isolate the directory path
    [pathstr, ~, ~] = fileparts(original_filepath);
    
    % Define the desired file name
    new_filename = 'beast_cp_mu_0.0_energy_0.0_run_num_1_total_num_target_2.mat';
    
    % Create the new full file path
    new_filepath = fullfile(pathstr, new_filename);
    

        load(new_filepath);
        CP=foo;

        if resample_flag==1
           if post_mortem_flag_18_5_21==1
            Total_energy2=pivot2.foo; 
            else

            Total_energy2= sum(cumsum(pivot2.foo),2);
           end
downsample_index=1;
        IL=floor(linspace(1,length(Total_energy2),length(Total_energy2)/downsample_index));
        Total_energy =Total_energy2(IL);
        elseif resample_flag==0
        Total_energy=ceil(pivot2.foo);
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
                Time=Aop*time_Factor;
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

%optional way to create figure
    
    nexttile
    tic
    %Total_energy= pivot2.foo; 
    %Total_energy=movingaverage(ceil(pivot2.foo));
%    o=beast(Total_energy, 'start', 0, 'tseg.min',floor(0.1*mimi(ii)),'season','none');
% o=beast(Total_energy, 'start', 0, 'tseg.min',0.01*length(Total_energy2),'season','none');
%    %here take the new Y
%    yyy=figure;xxx=plotbeast(o);
%    h = findobj(xxx(1),'Type','line');
%    Y=h(end-1,:).YData; 
%    close(yyy)
% %    o=beast(Total_energy, 'start', 0,'season','none');
% %    oo=beast(o.trend.order, 'start', 0,'season','none');
% % cp=sort(take_change_points(o.trend.cpOccPr,o.trend.ncp_median,floor(0.1*mimi(ii))));
% cp=sort(o.trend.cp(1:o.trend.ncp_median));
% 
% if sum(isnan(cp))~=0
%    cp=cp(1:(find (isnan(cp),1)-1));     
% end   
% cp=[1 ;cp];
% %taking the first one as well
% %     cp=[1;cp];
%    if ~isempty(find (cp==0))
%    cp=cp(find (cp~=0)); 
%    end     
% 
% if sum(isnan(cp))~=0
%    cp=cp(1:(find (isnan(cp),1)-1));     
% end   
% cp=[1 ;cp];
% %taking the first one as well
% %     cp=[1;cp];
% if ~isempty(find (cp==0))
% cp=cp(find (cp~=0)); 
% end     


cp=CP(:,1);
% noisy_signal = awgn(Total_energy, snr, 'measured');
        plot(Time*1000,Total_energy);
        ylimits=ylim;
        hold on;
        % L=(find(Distance(:,1)==0,1))*10^3;
        % cps=cp(cp<(min(find(Distance(:,2)==0,1),find(Distance(:,1)==0,1)))*fcc);
        % cps=cp;
        ylim([min(Total_energy)-20 0]);
        cps=cp(1:find(cp==0,1)-1);
        ylimits=ylim;

            %plot( [cp(i),cp(i)]*10000, get(gca,'Ylim'),'color',[0.8500 0.3250 0.0980]);
xxxxx=Aop(ceil((cps/1000)))*time_Factor*1000;
s=scatter(xxxxx,(ylimits(1)+7).*ones(1,length(xxxxx)),[],[0.8500 0.3250 0.0980],"|",'filled','MarkerEdgeColor',[0.8500 0.3250 0.0980]);
            s.SizeData = 20;
            %xlim([0 2000]*10000);
            %ylabel('$E$\ $[K_{B}T]$','Interpreter','latex');
             % ylabel( '{\boldmath${E$\ $[K_{B}T]}$}','Interpreter','latex')
if ii==1 || ii==4 || ii==7 || ii==10
         ylabel( '{\boldmath${E\ [k_{B}T]}$}','Interpreter','latex')
end
if ii ==10 || ii ==11 || ii ==12 
            % xx0=xlabel('Time [Cs]');
            xx0= xlabel( '{\boldmath${Time\ [s]}$}','Interpreter','latex');
end
%              xx0.Position=xx0.Position-[0 10 0];
            % ylim([min(Total_energy)-10 0])  
            xlim([0 Time(end)*1000])
               xlim([0 0.01]);
  % xlim([ 0 last_time]);
 % xlim([Time(ceil((cps(end-25)/fcc))) Time(2317)*fcc]);
set(gca,'FontSize',6);
yyaxis right
Trend=CP(1:find(cp==0,1)-1,7);
time_for_trend2=CP(1:find(cp==0,1)-1,3)*time_Factor*1000;
time_for_trend2=time_for_trend2(time_for_trend2>0);
Trend=Trend(time_for_trend2>0);
yaks=find(time_for_trend2<time_for_trend2(end)+1);
time_for_trend=time_for_trend2(yaks);
Trend=Trend(yaks);

times=time_for_trend; trend_values=Trend/time_Factor;time_durations = diff(time_for_trend);  % Difference between consecutive time points
% Initialize time axis and square wave
times=[0; times]/1000;
total_time = times(end);  % The total time span from the first to the last time point
time_axis = times(1):time_Factor:total_time;  % A finer time axis for smooth plotting
square_wave = zeros(1, length(time_axis));  % Initialize the square wave

% Fill the square wave based on trend values and time intervals
current_time_index = 1;
% trend_values(269)=-0.0003/time_Factor;

for i = 1:length(trend_values)-1
    % Find the time indices that correspond to the current time interval
    time_indices = time_axis >= times(i) & time_axis < times(i+1);
    
    % Assign the trend value to the corresponding time interval
    square_wave(time_indices) = trend_values(i);
end
% The last segment: fill in the last trend value until the last time
% square_wave(time_axis >= times(end-1)) = trend_values(end);
% xlim([ 0 last_time+5300*time_Factor]);
plot(time_axis, square_wave,'k'); 
% xlim([ 22500 23440]*time_Factor);
%  ylim([-8 1.5]/time_Factor);
% pp=ylim;
ax=gca;
xgca=gca;
ax.YAxis(2).Color = 'k';  
ax.YAxis(1).Color = [0, 0.4470, 0.7410];  
% set(gca, 'Position',[0.1300 0.2214 0.7338 0.7036]);
if ii==3 || ii==6 || ii==9 || ii==12
 ylabel( '{\boldmath${\ Trend\ [{k_{B}T}\cdot{s^{-1}}]}$}','Interpreter','latex')
end 

% ax = axes(figgg);
% han = gca;
% han.Visible = 'off';
% xlabel(t,'MC steps','FontSize',6)  ;    
% ylabel(t,'$E \ [K_{B}T]$','FontSize',6,'Interpreter','latex');
% set(gca,'FontSize',6);
hold on;
end

% savefig(figure(gcf),fullfile(Main_dir_hazirim,mydir,'Energy Trajectories.fig'),'compact');
% cd('C:\Users\admin\Pictures');
% saveas(figgg,horzcat('EnergyTrajectories' ,'.epsc'));
% saveas(figgg,horzcat('EnergyTrajectories' ,'.png'));
% saveas(figgg,horzcat('EnergyTrajectories' ,'.fig'));
% print ('EnergyTrajectories','-depsc')
cd(hallmark_cd);
print ('EnergyTrajectoriesKMC','-dpng','-r600');
% cd('C:\Users\admin\Documents\Run_Matlab_Fast_Folder\FASTKMC');
% print ('EnergyTrajectoriesKMC','-dpng','-r600');

% pop=median(mimi);
% figure; scatter(mu,mimi(:),60,'filled');
% ylabel('Tfas Distribution','FontSize',36);
% xlabel('Drive Mu','FontSize',36);
% title(horzcat('Tfas Distribution ' ,num2str(mu),'',num2str(energy_vec), ' Number of Targets ',str2),'FontSize',36);
% hold on;
% scatter(mu,pop,90,'d','k','filled');
ssz=size(A);
% cmin=length(Total_energy);
% cmax=0;
% for ii=1:1:ssz(1)
% c=min(A{ii,5});
% m=max(A{ii,5});
% cmin=min(c,cmin);
% cmax=max(m,cmax);
% end
% hold off;figure
% for ii=1:1:ssz(1)
%     x=A{ii,1};%mean_vec
%     y=A{ii,2};%std_vec
%     z=A{ii,4};%trend
%     c=A{ii,5};%time
%     
% TrajectoryPlotVecs(x,y,z,c,cmin,cmax);
% 
% end
% xlabel('Mean');
% ylabel('Variance');
% zlabel('Trend');
% title('Trajectories of Simulation CG States with Time Colorbar')
% % savefig('Stochastic 1');
% %mydir='yourfullyqualifiedpath';
% savefig(figure(gcf),fullfile(Main_dir_hazirim,mydir,'Stochastic 1.fig'),'compact');
% 
% %Now the remaining time to self assembly
% 
% cmin=length(Total_energy);
% cmax=0;
% for ii=1:1:ssz(1)
% c=min(A{ii,7}~=0);
% m=max(A{ii,7});
% cmin=min(c,cmin);
% cmax=max(m,cmax);
% end
% hold off;figure
% for ii=1:1:ssz(1)
%     x=A{ii,1};%mean_vec
%     y=A{ii,2};%std_vec
%     z=A{ii,4};%trend
%     c=A{ii,7};%time
%     %new index
%     INDEX_nan=find(isnan (c));
%     Index_Zero=find(c==0);
%     Index_correct=find(c~=0);
% TrajectoryPlotVecs4TfasMap(x,y,z,c,cmin,cmax,INDEX_nan,Index_Zero,Index_correct);
% 
% end
% xlabel('Mean');
% ylabel('Variance');
% zlabel('Trend');
% title('Trajectories of Simulation CG States with Time Remaining to Self Assemble Colorbar')
% % savefig('Stochastic 2');
% savefig(figure(gcf),fullfile(Main_dir_hazirim,mydir,'Stochastic 2.fig'),'compact');
% 
% %Time remanining to Self Assemble, just first one
% 
% cmin=length(Total_energy);
% cmax=0;
% 
% for ii=1:1:ssz(1)
%   AA=A{ii,7};  
%  ending_theme=find(AA==0,1);
%  AA=AA(1:ending_theme);
% c=min(AA(AA>0));
% m=max(AA);
% cmin=min(c,cmin);
% cmax=max(m,cmax);
% 
% end
% 
% hold off;figure
% for ii=1:1:ssz(1)
%     x=A{ii,1};%mean_vec
%     y=A{ii,2};%std_vec
%     z=A{ii,4};%trend
%     c=A{ii,7};%time
%     ending_theme=find(c==0,1);
%     x=x(1:ending_theme);
%     y=y(1:ending_theme);
%     z=z(1:ending_theme);
%     c=c(1:ending_theme);
% %     A_reduced=[A_reduced; x y z];
%     %new index
%     INDEX_nan=find(isnan (c));
%     Index_Zero=find(c==0);
%     Index_correct=find(c~=0);
% % TrajectoryPlotVecs4TfasMap(x,y,z,c,cmin,cmax,INDEX_nan,Index_Zero,Index_correct);
% TrajectoryPlotVecs(x,y,z,c,cmin,cmax);
% end
% xlabel('Mean','FontSize',28);
% ylabel('Variance','FontSize',28);
% zlabel('Trend','FontSize',28);
% title('Trajectories of Simulation CG States with Time Remaining to FIRST Self Assemble Colorbar')
% % savefig('Stochastic 3');
% savefig(figure(gcf),fullfile(Main_dir_hazirim,mydir,'Stochastic 3.fig'),'compact');
% 
% hold off;figure
% 
% for ii=1:1:ssz(1)
%     x=A{ii,1};%mean_vec
%     y=A{ii,2};%std_vec
%     z=A{ii,4};%trend
%     c=A{ii,6};%SA
%     
% TrajectoryPlotVecs(x,y,z,c,0,1);
% 
% end
% xlabel('Mean');
% ylabel('Variance');
% zlabel('Trend');
% title('Trajectories of Simulation CG States with Time Colorbar')
% % savefig('Stochastic 4');
% savefig(figure(gcf),fullfile(Main_dir_hazirim,mydir,'Stochastic 4.fig'),'compact');
% 
% 
% 
% hold off;figure
% for ii=1:1:ssz(1)
%     x=A{ii,1};%mean_vec
%     y=A{ii,2};%std_vec
%     z=A{ii,4};%trend
%     c=A{ii,6};%SA
%     ending_theme=find(c==1,1);
%     x=x(1:ending_theme);
%     y=y(1:ending_theme);
%     z=z(1:ending_theme);
%     c=c(1:ending_theme);
% 
% TrajectoryPlotVecs(x,y,z,c,0,1);
% 
% end
% xlabel('Mean');
% ylabel('Variance');
% zlabel('Trend');
% title('Trajectories of Simulation CG States with Time Colorbar,TFAS ONLY')
% savefig(figure(gcf),fullfile(Main_dir_hazirim,mydir,'Stochastic 5.fig'),'compact');
% 
% cd(xxy);
% A_reduced=[];
% C_reduced=[];
% for ii=1:1:ssz(1)
%     x=A{ii,1};%mean_vec
%     y=A{ii,2};%std_vec
%     z=A{ii,5};%trend
%     w=A{ii,3};%skew
%     v=A{ii,4};%total_trajectorytime
%     c=A{ii,7};%time to self assemble
%     ending_theme=find(c==0,1);
%     x=x(1:ending_theme);
%     y=y(1:ending_theme);
%     z=z(1:ending_theme);
%     c=c(1:ending_theme);
%     w=w(1:ending_theme);
%     v=v(1:ending_theme)
%     B_reduced=[ x; y; z; c]';
%     D_reduced=[ x; y; z; w; v; c]';
%     A_reduced=[A_reduced; B_reduced];
%     C_reduced=[C_reduced; D_reduced];
% end
% save('A_redcued.mat','A_reduced');
% save('C_redcued.mat','A_reduced');
% 
%     M_reduced=C_reduced(:,1:3);
% Mm_reduced=(M_reduced-mean(M_reduced,1))./std(M_reduced,1);
% [coeff,score,latent] = pca(Mm_reduced);
% cmin=0;
% cmax=3000;
% hold off;figure
% for ii=1:1:ssz(1)
%     x=score(:,1);%mean_vec
%     y=score(:,2);%std_vec
%     z=score(:,3);%trend
%     c=C_reduced(:,6);
% %     c=A{ii,7};%time
% %     ending_theme=find(c==0,1);
% %     x=x(1:ending_theme);
% %     y=y(1:ending_theme);
% %     z=z(1:ending_theme);
% %     c=c(1:ending_theme);
% %     A_reduced=[A_reduced; x y z];
%     %new index
% %     INDEX_nan=find(isnan (c));
% %     Index_Zero=find(c==0);
% %     Index_correct=find(c~=0);
% % TrajectoryPlotVecs4TfasMap(x,y,z,c,cmin,cmax,INDEX_nan,Index_Zero,Index_correct);
% TrajectoryPlotVecs(x,y,z,c,cmin,cmax);
% end
% xlabel('Mean','FontSize',28);
% ylabel('Variance','FontSize',28);
% zlabel('Trend','FontSize',28);
% title('NEW COORD,Trajectories of Simulation CG States with Time Remaining to FIRST Self Assemble Colorbar')
% % savefig('Stochastic 3');
% savefig(figure(gcf),fullfile(Main_dir_hazirim,mydir,'Stochastic 6.fig'),'compact');
% 

% %check matrix 
%  load('MetricMap.mat'); %should be h here
% 
%  for ii=1:1:size(A_reduced,1)
% % (ii)
% % 
%  end


