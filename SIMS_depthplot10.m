%Depth profiling analysis program to enable you to plot and compare depth
%profiles of different elements across different experiments, designed for SIMS analysis. 
%A broad range of normalisation options for different elements and
%combinations of different elements across one or multiple experiments.
%Shift correction enables you to automatically align depth profiles using a
%maximum peak in each profile, or a rising/decaying slope within the
%profile for a specific element. This enables you to align material interfaces and compare
%compositions across these interfaces. 

%%%%%%%%%%%%%To Do
%double axis plots, tidy up shift correct
%clean up plot refining parameters at end, and the individual line naming
%per plot per key, needs to be consistent....

%% Starting Parameters
% Keeping it clean and preventing mixing variable between experiments
clear all 
clc

%adding colorblind colourmaps to path and vector of markers for plotting 
addpath('/Users/horatiocox/Library/CloudStorage/OneDrive-UniversityCollegeLondon/B-PHD/MATLAB/Depth_profile_analysis/Colormaps')
markers = {'+','o','*','x','s','d','^','v','>','<','p','h','+','o','*','x','s','d','^','v','>','<','p','h','+','o','*','x','s','d','^','v','>','<','p','h'};

%data smoothing, set value is the filtering window (higher -> smoother data), set either to zero to turn off
mean_filt = 1;   %mean filter
med_filt = 7; %median filter


plot_mode = 3;
%plot_mode = 0: generate a separate plot for each selected pelement comparing across all the imported files
%plot_mode = 1: full experiments separately - separate plot for each file (depth profile) with all of the elements plotted
%plot_mode = 2: partial experiments separately - separate plot for each file (depth profile) with only pelements plotted
%plot_mode = 3: compares selected elements (pelement) across a number of experiments (pelement), combining all of the files (instances of i) into one plot
%plot_mode = 4: dual axis plot, select and split elements using pelement_double

%note - if you want to simultaniously plot normalised and un-normalised
%elements, you have to use plot mode 4 and split the elements using
%pelement_double, otherwise, if normalisation is on, all the pelements will
%be normalised


fig_save = 0; %set to 1 to save a .tif file in the current folder for each plot, set to zero to turn off

plotnorm = 0; %Set as 0 for plots to normalise to an element within themself, set as a number for them to normalise to a single "reference" plot in the selection
norm = 1;  
%norm = 0 is off (no normalisation)
%norm = 1 & plotnorm = 0 : normalise pelement (within each plot) to nelement (within the same plot)
%norm = 1 & plotnorm = x: normalise pelement (within each plot) to nelement (in plot x)
%norm = 2 & plotnorm = x: compare ratio of pelement/nelement between plots -> pelement(x)/nelement(x) // pelement(plotnorm)/nelement(plotnorm) 

nelement = [1]; %Change the element being used for normalisation(denominator), for multiple elements it will normalise to the sum of the nelements, Requires each element to be in the same colum for each file/experiment
pelement = [3 9 10 11]; %Change the plotted elements and elements being normalised(numerator)


%Plotting with separate left and right y axis scales
pelement_double = [[3 9];[4 4]]; %split plotted elements as follows: [[pelements-yaxis_left];[pelements-yaxis_right]]
%NOTE pelement_double are not fed into normalisation automatically
%for normalisation, put elements into pelement. 
%This way, normalised and unormalised elements can be plotted on the same graph


%transform data scaling
transform_scale = 0; % = 0 turns off all scaling, =1 and fractional difference in y from normalisation is converted to percentage,
% >0 to scale the x data from sputtering time to depth (e.g: if 173s sputtering time = 25nm depth then transform_scale=25/173)

%shift correction to offset and align different plots in the x direction
shift = 0; % aligns plots to the below element peak (Mo) so that all the peaks line up
Mo = 11; %Mo position for peak aligning

shift2 = 0;%0 is off, 1 is on and aligns the plots aligning the decay slope of Si signal or selected element in "silicon"
silicon = 8; %the selected element used to align the plots
scaling = 1; %dummy variable - don't change from 1
scale_factor = 10; % amount to scale the data for shift correction resolution improvement
range_multiplier = 10; % range of data over which to align for shift correction




%% Import the data from txt files
%this opens a search window from which you can select a number of data
%files to analyse and plot
files = uigetfile('*.*',  'All Files (*.*)','MultiSelect','on');     % Select a number of files and it putes them into a files array
if  ischar(files)== 1                                                   %If selecting only one file, makes sure it is cell type for later commands
    files = cellstr(files);
end

for i = 1:numel(files)
    
    
    filename = files{i};
    fileID = fopen(filename,'r');
    if fileID ~= -1                 %If file doesn't exist ID = -1 and following code isn't computed
        
        AA = textscan(fileID,'%s');     %Loads file into a vector
        n = AA{1,1};
        [row col] = find(contains(n,'#'));  %Creates an array with the cell collumns and rows in which the # symbol is found
        Q{i} = AA{1,1}(row(2)+1:row(3)-1,:);  %Selects the data inbetween # which gives only the element symbols
        StartRow = numel(row);                      %Determines start row for numbers as first without a #
        B{i} = dlmread(filename,'',StartRow,0);       %Reads numbers into a cell array
        
        fclose(fileID)
    end
end
%% Smoothing Functions
%B remains untouched, so that multiple filtering runs can be done without
%accumulating changes to B, Only M is changed each time. 
%Can apply neither, one or both filters. 
M = B; %keep B constant and modify M to prevent double smoothing when running multiple times. 

if med_filt ~= 0
    for i = 1:numel(files)
        M{i}(4:end,:) = medfilt1(B{i}(4:end,:),med_filt);
        if mean_filt ~= 0
            M{i}(:,2:end) = smoothdata(M{i}(:,2:end),'movmean',mean_filt); %(2:end)->not smoothing the sputtering time in colum 1
        end
    end
else
    if mean_filt ~= 0
        for i = 1:numel(files)
            M{i}(:,2:end) = smoothdata(B{i}(:,2:end),'movmean',mean_filt);
        end

    else
        M = B; %if both filters are off, this sets M to the raw data B
    end

end
%% Plot Shifting to Align Peak Signal
%section to shift plots to align them (assumption of a raised surface feature)
%This aligns x axis such that the max peaks in y overlap thus aligning the
%depth profiles such that the interface is aligned
% Note: this has to appear before any normalisation, so that the
%normalisation is applied to the "correctly" aligned plots
if shift == true
    
    ll = size(M{1})
    ini = 50
    for i = 1:numel(files)
        M{i} = [ones(ini,ll(2));M{i}(:,:);ones(ini,ll(2))];
        Y{i}(:,Mo+1) = smoothdata(M{i}(:,Mo+1),'movmean',25); %(46 for SIMS 13)(28 for sims 5 pt)`(38 for mo sims5) (50 for sims6 Ti)
        
    end
    for i = 1:numel(files)
        maxindex = find(Y{plotnorm}(:,Mo+1) == max(Y{plotnorm}(:,Mo+1)));
        maxindex2 = find(Y{i}(:,Mo+1) == max(Y{i}(:,Mo+1)));
        ran = maxindex-maxindex2;
        M{i}(ini+ran:end+ran-ini,:) = M{i}(ini:end-ini,:);
        M{i}(:,1) = M{i}(:,1)+(ran * (M{1}(101,1)- M{1}(100,1)));
    end
end


%% Plot Shifting to Signal Decay
%Aligns the decaying signal from one layer to another to align interfaces
%in the depth profile. 
if shift2 == true
    sil=silicon+1;
    ll = size(M{1});
    ini =150;
    scaling = scale_factor;
    for i = 1:numel(files)       
        M{i} = interp1(1:length(M{i}), M{i}, 1:1/scaling:length(M{i}), 'linear');
        M{i} = [ones(ini,ll(2));M{i}(:,:);ones(ini,ll(2))];
        L{i}(:,:) = M{i}(:,:);
        L{i}(:,sil) = smoothdata(M{i}(:,sil),'movmean',50); %(50 for sims 6 Ti)
    end
    
    
    a = L{1}(:,1);
    b = L{1}(:,sil);
    
    [M1 I ] = max(b);
    ST = a(I);
    I = I+(scaling*120); % (sims 5 Mo *80 interp range 10) ( sims 5 pt *10 interp,scale 10 range 4 ) (sims 6 ti interp *90 range 30)(*50 for tiox range 30)

    rang = range_multiplier*scaling;

    c = L{1}(I:I+rang,1);
    d = L{1}(I:I+rang,sil);
    
    figure;
    plot(L{1}(I-rang/2:I+rang*1.5,1),L{1}(I-rang/2:I+rang*1.5,sil),'LineWidth',2)
    hold on
    for j = 1:numel(files)-1
        for k = 1:(2*rang)
            sumv =  (L{1}(I:I+rang,sil)-L{j+1}(I+k-rang:I+k,sil))'*(L{1}(I:I+rang,sil)-L{j+1}(I+k-rang:I+k,sil)); %+(L{1}(I:I+rang,1)-L{j+1}(I+k-rang/2:I+k+rang/2,1))'*(L{1}(I:I+rang,1)-L{j+1}(I+k-rang/2:I+k+rang/2,1));

            dvec{j}(k) = sumv;

        end
        [M2(j) I2(j)] = min(dvec{j});
        I2(j)=I2(j)-rang;

        L{j+1}(:,1) = L{j+1}(:,1)-(I2(j) * (L{j+1}(ini+10,1)- L{j+1}(ini+9,1)));
        M{j+1}(:,1) = L{j+1}(:,1); %this and line above reduce timestamp by integer number of timesteps to shift curve
        M{j+1}(ini-I2(j):end-ini -I2(j),:) = M{j+1}(ini:end-ini,:); % this maps the data back so that it lines up for normalisation again (not needed for pristine as it isn't shifted)

        M{j+1}([1:ini-I2(j) (end-ini -I2(j)+1):end],:) = 1;
        M{j+1}(1:ini-I2(j),1) = M{j+1}(ini-I2(j)+1,1)-2; %takes starting values and moves timestamp 2s left of first datapoint (to move the data padding out the way)
        M{j+1}((end-ini -I2(j)+1):end,1) = M{j+1}((end-ini -I2(j)),1)+2; %takes end values and moves timestamp 2s right of first datapoint



        I3(j) = I2(j)+I;

        plot(L{j+1}(I3(j)-2*rang:I3(j)+rang*3,1),L{j+1}(I3(j)-2*rang:I3(j)+rang*3,sil),'--','LineWidth',1,'DisplayName', strcat(char(Q{j}(sil-1)),files{j}(1:6),'-',files{j}(8:end-4)))

    end
    cmap = plasma2(numel(files));
    legend('show','Location','northeastoutside')
    xline(L{1}(I,1),'-',{'Lower Reference','threshold'});
    xline(L{1}(I+rang,1),'-',{'Upper Reference','threshold'});
    xline(L{1}(I-rang,1),'-',{'Lower Shift','threshold'});
    xline(L{1}(I+rang*2,1),'-',{'Upper Shift','threshold'});

    set(gca,'color','none','ColorOrder', cmap,'FontSize',20,'FontName', 'Arial');
    hold off
    
   
    M{1}(1:ini,1) = M{j+1}(ini+1,1)-2; %for j = -1 takes starting values and moves timestamp 2s left of first datapoint
    M{1}((end-ini+1):end,1) = M{j+1}((end-ini),1)+2; %for j = -1 takes end values and moves timestamp 2s left of first datapoint
        
    
disp(strcat('The Shift Vector Is',num2str(I2)))
I2 = [0,I2]
end

%% Normalisation of Data

%Section to normalise the data
if norm ~= 0
    R = M; %an attempt to make R constant even when M is being normalised so that it continues to normalize as we move past the plot containing the element we want to normalise to

    x=size(M{1}); %this loop finds the smallest array within M and ensures that the normalisation is limited to those dimensions
    for k = 2:length(M) %which prevents errors if the input data have different sizes
        if length(M{k})< x(1)
            x(1) = length(M{k});
        end
        if width(M{k})< x(2)
            x(2) = width(M{k});
        end
    end

    for i = 1:numel(files)


        if length(nelement)>1 %normalise to sum of multiple elements if multiple are inputted
            for j = 2:length(nelement)
                R{i}(:,nelement(1)+1) = R{i}(:,nelement(1)+1) + R{i}(:,nelement(j)+1) 
            end
        end

        %Constant normalisation reference that doesn't change throughout the normalisation (as though M is changed, R is not)
        %normalise either within each depth profile or to a single reference plot 
        if plotnorm == 0|norm == 2 %normalise within a plot (applies to norm == 2 as for this, you normalise within a plot first then between plots)
            normel = R{i}(1:x(1),nelement(1)+1); 
        else %normalising everything to a specific file
            normel = R{plotnorm}(1:x(1),nelement(1)+1); %only normalising to one of the files (one of the depth_profile runs)
        end

        %normalisation of data
        M{i}(1:x(1),pelement+1) = (M{i}(1:x(1),pelement+1))./normel;
        M{i}(isinf(M{i})|isnan(M{i})) = 1.0;


    end

    if norm == 2 %normalising ratio of elements within plots, between plots
        if plotnorm == 0
            disp("error, you need to set a reference plot to normalise to")
        else
            O = M; %prevent the normalisation updating during the loop
            normel2 = O{plotnorm}(:,pelement+1); %now it is M not R as needs to be already normalised within plot
            for i = 1:numel(files) %this normalises each plot to the reference plot
                M{i}(1:x(1),pelement+1) = (M{i}(1:x(1),pelement+1))./normel2;
            end
        end
    end

elseif norm > 2
    disp("Norm error, must be set to 1, 2 or 3")

end


%% Transform Data Scale
if transform_scale ~= 0
    for i = 1:numel(files)
        M{i}(:,1) = M{i}(:,1)*(transform_scale); %change sputtering time to depth in sample
        M{i}(:,pelement+1) = (M{i}(:,pelement+1)-1)*100; %changes normalised data to show percentage increase/decrease

    end
end

%% 

% M{1}(:,1) = M{1}(:,1)-278 
% M{2}(:,1) = M{2}(:,1)-342 
% M{1}(:,3) = M{1}(:,3)*(1/1.77)
% M{2}(:,3) = M{2}(:,3)*(1/1.56)



%% Plotting
if plot_mode == 0 %All files compared in separate plot for each element so pelement needs to be incremented first
    for j = pelement
        figure();
        for i = 1:numel(files)
            x=j;
            plot(M{i}(:,1),M{i}(:,j+1),strcat('-',markers{j}),'MarkerSize',5,'LineWidth',1,'MarkerSize',3,'MarkerIndices',1:1:length(M{i}(:,1)),'DisplayName', strcat(char(Q{i}(j)),files{i}(1:6),'-',files{i}(8:end-4)))
            xlabel('Sputter time')
            ylabel('Counts')
            hold on
        end

        label_out(transform_scale) %function to label x and y depending on the axis transformation mode
        title(Q{1}(j))
        legend('show','Location','NorthEastOutside')
        set(gca,'color','none','ColorOrder', plasma2(numel(files)),'FontSize',20,'FontName', 'Arial');
        hold off

        if fig_save == 1
            figname = [char(Q{i}(j)),'.tif'];
            print('-dtiff',figname)
        end
    end

else %all other plot modes increment the files first
    figure();
    color_vector = plasma2(numel(files));
            
    for i = 1:numel(files)

         if plot_mode == 3 %combining all of the files (instances of i) into one plot

            if shift2 == 0
                for j = pelement
                    x=j;
                    % yyaxis left
                    plot(M{i}(:,1),M{i}(:,j+1),strcat('-','d'),'MarkerFaceColor','w','MarkerSize',4,'LineWidth',3,'MarkerIndices',1:scaling:length(M{i}(:,1)),'DisplayName', strcat(char(Q{i}(j)),files{i}(1:6),'-',files{i}(8:end-4)),'Color',color_vector(i,:))
                    hold on
                end
            else
                for j = pelement
                    x=j;
                    % yyaxis left
                    plot(M{i}(:,1),M{i}(:,j+1),strcat('-','d'),'MarkerFaceColor','w','MarkerSize',3.5,'LineWidth',2,'MarkerIndices',(1+ini-I2(i)):scaling:length(M{i}(:,1)),'DisplayName', strcat(char(Q{i}(j)),files{i}(1:6),'-',files{i}(8:end-4)),'Color',color_vector(i,:))
                    hold on
                end
            end

         elseif plot_mode == 4 %double axis to split elements
                count_dummy = 1;
                cmap_dual = plasma2(2);
                colororder([cmap_dual(1,:);cmap_dual(2,:)])
                yyaxis left
                for j = pelement_double(1,:)                 
                  plot(M{i}(:,1),M{i}(:,j+1),strcat('-',markers{count_dummy}),'Color',cmap_dual(1,:),'MarkerFaceColor','w','MarkerSize',3,'LineWidth',3,'MarkerIndices',1:scaling:length(M{i}(:,1)),'DisplayName', strcat(char(Q{i}(j)),files{i}(1:6),'-',files{i}(8:end-4)))
                    hold on
                    count_dummy = count_dummy +1;
                end

                yyaxis right
                for j =  pelement_double(2,:)
                    plot(M{i}(:,1),M{i}(:,j+1),strcat('-',markers{count_dummy}),'Color',cmap_dual(2,:),'MarkerFaceColor','w','MarkerSize',3,'LineWidth',3,'MarkerIndices',1:scaling:length(M{i}(:,1)),'DisplayName', strcat(char(Q{i}(j)),files{i}(1:6),'-',files{i}(8:end-4)))
                    hold on
                    count_dummy = count_dummy +1;
                end


        else  %plotting each file (instance of i) in a separate plot
            if plot_mode == 1
                for j = 1:numel(Q{i})
                    semilogy(M{i}(:,1),M{i}(:,j+1),strcat('-',markers{j}),'MarkerSize',5,'LineWidth',1,'MarkerSize',3,'MarkerIndices',1:1:length(M{i}(:,1)),'DisplayName',char(Q{i}(j)))
                    hold on
                end
            elseif plot_mode == 2
                for j = pelement
                    semilogy(M{i}(:,1),M{i}(:,j+1),strcat('-',markers{j}),'MarkerSize',5,'LineWidth',2,'MarkerSize',3,'MarkerIndices',1:1:length(M{i}(:,1)),'DisplayName',char(Q{i}(j)))
                    hold on
                end
            end
            hold off

        end

        label_out(transform_scale) %function to label x and y depending on the axis transformation mode
        legend('show','Location','NorthEastOutside')
        
        if plot_mode ~= 4|3
        set(gca,'LineWidth',1.1,'color','none','ColorOrder', plasma2(numel(files)),'FontSize',25,'FontName', 'Arial');
        else
           set(gca,'LineWidth',1.1,'FontSize',25,'FontName', 'Arial'); 
        end

        if norm == 2 & plotnorm ~= 0
            title(strcat(char(Q{i}(pelement)), '/', char(Q{i}(nelement)),'/','pristine'))
        elseif norm == 1 & plotnorm ~= 0
            title(strcat(char(Q{i}(pelement)), '/(', char(Q{i}(nelement)),'_plot',num2str(plotnorm)))
            title(strcat(files{i}(1:end-4), '-', files{i}(8)))
                    title(files{i}(1:end-4))
        elseif norm == 1 & plotnorm == 0
            title(strcat(char(Q{i}(pelement)), '/(', char(Q{i}(nelement))))
        elseif norm == 0 
            title('No Normalisation Applied')
        end


        %                 print('-dtiff',figname);
        mint = 0;
        for i = 1:numel(files)
            if mint>min(M{i}(:,1));
                mint=min(M{i}(:,1));
            end
        end
        xlim([mint max(M{1}(:,1))]); %Have changed from "si" to M{1}(end,1) as sputter time is greater than the element number so were cutting off ened of graph



        if fig_save == 1
            figname = [files{i}(1:end-4),'.tif'];
            print('-dtiff',figname)
        end

    end
hold off
end



%  a = (M{1}(20:80,5)-M{2}(20:80,5))./M{2}(20:80,5)
%  mean(a)
%  std(a)
%% 


    %% 
    
    % V = [M;Q]
    % plot (M{1,1}(:,1),M{1,1}(:,2))
    % filename = 'a26293_9.txt';
    % P{1} = dlmread(filename,'',[1 0 1 9] );
    % P{2} = dlmread(filename,'',4,0);
    %[row col] = find(cellfun(@(x) strcmp(x,'#'),n)) original find code for
    %cell boxes with only # in
    %sprintf('a26293_%i.txt',i); appends the i to the name
    %files = dir (uigetdir('C:')) %gets all filenames in a selected folder and puts them into files struct


         %                 figname = [files{i}(1:end-4),'.tif'];


        %tdata = load('topocmap.mat');  %load into structure

