% 
clear, clf, clc

%%% ! It is based on some internal Matlab functions. Could not work with some Matlab versions


%%%-----------------------------------------------------------------------%%%
%%% OPTIONS

%%% Expected peek properties (discrimination criteria):
%%% Wmax/Wmin - maximum and minimum linewidth (peak-to-peak) in Gauss

Wmax = 50; 
Wmin = 2; 

%%% HeightMin - minimum height of a peak in the units of intensity (a.u.)

HeightMin = 10; 


%%% Pre-processing of the spectra:
%%% smtwindow - number of points around each datapoint taking for averaging (1 - no smoothing)

smtwindow = 1;


%%%-----------------------------------------------------------------------%%%



    

% Searching for .spc files
listing = dir('*.spc');
numfiles = size(listing,1);

field = 1;
disp('The filename format is expected to be "field1_field2_field3_... .spc" ')
disp('Example of a current input:  ')
disp(listing(1).name)
disp('')
disp('Which field contains the information about the rotation angle?'  )
field = input('type: 1 or 2 or 3 ...:    '  )


% Trying to upload the parameters from file, if they were previously saved
parupload = 0;
pardata = [];
if exist('parameters.fitpar', 'file')
    pardata = load('parameters.fitpar');
    if size(pardata,1) == numfiles
        parupload = 1;
    end
end
    

for ispc = 1:numfiles
    % Trying to upload the parameters from file, if they were previously saved
    if parupload == 1 
        Wmax = pardata(ispc,1);
        Wmin = pardata(ispc,2);
        HeightMin = pardata(ispc,3); 
        smtwindow = pardata(ispc,4);
    end
    satisfied = 'n';
    while not(strcmp(satisfied,'y'))
           
        % Uploading a spectrum
        ifilename = listing(ispc).name;
        [x,y,Params] = eprload(ifilename,'G');
    
        % Picking a rotation angle from the filename 
        tmp = strsplit(ifilename(1:end-4), '_');
        output(ispc,1) = str2double(tmp(field));
    
   
        % Smoothing (datasmooth is a Matlab function)
        y = datasmooth(y, smtwindow,'savgol');



        % Finding peaks:
        % - finding local minima and maxima of a first derivative (findpeaks is a Matlab function)
        [pksmax,imax] = findpeaks( y,'MinPeakProminence', HeightMin );
        [pksmin,imin] = findpeaks(-y,'MinPeakProminence', HeightMin );
        locmax = x(imax);
        locmin = x(imin);
        % - searching for the lines which has widths within the interval [Wmin,Wmax]
        num = 1;
        loc = [];
        tmp1 = [];
        tmp2 = [];
        if not(isempty(locmax)) & not(isempty(locmin))
            for i = 1:max(size(locmax))
                for j = 1:max(size(locmin))
                    if (locmin(j) - locmax(i)) < Wmax & (locmin(j) - locmax(i)) > Wmin
                       tmp1 = findpeaks( y(imax(i):imin(j)),'MinPeakProminence', HeightMin );
                       tmp2 = findpeaks(-y(imax(i):imin(j)),'MinPeakProminence', HeightMin );
                       if isempty(tmp1) & isempty(tmp2)
                           loc(num) = (locmin(j)+locmax(i))/2;
                           locmin(j) = Inf;
                           locmax(i) = Inf;
                           num = num + 1;
                       end
                    end
                end
            end
        end
    
        
        % Showing the result and changing the parameters
        plot(x,y)
        for i = 1:max(size(loc))
            line([loc(i) loc(i)], [max(y) min(y)], 'LineStyle', '--', 'Color', 'r')
        end
    

        prompt = {'Minimum linewidth (G):','Maximum linewidth (G):', 'Minimum height (a.u.):', 'Smoothing window (points):', 'Do you like the results? (y/n)'};
        dlg_title = 'Change criteria';
        num_lines = 1;
        defaultans = {num2str(Wmin),num2str(Wmax),num2str(HeightMin),num2str(smtwindow),'y'};
        answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
    
        satisfied = answer(5);
        if strcmp(satisfied,'n')
            Wmin = str2double(answer(1)); 
            Wmax = str2double(answer(2)); 
            HeightMin = str2double(answer(3)); 
            smtwindow = str2double(answer(4));
        end
        
    end
    
    parameters(ispc,:) = [Wmax , Wmin , HeightMin , smtwindow];
    
    % Adding peak positions to the output
    output(ispc,2:size(loc,2)+1) = loc;
end


% Saving output and parameters
save ang_dep_peak.dat output -ascii
save parameters.fitpar parameters -ascii








