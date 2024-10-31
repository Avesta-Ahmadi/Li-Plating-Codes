function [data, sdata] = preprocess()
%%% -----------------------------------------------------------------
% Process raw data to be used for analysis
% Input: -
% creating a struct (data) for storing raw data
% Output: data
% R1 = C20 rate
% R2 = C3 rate
% R3 = C2 rate
% R4 = 1C rate
% R5 = 2C rate
% R6 = 3C rate
% save data in a user friendly .mat file for future processing

data = struct("R1", struct(), "R2", struct(), "R3", struct(), "R4", struct(), "R5", struct(), "R6", struct());
sheets = ["C20", "C3", "C2", "1C", "2C", "3C"];
fieldnames = ["R1", "R2", "R3", "R4", "R5", "R6"];
colnames = ["steps","time","V","I","capacity","charge"];
for i=1:6
    data.(fieldnames(i)).NMR = readtable('raw_data.xlsx','Sheet',sheets(i)+'-NMR', 'PreserveVariableNames', true);
    data.(fieldnames(i)).echem = readtable('raw_data.xlsx','Sheet',sheets(i)+'-echem', 'PreserveVariableNames', true);
    data.(fieldnames(i)).echem.Properties.VariableNames = colnames;
    data.(fieldnames(i)).echem.time = data.(fieldnames(i)).echem.time/3600;
    if i==1
        data.(fieldnames(i)).NMR = data.(fieldnames(i)).NMR(:,[2,3,6,10,14,17,18])
        data.(fieldnames(i)).NMR.Properties.VariableNames = ["charge","time","electrolyte","phase1","phase2","total_solid","total"];
    end
    if i==2
        data.(fieldnames(i)).NMR = data.(fieldnames(i)).NMR(:,[2,3,6,10,14,18,22,25,26])
        data.(fieldnames(i)).NMR.Properties.VariableNames = ["charge","time","electrolyte","phase1","phase2","phase3","plated","total_solid","total"];
    end
    if i==3
        data.(fieldnames(i)).NMR = data.(fieldnames(i)).NMR(:,[2,3,6,10,14,18,22,25,26])
        data.(fieldnames(i)).NMR.Properties.VariableNames = ["charge","time","electrolyte","phase1","phase2","phase3","plated","total_solid","total"];
    end
    if i==4
        data.(fieldnames(i)).NMR = data.(fieldnames(i)).NMR(:,[2,3,6,10,14,18,22,26,29,30])
        data.(fieldnames(i)).NMR.Properties.VariableNames = ["charge","time","electrolyte","phase1","phase2","phase3","plated","dendrite","total_solid","total"];
    end
    if i==5
        data.(fieldnames(i)).NMR = data.(fieldnames(i)).NMR(:,[2,3,6,10,14,18,22,26,29,30])
        data.(fieldnames(i)).NMR.Properties.VariableNames = ["charge","time","electrolyte","phase1","phase2","phase3","plated","dendrite","total_solid","total"];
        data.R5.NMR.phase3(isnan(data.R5.NMR.phase3)) = 0;
        data.R5.NMR.dendrite(isnan(data.R5.NMR.dendrite)) = 0;
    end
    if i==6
        data.(fieldnames(i)).NMR = data.(fieldnames(i)).NMR(:,[2,3,6,10,14,18,22,26,29,30])
        data.(fieldnames(i)).NMR.Properties.VariableNames = ["charge","time","electrolyte","phase1","phase2","phase3","plated","dendrite","total_solid","total"];
    end
end


ss = string(getfield(data, 'R1').echem.steps);
ne = find(ss~="");
data.R1.echem.steps(ne(1)) = {"OCV pre charge"};
data.R1.echem.steps(ne(2)) = {"charge"};
data.R1.echem.steps(ne(3)) = {"OCV post charge"};
data.R1.echem.steps(ne(4)) = {"discharge"};
data.R1.echem.steps(ne(5)) = {"OCV post discharge"};

ss = string(getfield(data, 'R2').echem.steps);
ne = find(ss~="");
data.R2.echem.steps(ne(1)) = {"OCV pre charge"};
data.R2.echem.steps(ne(2)) = {"charge"};
data.R2.echem.steps(ne(3)) = {"OCV post charge"};
data.R2.echem.steps(ne(4)) = {"discharge"};
data.R2.echem.steps(ne(5)) = {""};
data.R2.echem.steps(ne(6)) = {"OCV post discharge"};

ss = string(getfield(data, 'R3').echem.steps);
ne = find(ss~="");
data.R3.echem.steps(ne(1)) = {"OCV pre charge"};
data.R3.echem.steps(ne(2)) = {"charge"};
data.R3.echem.steps(ne(3)) = {""};
data.R3.echem.steps(ne(4)) = {"OCV post charge"};
data.R3.echem.steps(ne(5)) = {"discharge"};
data.R3.echem.steps(ne(6)) = {""};
data.R3.echem.steps(ne(7)) = {"OCV post discharge"};

ss = string(getfield(data, 'R4').echem.steps);
ne = find(ss~="");
data.R4.echem.steps(ne(1)) = {"OCV pre charge"};
data.R4.echem.steps(ne(2)) = {"charge"};
data.R4.echem.steps(ne(3)) = {""};
data.R4.echem.steps(ne(4)) = {"OCV post charge"};
data.R4.echem.steps(ne(5)) = {"discharge"};
data.R4.echem.steps(ne(6)) = {""};
data.R4.echem.steps(ne(7)) = {"OCV post discharge"};

ss = string(getfield(data, 'R5').echem.steps);
ne = find(ss~="");
data.R5.echem.steps(ne(1)) = {"OCV pre charge"};
data.R5.echem.steps(ne(2)) = {"charge"};
data.R5.echem.steps(ne(3)) = {""};
data.R5.echem.steps(ne(4)) = {"OCV post charge"};
data.R5.echem.steps(ne(5)) = {"discharge"};
data.R5.echem.steps(ne(6)) = {""};
data.R5.echem.steps(ne(7)) = {"OCV post discharge"};

ss = string(getfield(data, 'R6').echem.steps);
ne = find(ss~="");
data.R6.echem.steps(ne(1)) = {"OCV pre charge"};
data.R6.echem.steps(ne(2)) = {"charge"};
data.R6.echem.steps(ne(3)) = {"OCV post charge"};
data.R6.echem.steps(ne(4)) = {"discharge"};
data.R6.echem.steps(ne(5)) = {""};
data.R6.echem.steps(ne(6)) = {"OCV post discharge"};



save('data.mat', '-struct', 'data');









% prepare data and store in sdata: compute C_1, C_2, time, and label for
% state of cell: charge, OCV, discharge
steps = ["OCV pre charge"; "charge"; "OCV post charge"; "discharge"; "OCV post discharge"];
sdata = {};

% C20 cycle
tt = data.R1.echem.time(string(data.R1.echem.steps) == "discharge");       
ind_discharge = size(data.R1.NMR.time,1) - sum(data.R1.NMR.time>tt);
tt = data.R1.echem.time(string(data.R1.echem.steps) == "OCV post charge");       
ind_rest = size(data.R1.NMR.time,1) - sum(data.R1.NMR.time>tt);
c1 = data.R1.NMR.phase1 + data.R1.NMR.phase2;
t = data.R1.NMR.time*3600;
[~,indices] = unique(data.R1.echem.time);
J = interp1(data.R1.echem.time(indices), data.R1.echem.I(indices),data.R1.NMR.time,'spline');
J(J==0) = 1e-15; 
label = strings(size(t));
label(1:ind_rest) = 'charge';
label(ind_rest:ind_discharge) = 'OCV';
label(ind_discharge:end) = 'discharge';
sdata.R1 = table(t,J,c1,label);

% C3 cycle
tt = data.R2.echem.time(string(data.R2.echem.steps) == "discharge");       
ind_discharge = size(data.R2.NMR.time,1) - sum(data.R2.NMR.time>tt);
tt = data.R2.echem.time(string(data.R2.echem.steps) == "OCV post charge");       
ind_rest = size(data.R2.NMR.time,1) - sum(data.R2.NMR.time>tt);
data.R2.NMR.phase3(isnan(data.R2.NMR.phase3)) = 0;                         % replace NAN with zero
c1 = data.R2.NMR.phase1 + data.R2.NMR.phase2 + data.R2.NMR.phase3;
c2 = data.R2.NMR.plated;
t = data.R2.NMR.time*3600;
[~,indices] = unique(data.R2.echem.time);
J = interp1(data.R2.echem.time(indices), data.R2.echem.I(indices),data.R2.NMR.time,'spline');
J(J==0) = 1e-15; % to prevent division by zero.
label = strings(size(t));
label(1:ind_rest) = 'charge';
label(ind_rest:ind_discharge) = 'OCV';
label(ind_discharge:end) = 'discharge';
sdata.R2 = table(t,J,c1,c2, label);

% C2 cycle
tt = data.R3.echem.time(string(data.R3.echem.steps) == "discharge");       
ind_discharge = size(data.R3.NMR.time,1) - sum(data.R3.NMR.time>tt);
tt = data.R3.echem.time(string(data.R3.echem.steps) == "OCV post charge");       
ind_rest = size(data.R3.NMR.time,1) - sum(data.R3.NMR.time>tt);
data.R3.NMR.phase3(isnan(data.R3.NMR.phase3)) = 0;                         % replace NAN with zero
c1 = data.R3.NMR.phase1 + data.R3.NMR.phase2 + data.R3.NMR.phase3;
c2 = data.R3.NMR.plated;
t = data.R3.NMR.time*3600;
[~,indices] = unique(data.R3.echem.time);
J = interp1(data.R3.echem.time(indices), data.R3.echem.I(indices),data.R3.NMR.time,'spline');
J(J==0) = 1e-15; % to prevent division by zero.
label = strings(size(t));
label(1:ind_rest) = 'charge';
label(ind_rest:ind_discharge) = 'OCV';
label(ind_discharge:end) = 'discharge';
sdata.R3 = table(t,J,c1,c2,label);

% 1C cycle
tt = data.R4.echem.time(string(data.R4.echem.steps) == "discharge");       
ind_discharge = size(data.R4.NMR.time,1) - sum(data.R4.NMR.time>tt);
tt = data.R4.echem.time(string(data.R4.echem.steps) == "OCV post charge");       
ind_rest = size(data.R4.NMR.time,1) - sum(data.R4.NMR.time>tt);
data.R4.NMR.phase3(isnan(data.R4.NMR.phase3)) = 0;                         % replace NAN with zero
c1 = data.R4.NMR.phase1 + data.R4.NMR.phase2 + data.R4.NMR.phase3;
c2 = data.R4.NMR.plated + data.R4.NMR.dendrite;
t = data.R4.NMR.time*3600;
[~,indices] = unique(data.R4.echem.time);
J = interp1(data.R4.echem.time(indices), data.R4.echem.I(indices),data.R4.NMR.time,'spline');
J(J==0) = 1e-15; % to prevent division by zero.
label = strings(size(t));
label(1:ind_rest) = 'charge';
label(ind_rest:ind_discharge) = 'OCV';
label(ind_discharge:end) = 'discharge';
sdata.R4 = table(t,J,c1,c2,label);

% 2C cycle
tt = data.R5.echem.time(string(data.R5.echem.steps) == "discharge");       
ind_discharge = size(data.R5.NMR.time,1) - sum(data.R5.NMR.time>tt);
tt = data.R5.echem.time(string(data.R5.echem.steps) == "OCV post charge");       
ind_rest = size(data.R5.NMR.time,1) - sum(data.R5.NMR.time>tt);
data.R5.NMR.phase3(isnan(data.R5.NMR.phase3)) = 0;                         % replace NAN with zero
c1 = data.R5.NMR.phase1 + data.R5.NMR.phase2 + data.R5.NMR.phase3;
c2 = data.R5.NMR.plated + data.R5.NMR.dendrite;
t = data.R5.NMR.time*3600;
[~,indices] = unique(data.R5.echem.time);
J = interp1(data.R5.echem.time(indices), data.R5.echem.I(indices),data.R5.NMR.time,'spline');
J(J==0) = 1e-15; % to prevent division by zero.
label = strings(size(t));
label(1:ind_rest) = 'charge';
label(ind_rest:ind_discharge) = 'OCV';
label(ind_discharge:end) = 'discharge';
sdata.R5 = table(t,J,c1,c2,label);

% 3C cycle
tt = data.R6.echem.time(string(data.R6.echem.steps) == "discharge");       
ind_discharge = size(data.R6.NMR.time,1) - sum(data.R6.NMR.time>tt);
tt = data.R6.echem.time(string(data.R6.echem.steps) == "OCV post charge");       
ind_rest = size(data.R6.NMR.time,1) - sum(data.R6.NMR.time>tt);
data.R6.NMR.phase3(isnan(data.R6.NMR.phase3)) = 0;                         % replace NAN with zero
c1 = data.R6.NMR.phase1 + data.R6.NMR.phase2 + data.R6.NMR.phase3;
c2 = data.R6.NMR.plated + data.R6.NMR.dendrite;
t = data.R6.NMR.time*3600;
[~,indices] = unique(data.R6.echem.time);
J = interp1(data.R6.echem.time(indices), data.R6.echem.I(indices),data.R6.NMR.time,'spline');
J(J==0) = 1e-15; % to prevent division by zero.
label = strings(size(t));
label(1:ind_rest) = 'charge';
label(ind_rest:ind_discharge) = 'OCV';
label(ind_discharge:end) = 'discharge';
sdata.R6 = table(t,J,c1,c2,label);


% Normalize variables, rescale ODE system, define tm, cm, sf
% tm and cm are kinda chosen arbitrarily. cm is the max concentration seen
% in C3 cycle. tm is the time taken for charge and rest.
cm = max(sdata.R2.c1 + sdata.R2.c2);
tm = min(sdata.R2.t(sdata.R2.label=='discharge'));
sf = tm/cm;
cm1 = max(sdata.R3.c1);
cm2 = max(sdata.R6.c2);
im = max(sdata.R4.J);

% rescale variables
sdata.R2.t  = sdata.R2.t/tm;
sdata.R2.c1 = sdata.R2.c1/cm;
sdata.R2.c2 = sdata.R2.c2/cm;
sdata.R3.t  = sdata.R3.t/tm;
sdata.R3.c1 = sdata.R3.c1/cm;
sdata.R3.c2 = sdata.R3.c2/cm;
sdata.R4.t  = sdata.R4.t/tm;
sdata.R4.c1 = sdata.R4.c1/cm;
sdata.R4.c2 = sdata.R4.c2/cm;
sdata.R5.t  = sdata.R5.t/tm;
sdata.R5.c1 = sdata.R5.c1/cm;
sdata.R5.c2 = sdata.R5.c2/cm;
sdata.R6.t  = sdata.R6.t/tm;
sdata.R6.c1 = sdata.R6.c1/cm;
sdata.R6.c2 = sdata.R6.c2/cm;

sdata.R2.J  = sdata.R2.J/im;
sdata.R3.J  = sdata.R3.J/im;
sdata.R4.J  = sdata.R4.J/im;
sdata.R5.J  = sdata.R5.J/im;
sdata.R6.J  = sdata.R6.J/im;



% save this version of data for further use
save('sdata.mat', "sdata", "cm", "tm", "sf","cm1","cm2");


end
