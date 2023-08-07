function [Data] = ReadDNSData3(filename) % x is the output of the function[x] where the filename here is the input

tic % starting elapsed time to keep track of the run time
%creating the 3d grid
Nx = 330; % number of blocks in x direction
Ny = 165; % number of blocks in y direction
Nz = 165; % number of blocks in z direction

% format of using textscan start
fid = fopen(filename, 'rb');
formatSpec = '%*f%*f%*f%f%f%f%*f%f%*f%*f%*f'; % '%*f' = skips column and '%f' = reads column so here columns 1 2 3 7 9 10 11 are skipped while columns 4 5 6 8 are considered
%cac = textscan(fid, formatSpec, 'Delimiter', ' ', 'MultipleDelimsAsOne', true);
cac = textscan(fid, formatSpec);
% main function that reads the .dat data file
fclose(fid);
% format of using textscan end

% now the cac will have a size of 1 x 4 where the 4 columns will correspond to the 4 columns of interest
u = cac{1}; % 1st column of cac is the 4th of the data we are reading
v = cac{2}; % 2nd column of cac = 5th column of data
w = cac{3}; % 3rd column of cac = 6th column of data
vof = cac{4}; % 4th coluomn of cac = 8th column of data

% reshaping of the column vectors to 3d matrix
u = reshape(u,[Nx Ny Nz]); 
v = reshape(v,[Nx Ny Nz]); 
w = reshape(w,[Nx Ny Nz]); 
vof = reshape(vof,[Nx Ny Nz]); 

% changing orientation just cause
u = permute(u,[2 1 3]);
v = permute(v,[2 1 3]);
w = permute(w,[2 1 3]);
vof = permute(vof,[2 1 3]);

Data = struct('u',u,'v',v,'w',w,'vof',vof); % organizing the output in a struct array
toc % this is for the tic-toc of the function
end 
