% Author: Dominic Bett
% Bioinformatics Project: Diffusion Kernel(Matrix Exponential)

disp('Computing the diffusion kernel...')
[data, varnames, casenames] = tblread('/Users/dbett/Dropbox/PROJECTS/ERN/Experiment/Laplacian/lapl_avg_Network_950.txt', '\t');
diff_kernel = expm(data);
tblwrite(diff_kernel,varnames,casenames,'/Users/dbett/Dropbox/PROJECTS/ERN/Experiment/Diffusion/dk_avg_Network_950.txt','\t');
disp('Done!')