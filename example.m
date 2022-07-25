addpath('getLebedevSphere')
addpath('splinepak/splinepak')
addpath('tensor_toolbox-v3.1')

%load data
load(fullfile('data', 'HCP_scan1.mat'))

%compute mean
Y_bar = mean(Y_scan_1, 3);

N = size(Y_scan_1,3);
nv0 = length(X0);
nv1 = length(X1);

%cented the data and score in cell array
Y = {};
for i=1:N
    Y{i} = sparse(Y_scan_1(:,:,i) - Y_bar);
end 

clear Y_scan_1

%load spherical triangulations 
tri_file_0 = 'data/tessellations/tess_0'; 
tri_file_1 = 'data/tessellations/tess_1';
[m_0,x_0,y_0,z_0,nt_0,TRI_0] = sreadtri(tri_file_0);
[m_1,x_1,y_1,z_1,nt_1,TRI_1] = sreadtri(tri_file_1);

%create marginal basis objects 
b0 = SplineBasis([x_0,y_0,z_0]);
b1 = SplineBasis([x_1,y_1,z_1]);

%instantiate ConCon object 
K = 22; %rank 
CCB = ConConBasis(b0, b1, 'K', K);

%fit ConCon object to data (see documentation for hyper-parameter options)
CCB.Fit(Y,X0,X1);

%reprepent the data using the ConConBasis object
CCR1 = ConConSmooth(CCB);
CCR1.smooth(Y, X0, X1);

%load scan-2 data and represent with the CC-basis learned from scan-1
load(fullfile('data', 'HCP_scan2.mat'))

Y2 = {};
Y_bar = mean(Y_scan_2, 3);
for i=1:N
    Y2{i} = sparse(Y_scan_2(:,:,i) - Y_bar);
end

clear Y_scan_2

CCR2 = ConConSmooth(CCB);
CCR2.smooth(Y2, X0, X1);

%compute LOOCV and plot pairwise distance matrix 
Smat_scan1 = get(CCR1, 'CoefMat');
Smat_scan2 = get(CCR2, 'CoefMat');

norms_scan1 = vecnorm(Smat_scan1');
norms_scan2 = vecnorm(Smat_scan2');
Smat_scan1_norm = bsxfun(@rdivide, Smat_scan1, norms_scan1');
Smat_scan2_norm = bsxfun(@rdivide, Smat_scan2, norms_scan2');

Block_Smat = [];
for i=1:N
    Block_Smat = [Block_Smat; [Smat_scan1_norm(i,:); Smat_scan2_norm(i,:)]];
end 

D = squareform(pdist(Block_Smat));
imagesc(D)

% LOOCV 
subject_hat = zeros(size(D,1),1)';
subject_true = repelem(1:N,2);
for i=1:size(D,1)
    [out,idx] = sort(D(i,:));
    subject_hat(i) = subject_true(idx(2)); 
end 
sum(subject_true == subject_hat)/length(subject_true)

    