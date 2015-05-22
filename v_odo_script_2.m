% run each script one at a time
% Currell Berry, 4/1/2015
% this document uses cells -- type ctrl+enter to run a cell

%%this cell reads in the data
files = dir('berry_garg/run4/*.jpg')
images = cell(1,length(files))
for i = 1:length(files)
    images{i} = iread(strcat('run4/',files(i).name), 'mono', 'double')
end

%% this cell runs SURF on each image in images
surfdata = cell(1,length(files))
for i = 1:length(files)
    surfdata{i} = isurf(images{i})
end

%% this cell runs SIFT on each image in images
siftdata = cell(1,length(files))
for i = 1:length(files)
    siftdata{i} = isift(images{i})
end


%% now get all the correspondences
datatouse = siftdata
matches = cell(1,length(files)-1)
for i = 1:length(files)-1
    matches{i} = datatouse{i}.match(datatouse{i+1})
end

%% this cell runs RANSAC on each match in matches,
%  thereby generating a list of fundamental matrices
Fmatrices = cell(1,length(matches))
for i = 1:length(matches)
    Fmatrices{i} = matches{i}.ransac(@fmatrix, 1e-3, 'verbose');
end

%% this displays the ransac results for two adjacent images.
ili = 8 %ili = image left index
iri = 1+ili
idisp({images{ili},images{iri}})
matches{ili}.inlier.subset(100).plot('g')

%% now calculate essential matrices
 K = [1211.2959 0           657.15924
      0         1206.00512  403.17667
      0         0           1       ];
 Ematrices = cell(1,length(Fmatrices))
 for i = 1:length(Fmatrices)
     Ematrices{i} = K'*Fmatrices{i}*K;
 end
 
%% now calculate relative transformations
% formulas taken from wikipedia http://en.wikipedia.org/wiki/Essential_matrix
% note the TX matrices are the "cross product representation of t"
% not sure currently how to get t from them, will look into further.
Rmatrices = cell(1, length(Ematrices));
TXmatrices = cell(1, length(Ematrices));
TX = cell(1, length(Ematrices));
Theta = zeros(1,length(Ematrices))

W = [0 -1 0
     1  0 0
     0  0 1];

Z = [ 0 1 0
     -1 0 0
      0 0 0];

for i = 1:length(Ematrices)
    [U, S, V] = svd(Ematrices{i})
    R = U*W'*V' % note W' = inv(W)
    txmat = V*Z*V'
    Rmatrices{i}=R;
    TXmatrices{i}=txmat;
    TX{i} = [txmat(3,2) txmat(1,3) txmat(2,1)];
    Theta(i) = atan2(R(2,1),R(1,1))/10
end


%%
%im1 = iread('.\dataset\pic-1426087272.jpg', 'mono', 'double')
%im2 = iread('.\dataset\pic-1426087275.jpg', 'mono', 'double');
%s1 = isurf(im1)
%s2 = isurf(im2)
%m = s1.match(s2)
%idisp({im1, im2})
%m.subset(100).plot('w')
%F = m.ransac(@fmatrix, 1e-3, 'verbose')


%%
items = cell(1,length(Theta))
postuples = zeros(length(Theta),3)
for i = 1:length(Theta)
items{i} = [cos(Theta(i)) sin(Theta(i)) TX{i}(1)
-sin(Theta(i)) cos(Theta(i)) TX{i}(2)
0              0             1]
    postuples(i,1:3)=[TX{i}(1) TX{i}(2) Theta(i)]
end

%%

poses = cell(1,length(items)+1)
poses{1} = eye(3)
points = zeros(length(items)+1,2)
points(1,1:2)=[poses{1}(1,3) poses{1}(2,3)]
for i = 2:length(poses)
    poses{i} = poses{i-1}*items{i-1}
    points(i,1:2) = [poses{i}(1,3) poses{i}(2,3)]
end
plot(points(:,1),points(:,2))
