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
%surfdata = cell(1,length(files))
%for i = 1:length(files)
%    surfdata{i} = isurf(images{i})
%end

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
    Fmatrices{i} = matches{i}.ransac(@fmatrix, 1e-2, 'verbose');
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
EulerAngles = cell(1, length(Ematrices));
Theta = zeros(1,length(Ematrices));

W = [0 -1 0
     1  0 0
     0  0 1];

Z = [ 0 1 0
     -1 0 0
      0 0 0];

for i = 1:length(Ematrices)
    [U, S, V] = svd(Ematrices{i})
    R = U*W'*V' % note W' = inv(W)
    txmat = V*W*S*V'
    Rmatrices{i}=R;
    TXmatrices{i}=txmat;
    TX{i} = [txmat(3,2) txmat(1,3) txmat(2,1)];
    EulerAngles{i} = [atan2(R(3,2),R(3,3)), atan2(-R(3,1),sqrt(R(3,2)*R(3,2) + R(3,3)*R(3,3))), atan2(R(2,1),R(1,1))]
    Theta(i) = atan2(R(2,1),R(1,1));
end

xVals = zeros(1,length(Theta));
yVals = zeros(1,length(Theta));
for i = 1:length(Theta)
    xVals(i) = TX{i}(1)
    yVals(i) = TX{i}(2)
end
%plot(xVals, yVals)
for i = 1:length(EulerAngles)
    EulerAngles{i}
end
plotValues = cell(1,length(Theta))
for i=1:length(plotValues)
    plotValues{i} = se2(TX{i}(1),TX{i}(2),Theta(i))
end
%for i=1:length(plotValues)
%    axis([-400 150 -300 150])
%    newVal = plotValues{i}
%    trplot2(newVal, 'frame', int2str(i), 'color','b')
%    hold on
%end
%axis([0 5 0 5])
%trplot2(plotValues{1}, 'frame', 1, 'color','b')
%hold on
%trplot2(plotValues{2}, 'frame', 2, 'color','b')
%plot_point([1,2],'*')
%%
%im1 = iread('.\dataset\pic-1426087272.jpg', 'mono', 'double')
%im2 = iread('.\dataset\pic-1426087275.jpg', 'mono', 'double');
%s1 = isurf(im1)
%s2 = isurf(im2)
%m = s1.match(s2)
%idisp({im1, im2})
%m.subset(100).plot('w')
%F = m.ransac(@fmatrix, 1e-3, 'verbose')

% Simultaneous Localization and Mapping
%clear all;
%close all;

%path(path, 'aprilTag');

% Specify the locations of the April tags in the Map
%A =  [5,6; 7,8; 1,0];
%% and a robot with noisy odometry

V=diag([0.1, 1.1*pi/180].^2);
veh=GenericVehicle(V,'dt',0.1);
veh.add_driver(DeterministicPath('berry_garg/run4/log-1423547652.txt'));

% Creating the map. It places landmarks according to 'A' matrix.
%map = LandmarkMap(20, A);

% Creating the sensor.  We firstly define the covariance of the sensor measurements
% which report distance and bearing angle
W = diag([0.1, 1*pi/180].^2);

% and then use this to create an instance of the Sensor class.
%sensor = GenericRangeBearingSensor(veh, map, W, 'animate');
% Note that the sensor is mounted on the moving robot and observes the features
% in the world so it is connected to the already created Vehicle and Map objects.

% Create the filter.  First we need to determine the initial covariance of the
% vehicle, this is our uncertainty about its pose (x, y, theta)
P0 = diag([0.005, 0.005, 0.001].^2);

% Now we create an instance of the EKF filter class
ekf = GenericEKF(veh, V, P0);
% and connect it to the vehicle and the sensor and give estimates of the vehicle
% and sensor covariance (we never know this is practice).

% Now we will run the filter for 1000 time steps.  At each step the vehicle
% moves, reports its odometry and the sensor measurements and the filter updates
% its estimate of the vehicle's pose
ekf.run(1000);
% all the results of the simulation are stored within the EKF object

% First let's plot the map
%clf; %map.plot()
% and then overlay the path actually taken by the vehicle
%veh.plot_xy('b');
% and then overlay the path estimated by the filter
ekf.plot_xy('r');
% which we see are pretty close

% Now let's plot the error in estimating the pose
%ekf.plot_error()
% and this is overlaid with the estimated covariance of the error.

% Remember that the SLAM filter has not only estimated the robot's pose, it has
% simultaneously estimated the positions of the landmarks as well.  How well did it
% do at that task?  We will show the landmarks in the map again
% map.plot();
% and this time overlay the estimated landmark (with a +) and the 3sigma 
% uncertainty bounds as green ellipses
% ekf.plot_map(3,'g');
