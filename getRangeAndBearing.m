function [z, js] = getRangeAndBearing(filename)
[pathstr, name, ext] = fileparts(filename);
image_width = 1280;
focal_length_in_mm = 4;
ccd_width_in_mm = 4; % given 1/4" ccd width, http://www.dpreview.com/articles/8095816568/sensorsizes
focal_length = image_width * focal_length_in_mm/ccd_width_in_mm; % http://www.cs.cornell.edu/~snavely/bundler/focal.html
actual_distance = 0.1077631; %m distance of centroid to bottom left
switch ext
    case '.pnm'
        disp(name)
        stats = getAprilTag_mex(filename);
        
        % order = top-left top-right bottom-right bottom-left
        centroid =  [(stats(:,2) + stats(:,4) + stats(:,6) + stats(:,8))./4 (stats(:,3) + stats(:,5) + stats(:,7) + stats(:,9))./4];
        bottom_left = [stats(:,8) stats(:,9)];
        dist = centroid-bottom_left;
        image_distance = sqrt(sum(dist.^2,2));
        range = focal_length*actual_distance./image_distance;

        centroidX = centroid(:,1);
        left = stats(:,8);
        bearing = atan((centroidX - left)/focal_length);
        z = [range bearing ];   % range & bearing measurement       
        disp('passed bearing')
        js = stats(:,1);
        disp(js)
    otherwise
        disp(name)
        error('Accepts only pnm files');
end