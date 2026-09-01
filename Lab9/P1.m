%% Part 1

image = imread("thorax_t1.jpg");

image = image(:,:,1);

image = double(image);

figure;
imshow(image,[]);
title('Select The Organs');

[x_lungl,y_lungl] = ginput(1);
disp([x_lungl,y_lungl]);

[x_lungr,y_lungr] = ginput(1);
disp([x_lungr,y_lungr]);

[x_liver,y_liver] = ginput(1);
disp([x_liver,y_liver]);

lungl = round([y_lungl, x_lungl]);
lungl_intensity = image(lungl(1), lungl(2));

lungr = round([y_lungr, x_lungr]);
lungr_intensity = image(lungr(1), lungr(2));

liver = round([y_liver, x_liver]);
liver_intensity = image(liver(1), liver(2));

T1 = 20;
T2 = 25;

mask_lungl = SelectRegion(image, lungl, lungl_intensity, T1);
mask_lungr = SelectRegion(image, lungr, lungr_intensity, T1);
mask_lung = mask_lungl | mask_lungr;

mask_liver = SelectRegion(image, liver, liver_intensity, T2);


% Display
figure;
Overlay(image, mask_lung);
title('Lungs');

figure;
Overlay(image, mask_liver);
title('Liver');


%%

function mask = SelectRegion(image, centre_point, centre_point_intensity, threshold)

[m,n] = size(image);

mask = false(m,n);      
checked = false(m,n);  

queue = centre_point;         

while ~isempty(queue)

    pixel = queue(1,:);
    queue(1,:) = [];

    r = pixel(1);
    c = pixel(2);

    % check image boundries
    if r<1 || r>m || c<1 || c>n
        continue;
    end

    % skip the pixel if checked
    if checked(r,c)
        continue;
    end

    checked(r,c) = true;

    % compare intensity 
    if abs(image(r,c) - centre_point_intensity) <= threshold

        mask(r,c) = true;

        % add neighbour points
        queue = [queue;
                 r-1 c;
                 r+1 c;
                 r c-1;
                 r c+1];
    end

end

end


%%

function Overlay(f, mask)

    m = max(f(:));
    fr = f;
    fg = f + mask / max(mask(:)) * m/2;
    fb = f;
   
    imshow(reshape([fr fg fb],[size(f,1) size(f,2) 3])/m, []);
    
    drawnow;
end