function BZreaction4
% Belousov-Zhabotinsky Reaction animation
% This MATLAB code is converted from Processing code available in this link
% http://www.aac.bartlett.ucl.ac.uk/processing/samples/bzr.pdf

% version 2. Corrected the drift of pixels as suggested by Jonh.

xres = 100; %x resolution
yres = 100; %y resolution
p = 1;
q = 2;

a = rand(xres,yres,2);
b = rand(xres,yres,2);
c = rand(xres,yres,2);
img = zeros(xres,yres,3);
c_a = zeros(xres,yres);
c_b = zeros(xres,yres);
c_c = zeros(xres,yres);

mm = mod((1:xres+2)+xres,xres)+1;
nn = mod((1:yres+2)+yres,yres)+1;
[mm,nn] = meshgrid(mm,nn);
idx = sub2ind([yres xres],nn(:),mm(:));
idx = reshape(idx,[yres xres]+2);

figure
h1 = image(img);
axis equal off

k = 0;    
while k < 100 && ishandle(h1)
    c_a = 0*c_a;
    c_b = 0*c_b;
    c_c = 0*c_c;
    
    for m=1:xres
        for n=1:yres
            idx_temp = idx(m:m+2,:);
            idx_temp = idx_temp(:,n:n+2);
            idx_temp = idx_temp(:);
            if p == 2
                idx_temp = idx_temp+(xres+0)*(yres+0);
            end
            
            c_a(m,n) = c_a(m,n)+ sum(a(idx_temp));
            c_b(m,n) = c_b(m,n)+ sum(b(idx_temp));
            c_c(m,n) = c_c(m,n)+ sum(c(idx_temp));            
        end
    end
    
    %correction of pixel drift
    c_a = circshift(c_a,[2 2]);
    c_b = circshift(c_b,[2 2]);
    c_c = circshift(c_c,[2 2]);
    
    c_a = c_a / 9.0;
    c_b = c_b / 9.0;
    c_c = c_c / 9.0;
    
    a(:,:,q) = double(uint8(255*(c_a + c_a .* (c_b - c_c))))/255;
    b(:,:,q) = double(uint8(255*(c_b + c_b .* (c_c - c_a))))/255;
    c(:,:,q) = double(uint8(255*(c_c + c_c .* (c_a - c_b))))/255;
    
    img(:,:,1)=c(:,:,q);
    img(:,:,2)=b(:,:,q);
    img(:,:,3)=a(:,:,q);
    
    if p == 1
        p = 2; q = 1;
    else
        p = 1; q = 2;
    end
    
    set(h1, 'CData', hsv2rgb(img));
    pause(0.1);
    k = k + 1;
end
end
