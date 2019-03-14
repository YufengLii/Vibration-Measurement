function dimg = mipcentraldiff(img,direction)

% MIPCENTRALDIFF ? ? Finite difference calculations?

%

% ? DIMG = MIPCENTRALDIFF(IMG,DIRECTION)

%

% ?Calculates the central-difference for a given direction

% ?IMG ? ? ? : input image

% ?DIRECTION : 'dx' or 'dy'

% ?DIMG ? ? ?: resultant image

%

img = padarray(img,[1 1],'symmetric','both');

[row,col] = size(img);

dimg = zeros(row,col);

switch (direction)
    case 'dx'
        dimg(:,2:col-1) =(img(:,3:col)-img(:,1:col-2))/2;
    case 'dy'
        dimg(2:row-1,:) =(img(3:row,:)-img(1:row-2,:))/2;
    otherwise
        disp('Direction is unknown');
end

dimg = dimg(2:end-1,2:end-1);
