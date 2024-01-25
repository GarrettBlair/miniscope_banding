% This is a very rough script that provides an interactive method for
% correcting the static banding issue that can appear in miniscope
% recordings. This isn't meant to fix the 'moving' type which are
% relatively short lived. I recommend averaging those out for now.
%
% To run  -  Load in the data you want to fix (can be just one good and one bad frame)
% and identify the good and bad frames. All of the bad frames should be of
% the same type (i.e. the same static banding problem). The end of the
% script will apply the fixes to all of the bad frames, but currently
% doesn't save the output image file
%
% planned improvements -    
%       1. means of correcting different 'species'
%       2. a more automated quantification of correction quality using 
%       correlation with the good frame
%       3. python implementation and an actual GUI (please help)
% 
% 
% version = 1.0   -  Jan 25, 2024
% Garrett J. Blair, PhD
% NYU Center for Neural Science
% garrettjblair.com  |  gblair92@gmail.com

%%%%%%%%%%%%%%%%%%%%% Change these to your data %%%%%%%%%%%%%%%%%%%%%
tiffname = "C:/Users/Garrett/Desktop/test_stripes/test2.tiff";
good_frame_ind = 850; % first use a known good frame
bad_frames_ind = good_frame+1:nf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


bad_frame           = Y(:,:,bad_frames_ind(1));
good_frame          = Y(:,:,good_frame);

use_default_buffer  = true; % these values should work for most uncropped and non downsamples miniscope videos

% some experimental stuff
find_bad_frames     = false; % experimental way of finding bad frames
quick_detect        = false; % if true uses fewer pixels biased towards the center of frame when fidning bad frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load in the video
% this might be memory intensive if you have small RAM
Y = sub_imread_big(tiffname); % good fast tiff reader, see bottom fnc
[h, w, nf] = size(Y);
% you could instead load each frame when needed

%% A decent way of identifying banding frames, WIP
if find_bad_frames == true % enter the frame indices that are all of the same banding characteristic
    c =  zeros(nf,1);
    if quick_detect==true
        npx      = 2000; % increase for better detection, decrease for speed
        seeds    = randn(npx,2); % normally dist bias to center of frame
        %         seeds = rand(npx,2); % can use uniform alternatively
        seeds    = (seeds-min(seeds(:))) / max(max(seeds-min(seeds(:)))); % normalize
        randr    = floor(seeds(:,1)*(h-1))+1; % convert to subscripts
        randc    = floor(seeds(:,2)*(w-1))+1; % convert to subscripts
        randinds = sub2ind([h,w], randr, randc); % convert to indices
        clearvars seeds randr randc
    else
        ds = 4;
        randinds = 1:ds:h*w;
    end
    im = good_frame(:);
    testpx1 = double(im(randinds));
    for i = 1:nf
        %%
        bad_frame=Y(:,:,i);
        bad_frame=bad_frame(:);
        
        testpx2 = double(bad_frame(randinds));
        c(i) = corr(testpx1, testpx2);
    end
    figure;
    hold on
    plot(c)
    scatter(good_frame_ind, c(good_frame), 'go')
    title('Correlation with good frame')
    bad_frames_ind = find(c <  .99); % change the thresh here if needed
end
%%


if use_default_buffer==true
    minbuffer = 7000*2; % for full FOV v4 miniscope    2*xid*h;
    maxbuffer = 9000*2; % for full FOV v4 miniscope    4*(xid+1)*h;
    currentoffset = 8076*2; %     mean([minoffset,maxoffset]);
    currentbuffer = currentoffset; %  mean([minbuffer,maxbuffer]);
    
else % determine the buffer size(approx # of rows) by the row correlation with the good image
    cr =  zeros(h,1);
    for ind = 1:h
        r1 = single(good_frame(ind,:));
        r2 = single(bad_frame(ind,:));
        cr(ind) =  corr(r1', r2');
    end
    %%
    % find the correlation discontinuities
    x = abs( diff(cr) );
    % normalize the jumps
    x = (x + min(x)) / max(x+min(x));
    % threshold the small ones out
    x(x<.1) = 0;
    xi =  abs(diff(find(x>0)));
    xi = xi(xi>2); % remove the small variation if present, should be 10-16 rows typically
    xid =  mode( xi );
    xi_int =  xi/xid;
    
    minbuffer = 2*xid*h;
    maxbuffer = 4*(xid+1)*h;
    currentoffset = mean([minoffset,maxoffset]);
    currentbuffer = mean([minbuffer,maxbuffer]);
    
    figure;
    hold on
    plot(x>0)
    plot(cr)
    legend({'row corr.', 'Detected jump'})
    title(sprintf('Buffer size estimate: %d rows, about %d pixels', xid, currentbuffer))
    
end

minoffset = 0;
maxoffset = 2*maxbuffer;
swap_dist = 0; % distnace to swap the image data across bands
minswap = -5; % usually it'll be between -2 and 2
maxswap = 5;


linear_badframe = gb_mat2vec(bad_frame);

% can switch between these scales for more precise tuning, smallest should
% be a pixel
scales = [ .01 .001 .0001, .00001 ];
scaling_ind = 1;
current_scale = scales(1+mod(scaling_ind, length(scales)));

im_good = cat(3, good_frame, good_frame, good_frame);

resp = 0;
while resp ~= 113 % q to quit
    %%
    bufferguess = mod(linspace(0+currentoffset, h*w+currentoffset, h*w), currentbuffer)<=(currentbuffer/2);
    
    buff_id = bufferguess(2:end) + cumsum(diff(bufferguess)>0);
    buff_id(bufferguess==0)=0;
    figure(1);
    clf
    
    is_buffer = linear_badframe;
    fixed_im = linear_badframe;
    % there's prob a better way to identify and swap the interlacing tbh
    u_id = unique(buff_id(~isnan(buff_id)));
    u_id = u_id(u_id>0);
    for ii = 1:length(u_id)
        new_id = circshift(u_id, swap_dist);
        indsR = find(buff_id==new_id(ii));
        indsL = find(buff_id==u_id(ii));
        indsdiff = length(indsL) - length(indsR);
        if indsdiff<0 && abs(indsdiff)==1
            % unsure if this is corect, but only off by a pixel if wrong
            indsR = indsR(1:end+indsdiff);
        elseif indsdiff>0
            % unsure if this is corect, but only off by a pixel if wrong
            indsL = indsL(1:end-indsdiff);
        end
        try
            fixed_im(indsL) = linear_badframe(indsR);
        catch
%             disp(ii)
%             sum(buff_id==u_id(ii))
%             sum(buff_id==new_id(ii))
%             sum(buff_id==new_id(ii)) - sum(buff_id==u_id(ii))
        end
    end
    is_buffer(bufferguess==0) = 0;
    bufferim = reshape(is_buffer, [h, w]);
    bufferim = uint8(gb_vec2mat(bufferim, [h,w]));
    
    fixed_im = uint8(gb_vec2mat(fixed_im, [h,w]));
    fixed_im = cat(3, fixed_im, fixed_im, fixed_im);
    buffer_overlay = cat(3, bad_frame, bad_frame+bufferim/2, bad_frame);
    im = cat(2, im_good, buffer_overlay, fixed_im);
    
    image(im)
    axis image
    str = sprintf('Size (left/right): %d        Offset(up/dwn): %d        Scaling(x): %1.1e        Stripe offset(</>): %i',...
        currentbuffer, currentoffset, current_scale, swap_dist);
    title(str)
    str = sprintf('Controls:      Buffer Size ( up/dwn )      Buffer Offset ( left/right )      Scaling factor( x )      Stripe offset( </> )      (R)eset (Q)uit');
    xlabel(str)
    set(gca, 'YTick', [], 'XTick', [])
    drawnow
    
    [px,py,resp]=ginput(1);
    switch resp
        case 28 % left
            currentoffset = currentoffset - ceil(current_scale*maxoffset);
            currentoffset = max(minoffset, currentoffset);
        case 29 % right
            currentoffset = currentoffset + ceil(current_scale*maxoffset);
            currentoffset = min(maxoffset, currentoffset);
        case 30 % up
            currentbuffer = currentbuffer + ceil(current_scale*maxbuffer);
            currentbuffer = min(maxbuffer, currentbuffer);
        case 31 % down
            currentbuffer = currentbuffer - ceil(current_scale*maxbuffer);
            currentbuffer = max(minbuffer, currentbuffer);
        case 122 % z - reset vals
            currentoffset = 8076*2; %     mean([minoffset,maxoffset]);
            currentbuffer = currentoffset; %  mean([minbuffer,maxbuffer]);
            current_scale = scales(1);
            swap_dist     = 0;
        case 120 % x - change scale
            scaling_ind = scaling_ind+1;
            current_scale = scales(1+mod(scaling_ind, length(scales)));
        case 44 % < dec swap size
            swap_dist = swap_dist-1;
            swap_dist = max(minswap, swap_dist);
        case 46 % > increase swap size
            swap_dist = swap_dist+1;
            swap_dist = min(maxswap, swap_dist);
    end
end
disp(currentbuffer)
disp(currentoffset)



%% Apply these fixes to all of the bad frames
figure(1);
clf

bufferguess = mod(linspace(0+currentoffset, h*w+currentoffset, h*w), currentbuffer)<=(currentbuffer/2);
buff_id = bufferguess(2:end) + cumsum(diff(bufferguess)>0);
buff_id(bufferguess==0)=0;
Yfix = Y;
for i = 1:1:length(bad_frames_ind)
    im3=Y(:,:,bad_frames_ind(i));
    linear_badframe = gb_mat2vec(im3);

    fixed_im=linear_badframe;
    % does the same as above
    for ii = 1:length(u_id)
        new_id = circshift(u_id, swap_dist);
        indsR = find(buff_id==new_id(ii));
        indsL = find(buff_id==u_id(ii));
        indsdiff = length(indsL) - length(indsR);
        if indsdiff<0 && abs(indsdiff)==1
            indsR = indsR(1:end+indsdiff);
        elseif indsdiff>0
            indsL = indsL(1:end-indsdiff);
        end
        try
            fixed_im(indsL) = linear_badframe(indsR);
        catch
            %             disp(ii)
        end
    end
    fixed_im = uint8(gb_vec2mat(fixed_im, [h,w]));
    Yfix(:,:,bad_frames_ind(i)) = fixed_im;
    imshow([im3 fixed_im])
    drawnow
end

%%

function vec1 = gb_mat2vec(mat)
% matlab reshape wasn't working for me
[r,c] = size(mat);
vec1 = zeros(r*c,1);
for i=1:r
    for j=1:c
        vec1((i-1)*r + j) = mat(i,j);
    end
end
end

function mat1 = gb_vec2mat(vec, sz)
% matlab reshape wasn't working for me
r=sz(1);
c=sz(2);
mat1 = zeros(r,c);
for i=1:r
    for j=1:c
        mat1(i,j) = vec((i-1)*r + j);
    end
end
end

function [stack_out,Nframes] = sub_imread_big(stack_name,varargin)
% copied from https://www.mathworks.com/matlabcentral/fileexchange/61376-imread_big-read-in-tiff-stacks-larger-than-4gb
% by GJB on 2022_02_24
%Tristan Ursell
%Read large image stack (TIFF)
%May 2019
%
% This function can load image stacks larger than 4GB which is the typical
% limitation on image files.  Designed to work with single-channel 
% uncompressed TIFF stacks.  Also works for files smaller than 4GB.
%
% [stack_out,Nframes]= imread_big(stack_name);
% [stack_out,Nframes]= imread_big(stack_name,[i j]);
%
% stack_name = the path and file name to the image stack
%
% [i j] = optional frame number range to load (i = j is allowed)
%
% Nframes = number of frames as determined by total file size divided by
% estimated size of each frame data block
%
% stack_out = the output image stack of size [M N Nframes], where M x N is
% the size of each image.
%

%get data block size
info1 = imfinfo(stack_name);
stripOffset = info1(1).StripOffsets;
stripByteCounts = info1(1).StripByteCounts;
%get image size
sz_x=info1(1).Width;
sz_y=info1(1).Height;
if length(info1)<2
    Nframes=floor(info1(1).FileSize/stripByteCounts);
else
    Nframes=length(info1);
end
%read in arguments
if nargin>1
    temp1=varargin{1};
    ilow=temp1(1);
    ihigh=temp1(2);
else
    ilow=1;
    ihigh=Nframes;
end
%load file id tag
fID = fopen (stack_name, 'r');
if info1(1).BitDepth==32
    stack_out = zeros([sz_y sz_x ihigh-ilow+1],'uint32');
elseif info1(1).BitDepth==16
    stack_out = zeros([sz_y sz_x ihigh-ilow+1],'uint16');
else
    stack_out = zeros([sz_y sz_x ihigh-ilow+1],'uint8');
end
start_point = stripOffset(1) + (0:1:(Nframes-1)).*stripByteCounts + 1;
%start_point = stripOffset(1) + (0:1:(Nframes-1))*stripByteCounts + 1;
q=0;
for i = ilow:ihigh
    %fprintf ('loading image ... %d\n', i);
    fseek (fID, start_point(i), 'bof');
    
    if info1(1).BitDepth==32
        A = fread(fID, [sz_x sz_y], 'uint32=>uint32');
    elseif info1(1).BitDepth==16
        A = fread(fID, [sz_x sz_y], 'uint16=>uint16');
    else
        A = fread(fID, [sz_x sz_y], 'uint8=>uint8');
    end
    
    q=q+1;
    try
        stack_out(:,:,q) = A';
    catch
        disp(['Terminated at frame ' num2str(i) '.'])
        break
    end
end
%stack_out=stack_out(:,:,1:q);
Nframes=q;
fclose(fID);
end