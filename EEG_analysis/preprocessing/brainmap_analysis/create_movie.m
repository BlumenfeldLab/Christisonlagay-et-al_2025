video_name = '6_patients_common_m_StimR.mp4';
total_frames_num = 77;
image_suffix = 'stitched_';

outputVideo = VideoWriter(video_name, 'MPEG-4');
outputVideo.FrameRate = 3;
open(outputVideo)

imagenames = cell(total_frames_num, 1);
for i = 1 : total_frames_num
    imagenames{i} = strcat([image_suffix num2str(i) '.tiff']);
end

for i = 1 : length(imagenames)
    img = imread(imagenames{i});
    scale_factor = 1920 / size(img, 2);
    img = imresize(img, scale_factor);
    writeVideo(outputVideo, img)
end

close(outputVideo)