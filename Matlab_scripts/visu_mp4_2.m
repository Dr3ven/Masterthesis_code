clear

name =('Slow_Eruption');

%% avi or mp4
v = VideoWriter([name '.avi']);
% v = VideoWriter([name '.avi']);

v.FrameRate = 11;  % Default 30
%v.Quality   = 75; % Default 75

open(v)
for count=1:1:200
    A  = imread([name,'_it', num2str(count), '.png']);
%    A = imresize(A,[1744 1548]);
%     A2 = A(:,2:end-2);
    writeVideo(v,A)
end
close(v)

%% gif
% delt = 0.1;
% for count=1:118
%     im         = imread([name '_' num2str(count) '.png']);
%     [imind,cm] = rgb2ind(im,256);
%     if (count == 1), imwrite(imind,cm,[name '.gif'],'gif','Loopcount',inf,'DelayTime',delt);
%     else             imwrite(imind,cm,[name '.gif'],'gif','WriteMode','append','DelayTime',delt);
%     end
% end
