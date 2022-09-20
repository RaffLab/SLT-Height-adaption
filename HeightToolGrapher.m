%Matlab script to make 3d plots & dynamic plots of spatial parameters vs
%time & intensity
clear all;
inputData = readmatrix('/Users/thomas/Documents/PhD/Data/Project Laura/Ana2/Highspeed/Center vs outside test/Individual Z tracking/AndorS3rNb.csv');
% 
% %3d plot that isn't that useful
% figure('Name','Z vs time/intensity')
% plot3(inputData(:,1),inputData(:,2),inputData(:,3),'o');
% xlabel('Frame')
% ylabel('Z-axis')
% zlabel('Intensity')


%2d plots over time

%first generate the data orgnaised by frame
minFrame = min(inputData(:,1))+1;
maxFrame = max(inputData(:,1))+1;
maxIntensity = max(inputData(:,3));
maxZ = max(inputData(:,2));
frameCell = cell(maxFrame,1);
for i=1:length(inputData)
    frameCell{inputData(i,1)+1} = [frameCell{inputData(i,1)+1};inputData(i,:)];
end

v = VideoWriter('AS3rNbCentrioleHeight.avi');
open(v);
figure('Name','Intensity vs Z over time')
ylabel('Intensity')
xlabel('Z')
for i=1:length(frameCell)
    plot(frameCell{i}(:,2),frameCell{i}(:,3),'o');
    xlim([0 maxZ])
    ylim([0 maxIntensity])
    title(strcat('Frame: ',num2str(i),'/',num2str(maxFrame)," Time: ",num2str(i*12/60),"min"))
    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v);







