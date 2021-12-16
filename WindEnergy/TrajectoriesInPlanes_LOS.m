clc,clear all, close all

ref_plane  = 250;         % refence plane for trajectory planning Cartesian [m]
req_planes = 50:50:1500;  % requested planes to calculate the Cartesian points based on LOS [m]

trajectory= [-55,-55,0,0,55,55 ; 55,-55,-55,55,55,-55]; % 6 point rectangular grid  first line is y direction second is z. Origin is (0,0)


% loop over trajectory to find LOS angles
for i=1:length(trajectory)    
   angley(i) = atand(trajectory(1,i)/ref_plane) ;
   anglez(i)= atand(trajectory(2,i)/ref_plane) ; 
end

%loop over planes to get points
for i=1:length (req_planes)
    
    iplane= req_planes(i);
    for ii=1:length(trajectory)    
       plane_traj{i}(1,ii)=iplane*tand(angley(ii));
       plane_traj{i}(2,ii)=iplane*tand(anglez(ii));
    end
end

%% Plotting
 figure ,plot(trajectory(1,:),trajectory(2,:),'X','Linewidth',2),grid on
 title (['Trajectory as requested on reference plane' num2str(ref_plane)])

figure
for i=1:length(req_planes)
    plot3(req_planes(i)*ones(length(trajectory(1,:)),1),plane_traj{1, i}(1,:)',plane_traj{1, i}(2,:)','X','Linewidth',2)
    %              plot3(All_patern{1,ii}{1, i}(:,1),All_patern{1, ii}{1, i}(:,2),All_patern{1, ii}{1, i}(:,3))
    
    hold on
end
plot3(0,0,0,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1],...
    'Marker','o','LineStyle','none');
hold off
grid on
set(gca,'YGrid','on', 'FontSize', 14,'YMinorGrid','on','XMinorGrid','on','ZMinorGrid','on')
