%A = [0 0 1 2 2 1; 0 0 0 1 0 0 ;3 3 0 0 0 0; 0 0 0 0 1 1;1 1 1 1 0 0; 0 1 1 1 0 0];
RGB=zeros(6,6);
x=0;
for i=1:6
    for j=1:6
        A=zeros(6,6);
        A(i,j)=1;
        constMatrix = [1   1.01 1.02 1.03 1.04 1.05;
                         1.01 1   1.01 1.02 1.03 1.04;
                         1.02 1.01 1  1.01 1.02 1.03;
                         1.03 1.02 1.01 1  1.01 1.02;
                         1.04 1.03 1.02 1.01 1 1.01; 
                         1.05 1.04 1.03 1.02 1.01 1];
                     A=A.*constMatrix;
        RGB(i,j)=resultantForce(A);
%         x=A(i,j);
    end
end
%resultantForce(A)
 colormap summer
h=surf(RGB,'FaceAlpha',0.5)
% set(h,'edgecolor','r','facecolor',[.1 .9 .1])


% shading interp
%  hold on 
% image([setResultant(idx==2,1),setResultant(idx==2,2),setResultant(idx==2,3)],'CDataMapping','scaled')
% hold off;
%    image(RGB,'CDataMapping','scaled');
% colorbar

% destinationFolder=[pwd '\maps\Random\'];
% Files=dir(destinationFolder);
% for k=1:length(Files)
%     if Files(k).bytes>0
%         
%    FileName=Files(k).name;
%    [filepath name ext]=fileparts(FileName);
%    if (strcmp(ext,'.mat'))
%        if(itstheFile(name))
%        example = matfile([destinationFolder FileName]);
%     RGB = example.map;
%    image(RGB,'CDataMapping','scaled','Tag',name);
% colorbar
%        end
%  val=resultantForce(RGB);
%  val=round(val,3);
%  if (val>0)
%      type ='pos\'
%  else
%      type= 'neg\'
%  end
%  destinationFolder = [pwd '\maps\' type  num2str(abs(rdforce),'%02d') '\' ];
%  
% if ~exist(destinationFolder, 'dir')
%   mkdir(destinationFolder);
% else
% %     if(~strcmp(type,'\Random\'))
% %    destinationFolder 
% %     end
%  end
% % if (strcmp(type,'\Values\pos\3.535000e+00\'))
% %\Values\pos\5.921000e+00
% %\maps\Values\pos\5.630000e+00
% %     destinationFolder
% % end
% save([destinationFolder val '.mat'], 'map')
%    end
%     end
%     
% end
% destinationFolder=[pwd '\images\Random\'];
% Files=dir(destinationFolder);
% for k=1:length(Files)
%     if Files(k).bytes>0
%         
%    FileName=Files(k).name;
%    [filepath name ext]=fileparts(FileName);
%    if (strcmp(ext,'.mat'))
%        
%        example = matfile([destinationFolder FileName]);
%     RGB = example.map;
% %   image(RGB,'CDataMapping','scaled','Tag',name);
% % colorbar
function yes=itstheFile(name)
yes=0;

if strcmp(name,'04-146')
    yes=1;
elseif strcmp(name,'04-239')
    yes=1;
elseif strcmp(name,'04-249')
    yes=1;
else
    yes=0;
end
end
function force= resultantForce(A)

for i=1:size(A,1)
    for j=1:size(A,2)
       if(A(i,j)>0)
        if (j>i)
            A(i,j)=-A(i,j);
        end
       end
    end
end
rawSum= sum(A,2);
colSum=sum(A);
rawSum= getNorm(rawSum);
colSum=getNorm(colSum);
for i=1:numel(rawSum)
    rawSum(i)=getDrift(i,rawSum(i));
    colSum(i)=getDrift(i,colSum(i));
        
end
force=colSum*rawSum;
end

function X=getNorm(A)
mu=mean(A);
st=std(A);
X=(A-mu)/st;
end
function val=getDrift(indx,drift)
if (indx==1)%a
    if(drift<0)%atob
        val=-0.97-drift;
    else
        val=0;
    end
elseif (indx==2)%b
    if(drift>0)%b toa
        
         val=-0.97-drift;
    else %b to c
       val=-0.43-drift;
    end
elseif (indx==3)%c
    if(drift<0)
        val=0-drift;
    else
        val=-0.43-drift;
    end
elseif (indx==4)%d
    if(drift>0)
        val=0-drift;
    else
        val=0.43-drift;
    end
elseif (indx==5)%e
    if(drift>0)
        val=0.43-drift;
    else
        val=0.97-drift;
    end
elseif (indx==6)%f
    if(drift<0)
         val=0;   
    else
        
        val=0.97- drift;
    end
        
end
end

% [idx,C] = kmeans(setResultant,3,'Distance','cityblock',...
%     'Replicates',5000,'Options',opts);
% valueCluster1=setResultant(idx==1,1);
% valueCluster2=setResultant(idx==2,1);
%  valueCluster3=setResultant(idx==3,1);
% rdCluster1=setResultant(idx==1,2);
% rdCluster2=setResultant(idx==2,2);
%  rdCluster3=setResultant(idx==3,2);
% posCluster1=setResultant(idx==1,3);
% posCluster2=setResultant(idx==2,3);
% posCluster3=setResultant(idx==3,3);
% 
% c1Spread= sqrt((valueCluster1-C(1,1)).^2+(rdCluster1-C(2,1)).^2+(posCluster1-C(3,1)).^2);
% c2Spread= sqrt((valueCluster2-C(1,2)).^2+(rdCluster2-C(2,2)).^2+(posCluster2-C(3,2)).^2);
% c3Spread= sqrt((valueCluster3-C(1,3)).^2+(rdCluster3-C(2,3)).^2+(posCluster3-C(3,3)).^2);
% c1Mean =mean(c1Spread);
% c2Mean=mean(c2Spread);
% c3Mean= mean(c3Spread);
% c1std=std(c1Spread);
% c2std=std(c2Spread);
% c3std=std(c3Spread);
% figure;
% plot3(setResultant(idx==1,1),setResultant(idx==1,2),setResultant(idx==1,3),'r.','MarkerSize',12)
% hold on
% plot3(setResultant(idx==2,1),setResultant(idx==2,2),setResultant(idx==2,3),'b.','MarkerSize',12)
% hold on
% plot3(setResultant(idx==3,1),setResultant(idx==3,2),setResultant(idx==3,3),'g.','MarkerSize',12)
% hold on
% plot(setResultant(idx==4,1),setResultant(idx==4,2),'y.','MarkerSize',12)
%  hold on
%  plot(setResultant(idx==5,1),setResultant(idx==5,2),'m.','MarkerSize',12)
% plot3(C(:,1),C(:,2),C(:,3),'kx',...
%      'MarkerSize',15,'LineWidth',3) 
%  plot3(C(:,1),C(:,2),C(:,3),'kx',...
%      'MarkerSize',15,'LineWidth',3) 
%  plot3(C(:,1),C(:,2),C(:,3),'kx',...
%      'MarkerSize',15,'LineWidth',3) 
%  plot3(C(:,1),C(:,2),C(:,3),'kx',...
%      'MarkerSize',15,'LineWidth',3) 
% legend('Cluster 1','Cluster 2','Cluster 3','Centroids',...
%        'Location','NW')
% title 'Cluster Assignments and Centroids'
% hold off
