%A = [0 0 1 2 2 1; 0 0 0 1 0 0 ;3 3 0 0 0 0; 0 0 0 0 1 1;1 1 1 1 0 0; 0 1 1 1 0 0];
RGB=zeros(6,6);
x=0;
for i=1:6
    for j=1:6
        A=zeros(6,6);
        A(i,j)=1;
        RGB(i,j)=resultantForce(A);
%         x=A(i,j);
    end
end
%resultantForce(A)
h=surf(RGB,'FaceAlpha',0.5)
% set(h,'edgecolor','r','facecolor',[.1 .9 .1])

% colormap summer
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
        if (j>i)
            A(i,j)=-A(i,j);
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
