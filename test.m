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
