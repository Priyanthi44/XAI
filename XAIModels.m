% startRange = 2;
%  stdc= 1;
%  endRange = 512;
% 
%  table = cell(endRange-startRange,1);
%   for r=startRange:endRange
%     table{r-startRange+1} = norminv((1:r-1)/r,0,stdc);
%   end
%% LOAD DATA and NORMALISE
   % Position data
    A =fopen('position.txt','r');
    posArray = fscanf(A,'%f');
    zmean=mean(posArray);
    zstd=std(posArray);
    posArray = (posArray-zmean)/zstd;
    
%DNN data
  A =fopen('nfeedback.txt','r');
    valArray = fscanf(A,'%f');
    zmean=mean(valArray);
    zstd=std(valArray);
    valArray = (valArray-zmean)/zstd;
    
     % Random data
    A =fopen('random.txt','r');
    rdArray = fscanf(A,'%f');
    rdArray = rdArray(size(rdArray)-size(valArray,1)+1:end);
    zmean=mean(rdArray);
    zstd=std(rdArray);
    rdArray = (rdArray-zmean)/zstd;
    
 
%% CREATE PAA

  w =32;
  p=1;
  %q=size(valArray,1);
%   E_x=getEsumx(valArray,q)-getEsumx(valArray,p);
%   E_xx=getEsumxx(valArray,q)-getEsumxx(valArray,p);
  L=size(valArray,1);
%   mu=E_x/L;
%   sigma=sqrt((E_xx-E_x^2/L)/L-1);
mu_val=mean(valArray);
mu_rd =mean(rdArray);
mu_pos=mean(posArray);

sigma_val=std(valArray);
sigma_rd=std(rdArray);
sigma_pos=std(posArray);
  segments=fix(L/w);
  paa_val=zeros(1,segments);
  paa_rd=zeros(1,segments);
  paa_pos=zeros(1,segments);
  start=1;
  p_end=w-1;
  
  for i=1:segments
%       paa(i)=((getEsumx(valArray,p_end)-(getEsumx(valArray,start))/(L/w))-mu)/sigma;
paa_val(i)=(getSum(start,p_end,valArray)/w-mu_val)/sigma_val;
paa_rd(i)=(getSum(start,p_end,rdArray)/w-mu_rd)/sigma_rd;
paa_pos(i)=(getSum(start,p_end,posArray)/w-mu_pos)/sigma_pos;
      start=p_end+1;
      p_end=start+w-1;
  end
%   max(paa_val)
%   min(paa_val)
%     max(paa_rd)
%   min(paa_rd)
%     max(paa_pos)
%     min(paa_pos)
%   figure(1);
%        [counts, bins] = hist(paa);
%    plot(bins, counts);
%    t_vect=1:1:size(paa,2);
%    figure(2);
%  plot(t_vect,paa);
%% CREATE SAX

isax_val= convertPaatoSax(paa_val);
isax_rd= convertPaatoSax(paa_rd);
isax_pos= convertPaatoSax(paa_pos);

isax_val=char(isax_val);
isax_rd=char(isax_rd);
isax_pos=char(isax_pos);

Len=8;
wordlen=Len;

valLinkedlist= strings(fix(size(isax_val,2)/wordlen)-1,1);
rdLinkedlist= strings(fix(size(isax_rd,2)/wordlen)-1,1);
posLinkedlist= strings(fix(size(isax_pos,2)/wordlen)-1,1);

for i=1:fix(size(isax_val,2)/wordlen)-2
    valLinkedlist(i)=string(isax_val(i:i+wordlen-1));
    rdLinkedlist(i)=string(isax_rd(i:i+wordlen-1));
    posLinkedlist(i)=string(isax_pos(i:i+wordlen-1));
end

valLinkedlist=categorical(valLinkedlist);
posLinkedlist=categorical(posLinkedlist);
rdLinkedlist=categorical(rdLinkedlist);
%% IDENTIFY MOTIFS OF DNN

figure
h1 = histogram(valLinkedlist);
 xlabel("Words")
 ylabel("Frequency")
title("DNN Distribution")
valCounts = h1.BinCounts;
valNames = h1.Categories;
RGBval=zeros(6); 
RGBrd=zeros(6);
RGBpos=zeros(6); 
setResultant=[0 0 0];
val = max(valCounts);
maxval=max(val);
while(val>3)
    

idxMaxCounts =valCounts == val ;
frequentClasses = valNames(idxMaxCounts);

if ~isempty(frequentClasses)
    classes=1;
    k=1;
%% Remove adjecent classes and indicate the new word length by k

    if size(frequentClasses,2)>1
     idxfrequent = ismember(valLinkedlist,frequentClasses);
     [classes,k]=getClasses(idxfrequent,frequentClasses);
    end
%    for j=1:size(k)
%       if k(j)>0
%           %Check whether the adjecent classes were motifs
%           idx= ismember(valLinkedlist,frequentClasses(j)); 
%           potentialmofits= valLinkedlist(idx+1)
%           
%       end
%   end
% for i=1:Classes

  for i=1:size(frequentClasses,2)
     
      if(classes(i))
        idxfrequent = ismember(valLinkedlist,frequentClasses(i));
        [motifpos]=identifyMotif(idxfrequent,val);
        if(k(i)>1)
            for j=1:k(i)
        if (compareElements(valLinkedlist(motifpos+j),1)<8)
            
            k(i)=k(i)-j;
            classes(i+j)=1;
            if(k(i)==1)
                break;
            end
        end
            end
        end
        
        %% TODO eliminate overlaps
        [wordlenbefore,wordlenafter]=expandMotif(motifpos,k(i),Len,valLinkedlist);
        
            type='\Random\';
          RGBrd=compareSequences(rdLinkedlist,motifpos,k(i),wordlenbefore,wordlenafter,type,  [num2str(val,'%02d') '-' num2str(i,'%02d')],Len);%,RGBrd*(val/maxval));    
          % Get the resultant force
          %% TODO multiply with constant matrix to represent the difference between a-b and a-f
          constMatrix = [1   1.1 1.2 1.3 1.4 1.5;
                         1.1 1   1.1 1.2 1.3 1.4;
                         1.2 1.1 1  1.1 1.2 1.3;
                         1.3 1.2 1.1 1  1.1 1.2;
                         1.4 1.3 1.2 1.1 1 1.1; 
                         1.5 1.4 1.3 1.2 1.1 1];
                     RGBrd=RGBrd*constMatrix;
             
          rdforce=resultantForce(RGBrd);
          if rdforce<0
           type=['\Values\neg\' num2str(abs(rdforce),'%02d') '\'];   
          else
             type=['\Values\pos\'  num2str(abs(rdforce),'%02d') '\']; 
          end
          %Categorise the rest depending on that
          

          RGBval=compareSequences(valLinkedlist,motifpos,k(i),wordlenbefore,wordlenafter,type , [num2str(val,'%02d') '-' num2str(i,'%02d')],Len);%,RGBval*(val/maxval));
           
          RGBval=RGBval*constMatrix;
             
          valforce=resultantForce(RGBval);
          if rdforce<0
           type=['\Position\neg\'  num2str(abs(rdforce),'%02d') '\'];   
          else
                 type=['\Position\pos\'  num2str(abs(rdforce),'%02d') '\'];
          end

      
          RGBpos=compareSequences(posLinkedlist,motifpos,k(i),wordlenbefore,wordlenafter,type,  [num2str(val,'%02d') '-' num2str(i,'%02d')],Len);%,RGBpos*(val/maxval));
           
          RGBpos=RGBpos*constMatrix;
             
          posforce=resultantForce(RGBpos);
       
       setResultant =[setResultant ; [rdforce valforce posforce]];
       
      end

   end

end
    idxMaxCounts=0;
    val=val-1;
end
% image(setResultant,'CDataMapping','scaled');
% colorbar
% figure;
% plot(setResultant(:,1),setResultant(:,2),'r.');
% axis equal;
% title 'Random Environment vs DNN response';
% xlabel('Random Distribution')
% ylabel('DNN response')
% figure;
% plot(setResultant(:,1),setResultant(:,3),'r.');
% axis equal;
% title 'Random Environment vs Position';
% xlabel('Random Distribution')
% ylabel('Position')
% figure;
% plot(setResultant(:,2),setResultant(:,3), 'r.');
% axis equal;
% title 'DNN response vs Position';
% xlabel('DNN response')
% ylabel('Position')
% figure;
% plot3(setResultant(:,1),setResultant(:,2),setResultant(:,3),'.');
% axis equal;
% title 'Resultant Distribution';
% xlabel('Random Distribution')
% ylabel('DNN response')
% zlabel('Position')

% hold on;
% plot(setResultant(:,1),setResultant(:,3),'r.');
% axis equal;
% title 'Resultant Distribution';
% xlabel('Random Distribution')
% ylabel('DNN response')
%zlabel('Position')
opts = statset('Display','final');
[idx,C] = kmeans(setResultant,3,'Distance','cityblock',...
    'Replicates',50,'Options',opts);
rdCluster1=setResultant(idx==1,1);
rdCluster2=setResultant(idx==2,1);
  rdCluster3=setResultant(idx==3,1);
%  rdCluster4=setResultant(idx==4,1);
valCluster1=setResultant(idx==1,2);
valCluster2=setResultant(idx==2,2);
  valCluster3=setResultant(idx==3,2);
%  valCluster4=setResultant(idx==4,2);
posCluster1=setResultant(idx==1,3);
posCluster2=setResultant(idx==2,3);
 posCluster3=setResultant(idx==3,3);
% posCluster4=setResultant(idx==4,3);


c1Spread= sqrt((rdCluster1-C(1,1)).^2+(valCluster1-C(1,2)).^2+(posCluster1-C(1,3)).^2);
c2Spread= sqrt((rdCluster2-C(2,1)).^2+(valCluster2-C(2,2)).^2+(posCluster2-C(2,3)).^2);
 c3Spread= sqrt((rdCluster3-C(3,1)).^2+(valCluster3-C(3,2)).^2+(posCluster3-C(3,3)).^2);
% c4Spread= sqrt((rdCluster4-C(4,1)).^2+(valCluster4-C(4,2)).^2+(posCluster4-C(4,3)).^2);
c1Mean =mean(c1Spread);
c2Mean=mean(c2Spread);
 c3Mean= mean(c3Spread);
% c4Mean= mean(c4Spread);
c1std=std(c1Spread);
c2std=std(c2Spread);
 c3std=std(c3Spread);
% c4std=std(c4Spread);

figure;
plot3(setResultant(idx==1,1),setResultant(idx==1,2),setResultant(idx==1,3),'r.','MarkerSize',12)
hold on
plot3(setResultant(idx==2,1),setResultant(idx==2,2),setResultant(idx==2,3),'b.','MarkerSize',12)
hold on
plot3(setResultant(idx==3,1),setResultant(idx==3,2),setResultant(idx==3,3),'g.','MarkerSize',12)
  hold on
%  plot3(setResultant(idx==4,1),setResultant(idx==4,2),setResultant(idx==4,3),'m.','MarkerSize',12)
%  hold on
% plot(setResultant(idx==1,1),setResultant(idx==1,2),'y.','MarkerSize',12)
%  hold on
%  plot(setResultant(idx==2,1),setResultant(idx==2,2),'m.','MarkerSize',12)
plot3(C(1,1),C(1,2),C(1,3),'kx','MarkerSize',15,'LineWidth',3) 
plot3(C(2,1),C(2,3),C(2,3),'kx','MarkerSize',15,'LineWidth',3) 
 plot3(C(3,1),C(3,2),C(3,3),'kx','MarkerSize',15,'LineWidth',3) 
% plot3(C(4,1),C(4,2),C(4,3),'kx','MarkerSize',15,'LineWidth',3) 
legend('Cluster 1','Cluster 2','Cluster 3','Centroids','Location','NW')
%legend('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Centroids','Location','NW')
title 'Cluster Assignments and Centroids'
hold off
 figure;
      [counts, bins] = hist(setResultant(:,3));
    
 plot(bins, counts);

% surf([setResultant(idx==2,1),setResultant(idx==2,2),setResultant(idx==2,3)])
% hold on 
% image([setResultant(idx==2,1),setResultant(idx==2,2),setResultant(idx==2,3)],'CDataMapping','scaled')
% hold off;
% surf([setResultant(idx==1,1),setResultant(idx==1,2),setResultant(idx==1,3)])
% hold on 
% image([setResultant(idx==1,1),setResultant(idx==1,2),setResultant(idx==1,3)],'CDataMapping','scaled')
% hold off;
% surf([setResultant(idx==3,1),setResultant(idx==3,2),setResultant(idx==3,3)])
% hold on 
% image([setResultant(idx==3,1),setResultant(idx==3,2),setResultant(idx==3,3)],'CDataMapping','scaled')

%% ANALYSE DRIFT
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
force=round(force,3);
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

%% COMPARE SUBSEQUENCES
% function force= resultantForce(A,weight,rd)
% even=false;
% multiplier=2;
% for i=1:size(A,1)
%     for j=1:size(A,2)
%         if (j>i)
%             A(i,j)=-A(i,j);
%         end
%     end
% end
% forceSum= sum(A,2);
% for i=1:numel(forceSum)
%     forceSum(i)=weight*forceSum(i);
%     if(rd)
%     if(even)
%         forceSum(i)=-forceSum(i);
%     end
%     if(i==(numel(forceSum)/2))
%         even=true;
%         multiplier=1/2;
%     else
%         weight=weight*multiplier;
%     end
%     else
%          weight=weight*multiplier;
%     end
%         
% end
% force=sum(forceSum);
% end

function RGB=compareSequences(Linkedlist,motifpos,k,wordlenbefore,wordlenafter,type,name,Len)
before=char(Linkedlist(motifpos-1));
after=char(Linkedlist(motifpos+k-1+1));
arr=char(Linkedlist(motifpos));
if(k>1)
for i=2:k
   
        arr=strcat(arr,char(Linkedlist(motifpos+i)));  
end
end
list =strcat(before(Len-wordlenbefore+1:Len),arr,after(1:wordlenafter));



% compare rdList items with each other  and posList items with each other and identify motifs if any  
map=bitMap(list);
% map =squareform(pdist(map));
%% uncheck******uncomm
%image(map,'CDataMapping','scaled','Tag',name);
%colorbar
%%
destinationFolder = [pwd '\maps' type ];
if ~exist(destinationFolder, 'dir')
  mkdir(destinationFolder);
else
%     if(~strcmp(type,'\Random\'))
%    destinationFolder 
%     end
end
% if (strcmp(type,'\Values\neg\56\'))
%     destinationFolder
% end
%% uncheck******uncom
%save([destinationFolder name '.mat'], 'map')
RGB = map;%imread([destinationFolder name '.png']);
%%
% RGB =RGB/max(max(abs(RGB)));
% RGB=RGB+uint8(rgb);

end
%% EXPAND THE MOTIF

function [wordlenbefore,wordlenafter]=expandMotif(motifpos,k,Len,valLinkedlist)

%move forward the motif until a repeating pattern is not found 
if ~(k>1)
    k=0;
end
    
 wordlenbefore = Len-compareElements(valLinkedlist(motifpos-1),0);
 wordlenafter = compareElements(valLinkedlist(motifpos+k+1),1);

    
   
end
%% GET CLASSES

function [classes,k]=getClasses(idxfrequent,frequentClasses)
count=1;
match=false;
q=0;
cls=size(frequentClasses,2);
k=ones(cls,1); % holding information about the wordlength
classes=ones(cls,1); % holding information about the classes
p=0; % holding position of the classes array
for i=1:size(idxfrequent)
    if idxfrequent(i)==1
        j=1;
        p=p+1; 
            
           if idxfrequent(i+1)==1 %if next to eachother
%                className=valLinkedlist(idxfrequent(i+1));
%                idx=ismember(className,valLinkedlist);
              while(j<cls)
               if(~match)
                   q=i;
               count=p;
               match=true;
               end
               if idxfrequent(q+1)==1                 
               idxfrequent(q+1)=0;
               classes(p+1)=0;
               k(count)=k(count)+1;
               p=p+1;
               j=j+1;
               q=q+1;
               else
                   p=p-1;
                   break;
               end
              end
           else
               match=false;               
             
           end
           
    else
       if(p==cls)
           break;
       end
    end
        
     
        
end
end
%% IDENTIFY MOTIFS

function [motifpos]=identifyMotif(idxfrequent,val)
count=1;
motifpos=zeros(val,1);
for i=1:size(idxfrequent)
    if idxfrequent(i)==1
        motifpos(count)=i;
        count=count+1;
    end
    if(count==val+1)
        break;
    end
end
end
%% IDENTIFY MATCHING ELEMENTS IN ADJECENT ITEMS

function k= compareElements(List,order)
match=true;

len=strlength(char(List(1)));

if(order)
    k=1;
    
while(k<=len)
for j=1:size(List)-1
    
    el1=List(1);
    el2=List(j+1);
    
    
         el1=char(el1);
         el2=char(el2);
        
             if(el1(k)==el2(k))
                 match=true;
             else
                 match=false;
                 
             end
            
        
 
 
    if(~match)
        break;
    end
    
end
  if(~match)
      k=k-1;
      break;
        
   end
         k=k+1;
     
   
end


else
     k=len;
  
 while(k>=1)
    
for j=1:size(List)-1
    
    el1=List(1);
    el2=List(j+1);

    
    
         el1=char(el1);
         el2=char(el2);
        
             if(el1(k)==el2(k))
                 match=true;
             else
                 match=false;
             end
            
        
 
 
    if(~match)
        break;
    end
     
end
    if(~match)
        break;
    end
    k=k-1;
    
  
 end 
end

end
%% CONVERT PAA TO SAX

function isax= convertPaatoSax(paa)
isax=zeros(1,size(paa,2));
r=6;
range=norminv((1:r-1)/r,0,1);
for i=1:size(paa,2)
        if paa(i)>range(5)
        isax(i)='f';
        elseif paa(i)>range(4)
        isax(i)='e';
        elseif paa(i)>range(3)
        isax(i)='d';
        elseif paa(i)>range(2)
        isax(i)='c';
        elseif paa(i)>range(1)
        isax(i)='b';
        else 
        isax(i)='a';
        end
        
end
end
%% CREATE MATRIX

function bmap= bitMap(L)
bmap=zeros(6,6);
s=size(L);
List=L';
len=s(1)*s(2); % s1-number of sequences, s2-sequence length

for i=1:len
    
    if(mod(i,s(2))) % inside the sequence
       
               el=strcat(List(i),List(i+1));
               x=getIndex(el);
               bmap(x(1),x(2))=bmap(x(1),x(2))+1; 
%                if(x(2)<x(1))
%                    bmap(x(1),x(2))=-bmap(x(1),x(2));
%                end
            
    end
        
end

% nexteven=0;
% if(mod(s(2),2)) % odd number
%     nexteven=1;
% end
% if(~nexteven) %even number
%     for i=1:len
%     
%     if(mod(i,2))
%         
%         el=strcat(List(i),List(i+1));
%         x=getIndex(el);
%         bmap(x(1),x(2))=bmap(x(1),x(2))+1;
%     end
%        
%         
%     end
% else % odd number
% for i=1:len
%     
%     if(mod(i,s(2))) % inside the sequence
%         if(i>s(2))
%             if(nexteven)&&(~mod(i,2))
%                el=strcat(List(i),List(i+1));
%                x=getIndex(el);
%                bmap(x(1),x(2))=bmap(x(1),x(2))+1; 
%             elseif(mod(i,2))   
%                 el=strcat(List(i),List(i+1));
%                 x=getIndex(el);
%                 bmap(x(1),x(2))=bmap(x(1),x(2))+1;
%             end  
%         
%         elseif(mod(i,2))  % right before the end of sequence 
%         el=strcat(List(i),List(i+1));
%         x=getIndex(el);
%         bmap(x(1),x(2))=bmap(x(1),x(2))+1;
%         end
%       
%     else
%         el=strcat(List(i-1),List(i));
%         x=getIndex(el);
%         bmap(x(1),x(2))=bmap(x(1),x(2))+1;
%         if(mod(i/s(1),2))
%             nexteven=1;
%             
%         else
%             nexteven=0;
%             
%         end
%         
%     end
        
% end
% end

end
%% MATCH TWO LETTERS

function indx=getIndex(el)
indx= [0,0];
switch el
    case 'aa'
        indx=[1,1];
    case 'ab'
        indx=[1,2];
    case 'ac'
        indx= [1,3];
    case 'ad'
        indx=[1,4];
    case 'ae'
        indx=[1,5];
    case 'af'
        indx= [1,6];
     case 'ba'
        indx=[2,1];
    case 'bb'
        indx=[2,2];
    case 'bc'
        indx= [2,3];
    case 'bd'
        indx=[2,4];
    case 'be'
        indx=[2,5];
    case 'bf'
        indx= [2,6];
    case 'ca'
        indx=[3,1];
    case 'cb'
        indx=[3,2];
    case 'cc'
        indx= [3,3];
    case 'cd'
        indx=[3,4];
    case 'ce'
        indx=[3,5];
    case 'cf'
        indx= [3,6];  
          case 'da'
        indx=[4,1];
    case 'db'
        indx=[4,2];
    case 'dc'
        indx= [4,3];
    case 'dd'
        indx=[4,4];
    case 'de'
        indx=[4,5];
    case 'df'
        indx= [4,6];
    case 'ea'
        indx=[5,1];
    case 'eb'
        indx=[5,2];
    case 'ec'
        indx= [5,3];
    case 'ed'
        indx=[5,4];
    case 'ee'
        indx=[5,5];
    case 'ef'
        indx= [5,6];
      case 'fa'
        indx=[6,1];
    case 'fb'
        indx=[6,2];
    case 'fc'
        indx= [6,3];
    case 'fd'
        indx=[6,4];
    case 'fe'
        indx=[6,5];
    case 'ff'
        indx= [6,6];
    
end
end
 function sum=getSum(istart,iend, array)
 sum=0;
 for i=istart:iend
     sum=sum+ array(i);
 end
 end
  function Esum_x =getEsumx(valArray,iend)
  Esum_x=0;
  for i=1:iend
      Esum_x=Esum_x+valArray(i);
  end
  end 
  function Esum_xx =getEsumxx(valArray,iend)
  Esum_xx=0;
  for i=1:iend
      Esum_xx=Esum_xx+valArray(i)^2;
  end
  end