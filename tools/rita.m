function slask=rita(bild,st)
if nargin == 1,
 st='-';
end;
if size(bild)==0,
  slask=[];
else
  plot(bild(1,:),bild(2,:),st,'MarkerSize',10,'LineWidth',5);
  slask=[];
end;
