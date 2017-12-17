function  outputstring=plotmesh(Nodes,Elements),

Nel = size(Elements,1);

clf
plot3(Nodes(:,1),Nodes(:,2),Nodes(:,3),'g*')
for el=1:Nel,
   type   = Elements(el,1);
   node1  = Nodes(Elements(el,3),:);
   node2  = Nodes(Elements(el,4),:);
   hold on
   if type==1,
      plot3([node1(1);node2(1)],[node1(2);node2(2)],[node1(3);node2(3)],'r');
   elseif type==2,
      plot3([node1(1);node2(1)],[node1(2);node2(2)],[node1(3);node2(3)],'b');
   end
end
axis equal
view(-37.5,15)