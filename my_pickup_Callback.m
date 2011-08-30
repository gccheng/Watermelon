function my_pickup_Callback( src,evt,handles,width,height )
%MY_PICKUP_CALLBACK Response the click to pick up a position

pos=get(src,'CurrentPoint');
imgPos = [height-pos(2) pos(1)];
hold on; plot(imgPos(2),imgPos(1),'ro'); hold off;

set(handles.posi,'Value', imgPos);
strPosition = sprintf('(%d,%d)', imgPos(1),imgPos(2));
set(handles.posi, 'String', strPosition);

end

