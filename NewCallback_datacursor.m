function output_txt = myfunction(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).


pos = get(event_obj,'Position');


output_txt = {['X: ',num2str(pos(1),4)],...
    ['Y: ',num2str(pos(2),4)]};

title(gca, ['(X,Y) = (', num2str(pos(1)), ', ',num2str(pos(2)), ')']);

% If there is a Z-coordinate in the position, display it as well
if length(pos) > 2
    output_txt{end+1} = ['Z: ',num2str(pos(3),4)];
    title(gca, ['(X,Y,Z) = (', num2str(pos(1)), ', ',num2str(pos(2)),', ',num2str(pos(3)), ')']);
end

% plot([0, pos(1)],[0 pos(2)]);
% from start point with direction
% quiver(0,0,pos(1),pos(2),0,'k--','MaxHeadSize',0.5);

