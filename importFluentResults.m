function [x,y,u,v,p] = importFluentResults(filename)
    lines = readlines(filename);
    lines_split = split(lines(2:end-1),',',2);

    v = double(lines_split(:,6));
    u = double(lines_split(:,5));
    p = double(lines_split(:,4));
    y = double(lines_split(:,3));
    x = double(lines_split(:,2));
end
