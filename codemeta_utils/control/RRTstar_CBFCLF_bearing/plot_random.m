function plot_random(rnd,fatal_samples)
if ~isempty(rnd)
    c = 'b';
    if ~isempty(fatal_samples)
        flag = (fatal_samples == rnd);
        flag = any(logical(flag(1,:).*flag(2,:)));
        if flag
            c = 'r';
            figure(1)
            plot(rnd(1),rnd(2),'s','MarkerEdge',c,'MarkerFace','r','MarkerSize',9)
            hold on
        end
    end
end