function E=edgesSymmetrize(E)
E=unique([E;fliplr(E)],'rows');
