% definiram funkcijo mcc_pi
function [ktkrog,ktkvad] = mcc_pi(st_tock)

% x-koordinate naključnih točk
x = 2 * rand(st_tock,1) - 1;
y = 2 * rand(st_tock,1) - 1;
% logični vektor točk znotraj kroga
not_logic = x.^2 + y.^2 <= 1;% logični vektor točk znotraj kroga
% koordinate točk znotraj kroga in znotraj kvadrata
ktkrog = [x(not_logic),y(not_logic)];
ktkvad = [x,y];

end