%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%   cw函数，找出一位数组中相同元素及其位置    %%%%%%%%%%%%%%%%%%%
function [x] = cw_stage2(a)
 %a = [3 4 2 3 2 6 4 7 2 5 3 1 4 6 7 4 6];
A = unique(a);  %所有被分配的时延
b = tabulate(a);
c = b(find(b(:,2)>1),1);
c = c';     %行向量
cb = [];
for i = 1:length(c)
    d = [];
    d = find(a == c(i));
    cb = [cb d];
    x = setdiff(A,c);   %只分配1个link的时延
end
 %ca中第一列为相同的元素，其对应的行为他们所处的位置  
 %cb为所有相同元素的位置的集合，行向量