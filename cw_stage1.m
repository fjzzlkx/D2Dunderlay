%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%   cw_1�������ҳ�һλ��������ͬԪ�ؼ���λ��    %%%%%%%%%%%%%%%%%%%
function [Ulink,U1] = cw_stage1(a,x)
 %a = [3 4 2 3 2 6 4 7 2 5 3 1 4 6 7 4 6];
b = tabulate(a);
c = b(find(b(:,2)>1),1);
c = c';     %������
U1 = setdiff(c,0);
cb = [];
for i = 1:length(c)
    d = [];
    d = find(a == c(i));
    cb = [cb d];
    Ulink = setdiff(cb,x);
end
 %ca�е�һ��Ϊ��ͬ��Ԫ�أ����Ӧ����Ϊ����������λ��  
 %cbΪ������ͬԪ�ص�λ�õļ��ϣ�������