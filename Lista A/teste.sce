deff('[y]=test0(x)','y=x+x^2+sin(x*2*%pi)')
deff('[y]=test1(x)','y=-x+x^2+x^3')
deff('[y]=test2(x)','y=sqrt(x)')
x=-2:0.5:3;
a=1;
b=0;
t1=(a==1);
t2=(b>0.5);
if and([t1 t2]) then
 y=tert0(x);
 elseif or([t1 t2]) then
 y=test1(x);
 else
 y=test2(x);
end,
plot2d(x,y,-3)
set("current_figure",1)
xset('mark size', 2)
plot2d(x,y,-3)
set("current_figure",2)
xset('mark size', 4)
plot2d(x,y,-3)
set("current_figure",3)
xset('mark size', 5)
plot2d(x,y,-3)
