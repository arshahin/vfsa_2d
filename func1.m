function func=func1(x,y)
func=sign(sin(x)./x).*abs(sin(x)./x).^0.25.*sign(sin(y)./y).*abs(sin(y)./y).^0.25;
