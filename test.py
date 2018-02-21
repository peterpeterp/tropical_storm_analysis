

def circle(y,x,radius):
    y_min=max(0,y-radius)
    y_max=min(len(U.lat),y+radius+1)
    x_min=max(0,x-radius)
    x_max=min(len(U.lon),x+radius+1)
    y_,x_=[],[]
    for j in [x_min,x_max]:
        for i in range(y_min,y_max+1):
            y_.append(i)
            x_.append(j)
    for i in [y_min,y_max]:
        for j in range(x_min,x_max+1):
            y_.append(i)
            x_.append(j)
    return y_,x_

y,x=circle(10,12,4)
a=U.ix[20,0,y,x]
