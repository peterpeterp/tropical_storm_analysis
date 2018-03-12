
self=found_tcs
postions=self._detected.values
used_pos=[]
for p in postions[postions[:,-1].argsort()[::-1],:].tolist():
    if p not in used_pos:
        track=[p]

        running=True
        #go backwards
        while True:
            p=track[0]
            if len(track)==1:
                y_exp,x_exp=p[1],p[2]
            else:
                y_exp,x_exp=2*p[1]-track[1][1],2*p[2]-track[1][2]
            candidates=[]
            for p_1 in postions[postions[:,0]==p[0]-1,:].tolist():
                if ((y_exp-p_1[1])**2+(x_exp-p_1[2])**2)**0.5<self._search_radius:
                    candidates.append(p_1)
                    end=False
            if len(candidates)>0:
                track=[candidates[np.array(candidates)[:,-2].argmin()]]+track
            else:
                break

        #go forewards
        while True:
            p=track[-1]
            if len(track)==1:
                y_exp,x_exp=p[1],p[2]
            else:
                y_exp,x_exp=2*p[1]-track[-2][1],2*p[2]-track[-2][2]
            candidates=[]
            for p_1 in postions[postions[:,0]==p[0]+1,:].tolist():
                if ((y_exp-p_1[1])**2+(x_exp-p_1[2])**2)**0.5<self._search_radius:
                    candidates.append(p_1)
                    end=False
            if len(candidates)>0:
                track=track+[candidates[np.array(candidates)[:,-2].argmin()]]
            else:
                break

        asdas



a=[2,3,4,6,7,8,9,10,23,24]
b=np.diff(a,1)==1
def consecutive_sequence(in):
    i,su,out=0,0,[]
    while i <len(b):
        if b[i]:
            su+=1
        else:
            out.append((i-su,su))
            su=0
        i+=1
    out.append((i-su,su))
    return out

d=[sum(g) for b, g in itertools.groupby(b) if b]





#
