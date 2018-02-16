import itertools

def path_length(points):
    return np.sum(np.sqrt(np.sum(np.diff(np.array(points), axis=0)**2, axis=1)))


track=backup[:,:]
ordered_track=[]

segments=[]
permu=[]
permu_n=[]
for tt in set(track[:,0]):
    print tt
    permu.append(list(set(permutations(range(len(np.where(track[:,0]==tt)[0]))))))
    segments.append(track[track[:,0]==tt])
    permu_n.append(len(permu[-1]))

best_permu=(0,0)
length=path_length(np.concatenate((segments[0][permu[0][best_permu[0]],1:3],segments[1][permu[1][best_permu[1]],1:3])))
for per1,per2 in itertools.product(range(permu_n[0]),range(permu_n[1])):
    length_new=path_length(np.concatenate((segments[0][permu[0][per1],1:3],segments[1][permu[1][per2],1:3])))
    if length>length_new:
        best_permu=(per1,per2)
        length=length_new

ordered_track=np.concatenate((segments[0][permu[0][best_permu[0]],:],segments[1][permu[1][best_permu[1]],:]))
start_point=ordered_track[-1,1:3]

for step in range(len(permu_n))[2:]:
    best_permu=0
    length=path_length(np.vstack((start_point,segments[step][permu[step][best_permu],1:3])))
    for per in range(permu_n[step]):
        length_new=path_length(np.vstack((start_point,segments[step][permu[step][per],1:3])))
        if length>length_new:
            best_permu=per
            length=length_new

    ordered_track=np.concatenate((ordered_track,segments[step][permu[step][best_permu],:]))




    #asa
