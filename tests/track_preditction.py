import os,sys,glob,time,collections,gc,calendar
from datetime import datetime as datetime
import numpy as np
import dimarray as da
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.path import Path
from matplotlib.collections import PatchCollection
import matplotlib.patches as patches
import seaborn as sn
import webcolors

def projected_polygon(map,points):
    verts=[]
    for lon,lat in points:
        verts.append(m(lon,lat))

    codes2 = [Path.MOVETO,
             Path.LINETO,
             Path.LINETO,
             Path.LINETO,
             Path.CLOSEPOLY,
             ]

    path = Path(verts, codes2)
    return path

#https://stackoverflow.com/questions/31681373/making-svm-run-faster-in-python
from sklearn.ensemble import BaggingClassifier, RandomForestClassifier
from sklearn import datasets
from sklearn.multiclass import OneVsRestClassifier
from sklearn.svm import SVC

os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')

X_all=np.load('detection/training_1979_X.npy')
y_all=np.load('detection/training_1979_Y.npy')

indices_tc=np.where(y_all!=0)[0]
indices_no=np.where(y_all==0)[0]
indices=np.concatenate((np.where(y_all!=0)[0],indices_no[np.random.randint(len(indices_no), size=100000)]))

y=y_all[indices]
X=X_all[indices,:]

start = time.time()
clf = RandomForestClassifier(min_samples_leaf=20,class_weight={1:0.9,0:0.4})
clf.fit(X, y)
end = time.time()
print "Random Forest", end - start, clf.score(X,y)
proba = clf.predict_proba(X)

n_estimators = 10
start = time.time()
clf = OneVsRestClassifier(BaggingClassifier(SVC(kernel='linear', probability=True, class_weight={1:0.9,0:0.4}), max_samples=1.0 / n_estimators, n_estimators=n_estimators))
clf.fit(X, y)
end = time.time()
print "Bagging SVC", end - start, clf.score(X,y)
proba = clf.predict_proba(X)



indices=np.concatenate((np.where(y_all!=0)[0],indices_no[np.random.randint(len(indices_no), size=1000)]))
X=X_all[indices,:]

xx=X.reshape(X.shape[0],116/4,4)
xx[xx==-99]=np.nan

fig,axes=plt.subplots(nrows=5,ncols=3)
axes=axes.flatten()
for ax,var in zip(axes,[0,2,6,7,9,13,14,16,20,21,23,27,28,0,0]):
    for i in range(4):
        ax.scatter(range(xx.shape[0]),xx[:,var,i],marker='^')
        ax.get_xaxis().set_visible(False)

plt.savefig('plots/test.png')







#sad
