

dat['Year']==2016
plt.close('all')
plt.scatter(dat[dat['Year']==2016]['polstab.in'],dat[dat['Year']==2016]['governance'],color='magenta',marker='*')
plt.savefig('../test.png')
