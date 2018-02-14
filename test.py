import os,sys,glob,time,collections,gc,calendar
from subprocess import Popen

def foo(z=10,**kwargs):
    print z

def asdas(x,y,**kwargs):
    print x,y
    print kwargs
    foo(**kwargs)


asdas(2,3,rr=100)
