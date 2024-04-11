import sys
path = sys.argv[1]
from sklearn import metrics
import numpy as np
a = np.loadtxt('temp.txt',dtype=int,delimiter=' ')
M = a.T
labels_true = M[0,:]
labels_pred = M[1,:]
with open(path, "a") as filewrite:
  filewrite.write("调和兰德指数ARI:"+str(metrics.adjusted_rand_score(labels_true, labels_pred))+"\n")
  filewrite.write("调整互信息AMI:"+str(metrics.adjusted_mutual_info_score(labels_true, labels_pred))+"\n")
  filewrite.write("V测量:"+str(metrics.v_measure_score(labels_true, labels_pred))+"\n")
  filewrite.write("福尔克斯-锦葵指数FMI:"+str(metrics.fowlkes_mallows_score(labels_true, labels_pred))+"\n")
filewrite.close()
quit
