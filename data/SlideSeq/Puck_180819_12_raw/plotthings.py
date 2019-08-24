import matplotlib.pyplot as plt
import pandas as pd


pos = pd.read_csv("BeadLocationsForR.csv", index_col=0)
posf = pd.read_csv("BeadLocationsForR_filtered.csv", index_col=0)
pos2 = pd.read_csv("BeadMapping_10-17_1014/BeadLocationsForR.csv", index_col=0)
pos3 = pd.read_csv("BeadMapping_10-17_1014/BeadLocationsForRNoNuclei.csv", index_col=0)
pos4 = pd.read_csv("BeadMapping_10-17_1014/BeadLocationsForRUncropped.csv", index_col=0)

dge = pd.read_csv("MappedDGEForR.csv", index_col=0)

dge.shape


plt.figure()
plt.plot(pos.xcoord, pos.ycoord, 'o', ms=.1, alpha=0.5)
plt.show()

pos3.index
pos2.index
pos3.index & pos.index

pos2.index

(pos2.index & pos4.index).size
posf.index

pos4.shape
(pos4.index & pos.index).size

# posf & pos = 15
# pos2 is subset of pos
# pos2 is subset of pos4
# pos3 is subset of pos2
# pos4 is subset of pos (though almost all)

ii = pos3.index[0]
pos4.loc[ii]
pos2.loc[ii]
pos3.loc[ii]
pos.loc[ii]

# Ok, so everything matches EXCEPT that clearlly the _filtered dataset in the root directory is from a different experiment here
plt.figure()
plt.plot(pos2.xcoord, pos2.ycoord, 'o', ms=.5, alpha=0.5)
plt.show()

plt.figure()
plt.scatter(pos2.ycoord, pos2.xcoord*-1, c=dge.loc['Pcp4'][pos2.index], s=1, cmap='gray_r')
plt.show()

# OK, so the right data to use is the pos2 data.
# Then we swap x'=y and y'=-x to match the slides in the paper
# Note that this does not alter cell-cell distances
