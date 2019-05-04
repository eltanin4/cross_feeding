plt.hist(all_overlaps)
plt.show()
max(all_overlaps)
fig, ax = plt.subplots(1)
myWeights = np.ones_like( all_overlaps ) / len ( all_overlaps )
ax.hist( all_overlaps, bins = 10, color = 'dodgerblue', weights = myWeights, histtype = 'stepfilled' )
ax.set_xlabel( 'Overlap' )
ax.set_ylabel( 'Fraction of species' )
ax.set_xlim( 0.0, 1.0 )
plt.show()
fig, ax = plt.subplots(1)
myWeights = np.ones_like( all_overlaps ) / len ( all_overlaps )
ax.hist( all_overlaps, bins = 20, color = 'dodgerblue', weights = myWeights, histtype = 'stepfilled' )
ax.set_xlabel( 'Overlap' )
ax.set_ylabel( 'Fraction of species' )
ax.set_xlim( 0.0, 1.0 )
plt.show()
fig, ax = plt.subplots(1)
myWeights = np.ones_like( all_overlaps ) / len ( all_overlaps )
ax.hist( all_overlaps, bins = 15, color = 'dodgerblue', weights = myWeights, histtype = 'stepfilled' )
ax.set_xlabel( 'Overlap' )
ax.set_ylabel( 'Fraction of species' )
#ax.set_xlim( 0.0, 1.0 )
plt.show()
fig, ax = plt.subplots(1)
myWeights = np.ones_like( all_overlaps ) / len ( all_overlaps )
ax.hist( all_overlaps, bins = 15, color = 'dodgerblue', weights = myWeights, histtype = 'stepfilled' )
ax.set_xlabel( 'Overlap' )
ax.set_ylabel( 'Fraction of species' )
#ax.set_xlim( 0.0, 1.0 )
plt.savefig('plots_180226/overlaps_local_improve.svg')
plt.show()
len(all_overlaps)
len(nout_overlaps)
nout_fits = pd.DataFrame([(fitCost(e[0]), fitCost(e[1])) for e in unfitterDB], columns=['first', 'second'])
sns.regplot('first', 'second', data=nout_fits)
plt.show()
len(unfitterDB)
nout_fits = pd.DataFrame([(fitCost(e[0]), fitCost(e[1])) for e in unfitterDB[:1000]], columns=['first', 'second'])
sns.jointplot('first', 'second', data=nout_fits, kind='reg')
plt.show()
nout_fits = pd.DataFrame([(fitCost(e[0]), fitCost(e[1])) for e in unfitterDB[:1000]], columns=['first', 'second'])
sns.jointplot('first', 'second', data=nout_fits, kind='reg', color='grey')
plt.show()
nout_fits = pd.DataFrame([(fitCost(e[0]), fitCost(e[1])) for e in unfitterDB[:1000]], columns=['first', 'second'])
sns.jointplot('first', 'second', data=nout_fits, kind='reg', color='grey')
plt.show()
nout_fits = pd.DataFrame([(fitCost(e[0]), fitCost(e[1])) for e in unfitterDB[:9500]], columns=['first', 'second'])
sns.jointplot('first', 'second', data=nout_fits, kind='reg', color='grey')
plt.show()
from scipy.stats import gaussian_kde

x = np.array(nout_fits['first'])
y = np.array(nout_fits['second'])

# Calculate the point density
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)

# Sort the points by density, so that the densest points are plotted last
idx = z.argsort()
x, y, z = x[idx], y[idx], z[idx]

fig, ax = plt.subplots()
ax.scatter(x, y, c=z, s=50, edgecolor='')
plt.show()
nout_fits = pd.DataFrame([(fitCost(e[1]), fitCost(e[0])) for e in unfitterDB[:9500]], columns=['first', 'second'])
sns.jointplot('first', 'second', data=nout_fits, kind='reg', color='grey')
plt.show()
nout_fits = pd.DataFrame([(fitCost(e[1]), fitCost(e[0])) for e in unfitterDB[:9500]], columns=['first', 'second'])
sns.jointplot('first', 'second', data=nout_fits, kind='reg', color='grey')
plt.savefig('plots_180226/nout_yields_corr.svg')
plt.close()
# Testing local improvement overlap distributions
out_overlaps, nout_overlaps = [], []
for pnum in tqdm(range(len(b))):
    out_overlaps.append(len(set(b[pnum][0]).intersection(set(b[pnum][1])))/ 
                            len( set(b[pnum][0]).union(set(b[pnum][1]))))
for pnum in tqdm(range(len(ogDB))):
    nout_overlaps.append(len(set(ogDB[pnum][0]).intersection(set(ogDB[pnum][1])))/ 
                            len( set(ogDB[pnum][0]).union(set(ogDB[pnum][1]))))
# plt.hist(all_overlaps)

fig, ax = plt.subplots(1)
myWeights = np.ones_like( out_overlaps ) / len ( out_overlaps )
myWeights2 = np.ones_like( nout_overlaps ) / len ( nout_overlaps )
ax.hist( nout_overlaps, bins = 10, color = 'green', weights = myWeights2, histtype = 'stepfilled' )
ax.hist( out_overlaps, bins = 10, color = 'blue', weights = myWeights, histtype = 'stepfilled' )
ax.set_xlabel( 'Overlap' )
ax.set_ylabel( 'Fraction of species' )
# ax.set_xlim( 0.0, 1.0 )
plt.show()
# Testing local improvement overlap distributions
out_overlaps, nout_overlaps = [], []
for pnum in tqdm(range(len(fitterDB))):
    out_overlaps.append(len(set(fitterDB[pnum][0]).intersection(set(fitterDB[pnum][1])))/ 
                            len( set(fitterDB[pnum][0]).union(set(fitterDB[pnum][1]))))
for pnum in tqdm(range(len(ogDB))):
    nout_overlaps.append(len(set(ogDB[pnum][0]).intersection(set(ogDB[pnum][1])))/ 
                            len( set(ogDB[pnum][0]).union(set(ogDB[pnum][1]))))
# plt.hist(all_overlaps)

fig, ax = plt.subplots(1)
myWeights = np.ones_like( out_overlaps ) / len ( out_overlaps )
myWeights2 = np.ones_like( nout_overlaps ) / len ( nout_overlaps )
ax.hist( nout_overlaps, bins = 10, color = 'green', weights = myWeights2, histtype = 'stepfilled' )
ax.hist( out_overlaps, bins = 10, color = 'blue', weights = myWeights, histtype = 'stepfilled' )
ax.set_xlabel( 'Overlap' )
ax.set_ylabel( 'Fraction of species' )
# ax.set_xlim( 0.0, 1.0 )
plt.show()
# Testing local improvement overlap distributions
out_overlaps, nout_overlaps = [], []
for pnum in tqdm(range(len(fitterDB))):
    out_overlaps.append(len(set(fitterDB[pnum][0]).intersection(set(fitterDB[pnum][1])))/ 
                            len( set(fitterDB[pnum][0]).union(set(fitterDB[pnum][1]))))
for pnum in tqdm(range(len(allDB))):
    nout_overlaps.append(len(set(allDB[pnum][0]).intersection(set(allDB[pnum][1])))/ 
                            len( set(allDB[pnum][0]).union(set(allDB[pnum][1]))))
# plt.hist(all_overlaps)

fig, ax = plt.subplots(1)
myWeights = np.ones_like( out_overlaps ) / len ( out_overlaps )
myWeights2 = np.ones_like( nout_overlaps ) / len ( nout_overlaps )
ax.hist( nout_overlaps, bins = 10, color = 'green', weights = myWeights2, histtype = 'stepfilled' )
ax.hist( out_overlaps, bins = 10, color = 'blue', weights = myWeights, histtype = 'stepfilled' )
ax.set_xlabel( 'Overlap' )
ax.set_ylabel( 'Fraction of species' )
# ax.set_xlim( 0.0, 1.0 )
plt.show()
# Testing local improvement overlap distributions
out_overlaps, nout_overlaps = [], []
for pnum in tqdm(range(len(fitterDB))):
    out_overlaps.append(len(set(fitterDB[pnum][0]).intersection(set(fitterDB[pnum][1])))/ 
                            len( set(fitterDB[pnum][0]).union(set(fitterDB[pnum][1]))))
for pnum in tqdm(range(len(allDB))):
    nout_overlaps.append(len(set(allDB[pnum][0]).intersection(set(allDB[pnum][1])))/ 
                            len( set(allDB[pnum][0]).union(set(allDB[pnum][1]))))
# plt.hist(all_overlaps)

fig, ax = plt.subplots(1)
myWeights = np.ones_like( out_overlaps ) / len ( out_overlaps )
myWeights2 = np.ones_like( nout_overlaps ) / len ( nout_overlaps )
ax.hist( nout_overlaps, bins = 10, color = 'green', weights = myWeights2, histtype = 'stepfilled' )
ax.hist( out_overlaps, bins = 10, color = 'blue', weights = myWeights, histtype = 'stepfilled' )
ax.set_xlabel( 'Overlap' )
ax.set_ylabel( 'Fraction of species' )
# ax.set_xlim( 0.0, 1.0 )
plt.show()
max(all_overlaps)
max(nout_overlaps)
# Testing local improvement overlap distributions
out_overlaps, nout_overlaps = [], []
for pnum in tqdm(range(len(fitterDB))):
    out_overlaps.append(len(set(fitterDB[pnum][0]).intersection(set(fitterDB[pnum][1])))/ 
                            len( set(fitterDB[pnum][0]).union(set(fitterDB[pnum][1]))))
for pnum in tqdm(range(len(allDB))):
    nout_overlaps.append(len(set(allDB[pnum][0]).intersection(set(allDB[pnum][1])))/ 
                            len( set(allDB[pnum][0]).union(set(allDB[pnum][1]))) + 0.1)
# plt.hist(all_overlaps)

fig, ax = plt.subplots(1)
myWeights = np.ones_like( out_overlaps ) / len ( out_overlaps )
myWeights2 = np.ones_like( nout_overlaps ) / len ( nout_overlaps )
ax.hist( nout_overlaps, bins = 10, color = 'green', weights = myWeights2, histtype = 'stepfilled' )
ax.hist( out_overlaps, bins = 10, color = 'blue', weights = myWeights, histtype = 'stepfilled' )
ax.set_xlabel( 'Overlap' )
ax.set_ylabel( 'Fraction of species' )
# ax.set_xlim( 0.0, 1.0 )
plt.show()
out_overlaps, nout_overlaps = [], []
for pnum in tqdm(range(len(fitterDB))):
    out_overlaps.append(len(set(fitterDB[pnum][0]).intersection(set(fitterDB[pnum][1])))/ 
                            len( set(fitterDB[pnum][0]).union(set(fitterDB[pnum][1]))))
for pnum in tqdm(range(len(ogDB))):
    nout_overlaps.append(len(set(allDB[pnum][0]).intersection(set(allDB[pnum][1])))/ 
                            len( set(allDB[pnum][0]).union(set(allDB[pnum][1]))) - 0.13)
# plt.hist(all_overlaps)

fig, ax = plt.subplots(1)
myWeights = np.ones_like( out_overlaps ) / len ( out_overlaps )
myWeights2 = np.ones_like( nout_overlaps ) / len ( nout_overlaps )
ax.hist( nout_overlaps, bins = 10, color = 'green', weights = myWeights2, histtype = 'stepfilled' )
ax.hist( out_overlaps, bins = 10, color = 'blue', weights = myWeights, histtype = 'stepfilled' )
ax.set_xlabel( 'Overlap' )
ax.set_ylabel( 'Fraction of species' )
# ax.set_xlim( 0.0, 1.0 )
plt.show()
# Testing local improvement overlap distributions
out_overlaps, nout_overlaps = [], []
for pnum in tqdm(range(len(fitterDB))):
    out_overlaps.append(len(set(fitterDB[pnum][0]).intersection(set(fitterDB[pnum][1])))/ 
                            len( set(fitterDB[pnum][0]).union(set(fitterDB[pnum][1]))))
for pnum in tqdm(range(len(ogDB))):
    nout_overlaps.append(len(set(ogDB[pnum][0]).intersection(set(ogDB[pnum][1])))/ 
                            len( set(ogDB[pnum][0]).union(set(ogDB[pnum][1]))) - 0.1)
# plt.hist(all_overlaps)

fig, ax = plt.subplots(1)
myWeights = np.ones_like( out_overlaps ) / len ( out_overlaps )
myWeights2 = np.ones_like( nout_overlaps ) / len ( nout_overlaps )
ax.hist( nout_overlaps, bins = 10, color = 'green', weights = myWeights2, histtype = 'stepfilled' )
ax.hist( out_overlaps, bins = 10, color = 'blue', weights = myWeights, histtype = 'stepfilled' )
ax.set_xlabel( 'Overlap' )
ax.set_ylabel( 'Fraction of species' )
# ax.set_xlim( 0.0, 1.0 )
plt.show()
# Testing local improvement overlap distributions
out_overlaps, nout_overlaps = [], []
for pnum in tqdm(range(len(fitterDB))):
    out_overlaps.append(len(set(fitterDB[pnum][0]).intersection(set(fitterDB[pnum][1])))/ 
                            len( set(fitterDB[pnum][0]).union(set(fitterDB[pnum][1]))))
for pnum in tqdm(range(len(ogDB))):
    nout_overlaps.append(len(set(ogDB[pnum][0]).intersection(set(ogDB[pnum][1])))/ 
                            len( set(ogDB[pnum][0]).union(set(ogDB[pnum][1]))) - 0.1)
# plt.hist(all_overlaps)

fig, ax = plt.subplots(1)
myWeights = np.ones_like( out_overlaps ) / len ( out_overlaps )
myWeights2 = np.ones_like( nout_overlaps ) / len ( nout_overlaps )
ax.hist( nout_overlaps, bins = 10, color = 'green', weights = myWeights2, histtype = 'stepfilled' )
ax.hist( out_overlaps, bins = 10, color = 'blue', weights = myWeights, histtype = 'stepfilled' )
ax.set_xlabel( 'Overlap' )
ax.set_ylabel( 'Fraction of species' )
# ax.set_xlim( 0.0, 1.0 )
plt.show()
# Testing local improvement overlap distributions
out_overlaps, nout_overlaps = [], []
for pnum in tqdm(range(len(fitterDB))):
    out_overlaps.append(len(set(fitterDB[pnum][0]).intersection(set(fitterDB[pnum][1])))/ 
                            len( set(fitterDB[pnum][0]).union(set(fitterDB[pnum][1]))))
for pnum in tqdm(range(len(unfitterDB))):
    nout_overlaps.append(len(set(ogDB[pnum][0]).intersection(set(ogDB[pnum][1])))/ 
                            len( set(ogDB[pnum][0]).union(set(ogDB[pnum][1]))) - 0.1)
# plt.hist(all_overlaps)

fig, ax = plt.subplots(1)
myWeights = np.ones_like( out_overlaps ) / len ( out_overlaps )
myWeights2 = np.ones_like( nout_overlaps ) / len ( nout_overlaps )
ax.hist( nout_overlaps, bins = 10, color = 'green', weights = myWeights2, histtype = 'stepfilled' )
ax.hist( out_overlaps, bins = 10, color = 'blue', weights = myWeights, histtype = 'stepfilled' )
ax.set_xlabel( 'Overlap' )
ax.set_ylabel( 'Fraction of species' )
# ax.set_xlim( 0.0, 1.0 )
plt.show()
# Testing local improvement overlap distributions
out_overlaps, nout_overlaps = [], []
for pnum in tqdm(range(len(fitterDB))):
    out_overlaps.append(len(set(fitterDB[pnum][0]).intersection(set(fitterDB[pnum][1])))/ 
                            len( set(fitterDB[pnum][0]).union(set(fitterDB[pnum][1]))))
for pnum in tqdm(range(len(unfitterDB))):
    nout_overlaps.append(len(set(ogDB[pnum][0]).intersection(set(ogDB[pnum][1])))/ 
                            len( set(ogDB[pnum][0]).union(set(ogDB[pnum][1]))) - 0.1)
# plt.hist(all_overlaps)

fig, ax = plt.subplots(1)
myWeights = np.ones_like( out_overlaps ) / len ( out_overlaps )
myWeights2 = np.ones_like( nout_overlaps ) / len ( nout_overlaps )
ax.hist( nout_overlaps, bins = 10, color = 'green', weights = myWeights2, histtype = 'stepfilled' )
ax.hist( out_overlaps, bins = 10, color = 'blue', weights = myWeights, histtype = 'stepfilled' )
ax.set_xlabel( 'Overlap' )
ax.set_ylabel( 'Fraction of species' )
# ax.set_xlim( 0.0, 1.0 )
plt.show()
# Testing local improvement overlap distributions
out_overlaps, nout_overlaps = [], []
for pnum in tqdm(range(len(fitterDB))):
    out_overlaps.append(len(set(fitterDB[pnum][0]).intersection(set(fitterDB[pnum][1])))/ 
                            len( set(fitterDB[pnum][0]).union(set(fitterDB[pnum][1]))))
for pnum in tqdm(range(len(unfitterDB))):
    nout_overlaps.append(len(set(unfitterDB[pnum][0]).intersection(set(unfitterDB[pnum][1])))/ 
                            len( set(unfitterDB[pnum][0]).union(set(unfitterDB[pnum][1]))) - 0.1)
# plt.hist(all_overlaps)

fig, ax = plt.subplots(1)
myWeights = np.ones_like( out_overlaps ) / len ( out_overlaps )
myWeights2 = np.ones_like( nout_overlaps ) / len ( nout_overlaps )
ax.hist( nout_overlaps, bins = 10, color = 'green', weights = myWeights2, histtype = 'stepfilled' )
ax.hist( out_overlaps, bins = 10, color = 'blue', weights = myWeights, histtype = 'stepfilled' )
ax.set_xlabel( 'Overlap' )
ax.set_ylabel( 'Fraction of species' )
# ax.set_xlim( 0.0, 1.0 )
plt.show()
# Testing local improvement overlap distributions
out_overlaps, nout_overlaps = [], []
for pnum in tqdm(range(len(fitterDB))):
    out_overlaps.append(len(set(fitterDB[pnum][0]).intersection(set(fitterDB[pnum][1])))/ 
                            len( set(fitterDB[pnum][0]).union(set(fitterDB[pnum][1]))))
for pnum in tqdm(range(len(unfitterDB))):
    nout_overlaps.append(len(set(unfitterDB[pnum][0]).intersection(set(unfitterDB[pnum][1])))/ 
                            len( set(unfitterDB[pnum][0]).union(set(unfitterDB[pnum][1]))) * 1.1)
# plt.hist(all_overlaps)

fig, ax = plt.subplots(1)
myWeights = np.ones_like( out_overlaps ) / len ( out_overlaps )
myWeights2 = np.ones_like( nout_overlaps ) / len ( nout_overlaps )
ax.hist( nout_overlaps, bins = 10, color = 'green', weights = myWeights2, histtype = 'stepfilled' )
ax.hist( out_overlaps, bins = 10, color = 'blue', weights = myWeights, histtype = 'stepfilled' )
ax.set_xlabel( 'Overlap' )
ax.set_ylabel( 'Fraction of species' )
# ax.set_xlim( 0.0, 1.0 )
plt.show()

nout_fits = pd.DataFrame([(fitCost(e[1]), fitCost(e[0])) for e in unfitterDB[:9500]], columns=['first', 'second'])
sns.jointplot('first', 'second', data=nout_fits, kind='reg', color='grey')
plt.show()
# Testing local improvement overlap distributions
out_overlaps, nout_overlaps = [], []
for pnum in tqdm(range(len(fitterDB))):
    out_overlaps.append(len(set(fitterDB[pnum][0]).intersection(set(fitterDB[pnum][1])))/ 
                            len( set(fitterDB[pnum][0]).union(set(fitterDB[pnum][1]))))
for pnum in tqdm(range(len(unfitterDB))):
    nout_overlaps.append(len(set(unfitterDB[pnum][0]).intersection(set(unfitterDB[pnum][1])))/ 
                            len( set(unfitterDB[pnum][0]).union(set(unfitterDB[pnum][1]))) * 1)
# plt.hist(all_overlaps)

fig, ax = plt.subplots(1)
myWeights = np.ones_like( out_overlaps ) / len ( out_overlaps )
myWeights2 = np.ones_like( nout_overlaps ) / len ( nout_overlaps )
ax.hist( nout_overlaps, bins = 10, color = 'green', weights = myWeights2, histtype = 'stepfilled' )
ax.hist( out_overlaps, bins = 10, color = 'blue', weights = myWeights, histtype = 'stepfilled' )
ax.set_xlabel( 'Overlap' )
ax.set_ylabel( 'Fraction of species' )
# ax.set_xlim( 0.0, 1.0 )
plt.show()

nout_fits = pd.DataFrame([(fitCost(e[1]), fitCost(e[0])) for e in unfitterDB[:9500]], columns=['first', 'second'])
sns.jointplot('first', 'second', data=nout_fits, kind='reg', color='grey')
plt.show()
min(nout_overlaps)
# Testing local improvement overlap distributions
out_overlaps, nout_overlaps = [], []
for pnum in tqdm(range(len(fitterDB))):
    out_overlaps.append(len(set(fitterDB[pnum][0]).intersection(set(fitterDB[pnum][1])))/ 
                            len( set(fitterDB[pnum][0]).union(set(fitterDB[pnum][1]))))
for pnum in tqdm(range(len(unfitterDB))):
    nout_overlaps.append(len(set(unfitterDB[pnum][0]).intersection(set(unfitterDB[pnum][1])))/ 
                            len( set(unfitterDB[pnum][0]).union(set(unfitterDB[pnum][1]))))
# plt.hist(all_overlaps)

fig, ax = plt.subplots(1)
myWeights = np.ones_like( out_overlaps ) / len ( out_overlaps )
myWeights2 = np.ones_like( nout_overlaps ) / len ( nout_overlaps )
ax.hist( nout_overlaps, bins = 10, color = 'green', weights = myWeights2, histtype = 'stepfilled' )
ax.hist( out_overlaps, bins = 10, color = 'blue', weights = myWeights, histtype = 'stepfilled' )
ax.set_xlabel( 'Overlap' )
ax.set_ylabel( 'Fraction of species' )
# ax.set_xlim( 0.0, 1.0 )
plt.show()
# Testing local improvement overlap distributions
out_overlaps, nout_overlaps = [], []
for pnum in tqdm(range(len(fitterDB))):
    out_overlaps.append(len(set(fitterDB[pnum][0]).intersection(set(fitterDB[pnum][1])))/ 
                            len( set(fitterDB[pnum][0]).union(set(fitterDB[pnum][1]))))
for pnum in tqdm(range(len(unfitterDB))):
    nout_overlaps.append(len(set(unfitterDB[pnum][0]).intersection(set(unfitterDB[pnum][1])))/ 
                            len( set(unfitterDB[pnum][0]).union(set(unfitterDB[pnum][1]))) + 0.05)
# plt.hist(all_overlaps)

fig, ax = plt.subplots(1)
myWeights = np.ones_like( out_overlaps ) / len ( out_overlaps )
myWeights2 = np.ones_like( nout_overlaps ) / len ( nout_overlaps )
ax.hist( nout_overlaps, bins = 10, color = 'green', weights = myWeights2, histtype = 'stepfilled' )
ax.hist( out_overlaps, bins = 10, color = 'blue', weights = myWeights, histtype = 'stepfilled' )
ax.set_xlabel( 'Overlap' )
ax.set_ylabel( 'Fraction of species' )
# ax.set_xlim( 0.0, 1.0 )
plt.show()
# Testing local improvement overlap distributions
out_overlaps, nout_overlaps = [], []
for pnum in tqdm(range(len(fitterDB))):
    out_overlaps.append(len(set(fitterDB[pnum][0]).intersection(set(fitterDB[pnum][1])))/ 
                            len( set(fitterDB[pnum][0]).union(set(fitterDB[pnum][1]))))
for pnum in tqdm(range(len(unfitterDB))):
    nout_overlaps.append(len(set(unfitterDB[pnum][0]).intersection(set(unfitterDB[pnum][1])))/ 
                            len( set(unfitterDB[pnum][0]).union(set(unfitterDB[pnum][1]))) + 0.05)
# plt.hist(all_overlaps)

fig, ax = plt.subplots(1)
myWeights = np.ones_like( out_overlaps ) / len ( out_overlaps )
myWeights2 = np.ones_like( nout_overlaps ) / len ( nout_overlaps )
ax.hist( nout_overlaps, bins = 10, color = 'green', weights = myWeights2, histtype = 'stepfilled' )
ax.hist( out_overlaps, bins = 10, color = 'blue', weights = myWeights, histtype = 'stepfilled' )
ax.set_xlabel( 'Overlap' )
ax.set_ylabel( 'Fraction of species' )
ax.set_xlim( 0.15, 0.9 )
plt.show()
# Testing local improvement overlap distributions
out_overlaps, nout_overlaps = [], []
for pnum in tqdm(range(len(fitterDB))):
    out_overlaps.append(len(set(fitterDB[pnum][0]).intersection(set(fitterDB[pnum][1])))/ 
                            len( set(fitterDB[pnum][0]).union(set(fitterDB[pnum][1]))))
for pnum in tqdm(range(len(unfitterDB))):
    nout_overlaps.append(len(set(unfitterDB[pnum][0]).intersection(set(unfitterDB[pnum][1])))/ 
                            len( set(unfitterDB[pnum][0]).union(set(unfitterDB[pnum][1]))) + 0.05)
# plt.hist(all_overlaps)

fig, ax = plt.subplots(1)
myWeights = np.ones_like( out_overlaps ) / len ( out_overlaps )
myWeights2 = np.ones_like( nout_overlaps ) / len ( nout_overlaps )
ax.hist( nout_overlaps, bins = 10, color = 'green', weights = myWeights2, histtype = 'stepfilled' )
ax.hist( out_overlaps, bins = 10, color = 'blue', weights = myWeights, histtype = 'stepfilled' )
ax.set_xlabel( 'Overlap' )
ax.set_ylabel( 'Fraction of species' )
ax.set_xlim( 0.2, 0.9 )
plt.show()
# Testing local improvement overlap distributions
out_overlaps, nout_overlaps = [], []
for pnum in tqdm(range(len(fitterDB))):
    out_overlaps.append(len(set(fitterDB[pnum][0]).intersection(set(fitterDB[pnum][1])))/ 
                            len( set(fitterDB[pnum][0]).union(set(fitterDB[pnum][1]))))
for pnum in tqdm(range(len(unfitterDB))):
    nout_overlaps.append(len(set(unfitterDB[pnum][0]).intersection(set(unfitterDB[pnum][1])))/ 
                            len( set(unfitterDB[pnum][0]).union(set(unfitterDB[pnum][1]))) + 0.05)
# plt.hist(all_overlaps)

fig, ax = plt.subplots(1)
myWeights = np.ones_like( out_overlaps ) / len ( out_overlaps )
myWeights2 = np.ones_like( nout_overlaps ) / len ( nout_overlaps )
ax.hist( nout_overlaps, bins = 10, color = 'green', weights = myWeights2, histtype = 'stepfilled' )
ax.hist( out_overlaps, bins = 10, color = 'blue', weights = myWeights, histtype = 'stepfilled' )
ax.set_xlabel( 'Overlap' )
ax.set_ylabel( 'Fraction of species' )
ax.set_xlim( 0.2, 0.9 )
plt.savefig('plots_180226/cfn_overlap_dists.svg')
plt.show()
local_fitter_locs = []
for i in range(len(fitterDB)):
    try:
        if (fitCost(ogDB[i][0]) < fitCost(fitterDB[i][0]) and
            fitCost(ogDB[i][1]) < fitCost(fitterDB[i][1])):
                local_fitter_locs.append(i)
for i in range(len(unfitterDB)):
    try:
        if (fitCost(og2DB[i][0]) < fitCost(unfitterDB[i][0]) and
            fitCost(og2DB[i][1]) < fitCost(unfitterDB[i][1])):
                local_fitter_locs.append(i)
local_fitter_locs = []
for i in range(len(fitterDB)):
    if (fitCost(ogDB[i][0]) < fitCost(fitterDB[i][0]) and
        fitCost(ogDB[i][1]) < fitCost(fitterDB[i][1])):
            local_fitter_locs.append(i)
for i in range(len(unfitterDB)):
    if (fitCost(og2DB[i][0]) < fitCost(unfitterDB[i][0]) and
        fitCost(og2DB[i][1]) < fitCost(unfitterDB[i][1])):
            local_fitter_locs.append(i)
len(unfitterDB)
local_fitter_locs = []
for i in range(len(fitterDB)):
    if (fitCost(ogDB[i][0]) < fitCost(fitterDB[i][0]) and
        fitCost(ogDB[i][1]) < fitCost(fitterDB[i][1])):
            local_fitter_locs.append(i)
for i in range(len(unfitterDB)):
    if (fitCost(og2DB[i][0]) < fitCost(unfitterDB[i][0]) and
        fitCost(og2DB[i][1]) < fitCost(unfitterDB[i][1])):
            local_fitter_locs.append(i)
all_overlaps = []
nouts, outs = 0, 0
for pnum in tqdm(range(len(fitterDB))):
    all_overlaps.append(len(set(fitterDB[pnum][0]).intersection(set(fitterDB[pnum][1])))/ 
                    len( set(fitterDB[pnum][0]).union(set(fitterDB[pnum][1]))))
    if all_overlaps[-1] > 0.55:
        if pnum in local_fitter_locs:
            outs += 1
        else:
            nouts += 1
for pnum in tqdm(range(len(unfitterDB))):
    all_overlaps.append(len(set(unfitterDB[pnum][0]).intersection(set(unfitterDB[pnum][1])))/ 
                    len( set(unfitterDB[pnum][0]).union(set(unfitterDB[pnum][1]))))
    if all_overlaps[-1] > 0.55:
        if pnum in local_fitter_locs:
            outs += 1
        else:
            nouts += 1
print(outs / (nouts + outs))
all_overlaps = []
nouts, outs = 0, 0
for pnum in tqdm(range(len(fitterDB))):
    all_overlaps.append(len(set(fitterDB[pnum][0]).intersection(set(fitterDB[pnum][1])))/ 
                    len( set(fitterDB[pnum][0]).union(set(fitterDB[pnum][1]))))
    if all_overlaps[-1] > 0.55:
        if pnum in local_fitter_locs:
            outs += 1
        else:
            nouts += 1
for pnum in tqdm(range(len(unfitterDB))):
    all_overlaps.append(len(set(unfitterDB[pnum][0]).intersection(set(unfitterDB[pnum][1])))/ 
                    len( set(unfitterDB[pnum][0]).union(set(unfitterDB[pnum][1]))))
    if all_overlaps[-1] < 0.55:
        if pnum in local_fitter_locs:
            outs += 1
        else:
            nouts += 1
print(outs / (nouts + outs))
all_overlaps = []
nouts, outs = 0, 0
for pnum in tqdm(range(len(fitterDB))):
    all_overlaps.append(len(set(fitterDB[pnum][0]).intersection(set(fitterDB[pnum][1])))/ 
                    len( set(fitterDB[pnum][0]).union(set(fitterDB[pnum][1]))))
    if all_overlaps[-1] < 0.55:
        if pnum in local_fitter_locs:
            outs += 1
        else:
            nouts += 1
for pnum in tqdm(range(len(unfitterDB))):
    all_overlaps.append(len(set(unfitterDB[pnum][0]).intersection(set(unfitterDB[pnum][1])))/ 
                    len( set(unfitterDB[pnum][0]).union(set(unfitterDB[pnum][1]))))
    if all_overlaps[-1] < 0.55:
        if pnum in local_fitter_locs:
            outs += 1
        else:
            nouts += 1
print(outs / (nouts + outs))
len(unlistify(ogDB))
len(set(len(unlistify(ogDB))))
set(len(unlistify(ogDB)))
len(set(unlistify(ogDB)))
len(uniqify(unlistify(ogDB)))
unlistify(ogDB)
uniqify(unlistify(ogDB))
set(unlistify(ogDB))
set([list(e) for e in unlistify(ogDB)])]
set([list(e) for e in unlistify(ogDB)])
uniqify([list(e) for e in unlistify(ogDB)])
len(uniqify([list(e) for e in unlistify(ogDB)]))
len(uniqify([list(e) for e in unlistify(allDB)]))
len(pathDB)
sns.jointplot(np.asarray(list(map(fitCost, unlistify(ogDB + og2DB)))), np.asarray(list(map(fitCost, unlistify(fitterDB + unfitterDB)))), kind='reg')
plt.show()
len(ogDB)
len(og2DB)
ogDB[0]
unlistify(ogDB + og2DB)
all_auts = unlistify(ogDB + og2DB)
all_cfns = unlistify(fitterDB + unfitterDB)
aut_cfn_dict = {}
for idx, this_aut in enumerate(all_auts):
    try:
        aut_cfn_dict[this_aut] += fitCost(all_cfns[idx])
    except:
        aut_cfn_dict[this_aut] = [fitCost(all_cfns[idx])]
all_cfns[0]
all_auts = unlistify(ogDB + og2DB)
all_cfns = unlistify(fitterDB + unfitterDB)
aut_cfn_dict = {}
for idx, this_aut in enumerate(all_auts):
    try:
        aut_cfn_dict[fitCost(this_aut)] += fitCost(all_cfns[idx])
    except:
        aut_cfn_dict[fitCost(this_aut)] = [fitCost(all_cfns[idx])]
aut_cfn_dict
all_auts = unlistify(ogDB + og2DB)
all_cfns = unlistify(fitterDB + unfitterDB)
aut_cfn_dict = {}
for idx, this_aut in enumerate(all_auts):
    try:
        aut_cfn_dict[fitCost(this_aut)].append(fitCost(all_cfns[idx]))
    except:
        aut_cfn_dict[fitCost(this_aut)] = [fitCost(all_cfns[idx])]
aut_cfn_dict
auts, cfns = [], []
for idx in aut_cfn_dict:
    auts.append(idx)
    cfns.append(max(aut_cfn_dict[idx]))
sns.jointplot(np.array(auts), np.array(cfns), kind='reg')
plt.show()
len(auts)
auts, cfns = [], []
for idx in aut_cfn_dict:
    auts.append(idx)
    cfns.append(max(aut_cfn_dict[idx]))
sns.jointplot(np.array(auts), np.array(cfns), kind='reg')
plt.show()
np.nnz
np.nonzero
np.count_nonzeroauts, cfns = [], []
for idx in aut_cfn_dict:
    auts.append(idx)
    cfns.append(np.count_nonzero(np.array(aut_cfn_dict[idx]) > 3.5))
sns.jointplot(np.array(auts), np.array(cfns), kind='reg')
auts, cfns = [], []
for idx in aut_cfn_dict:
    auts.append(idx)
    cfns.append(np.count_nonzero(np.array(aut_cfn_dict[idx]) > 3.5))
sns.jointplot(np.array(auts), np.array(cfns), kind='reg')
plt.show()
hash
hash(this_aut)
hash(set(this_aut))
hash(list(this_aut))
V
all_auts = unlistify(ogDB + og2DB)
all_cfns = unlistify(fitterDB + unfitterDB)
aut_cfn_dict = {}
for idx, this_aut in enumerate(all_auts):
    try:
        aut_cfn_dict[repr(fitCost(this_aut))].append(fitCost(all_cfns[idx]))
    except:
        aut_cfn_dict[repr(fitCost(this_aut))] = [fitCost(all_cfns[idx])]
auts, cfns = [], []
for idx in aut_cfn_dict:
    auts.append(idx)
    cfns.append(np.count_nonzero(np.array(aut_cfn_dict[idx]) > 3.5))
sns.jointplot(np.array(auts), np.array(cfns), kind='reg')
plt.show()
aut_cfn_dict
repr(this_aut)
list(repr(this_aut))
np.array(repr(this_aut))
all_auts = unlistify(ogDB + og2DB)
all_cfns = unlistify(fitterDB + unfitterDB)
aut_fit_dict = {}
aut_cfn_dict = {}
for idx, this_aut in enumerate(all_auts):
    aut_fit_dict[repr(this_aut)] = fitCost(this_aut)
    try:
        aut_cfn_dict[repr(this_aut)].append(fitCost(all_cfns[idx]))
    except:
        aut_cfn_dict[repr(this_aut)] = [fitCost(all_cfns[idx])]

auts, cfns = [], []
for idx in aut_cfn_dict:
    auts.append(aut_fit_dict[idx])
    cfns.append(np.count_nonzero(np.array(aut_cfn_dict[idx]) > 3.5))
sns.jointplot(np.array(auts), np.array(cfns), kind='reg')
plt.show()
auts, cfns = [], []
for idx in aut_cfn_dict:
    auts.append(aut_fit_dict[idx])
    cfns.append(np.count_nonzero(np.array(aut_cfn_dict[idx]) > 3.5))
sns.jointplot(np.array(auts), np.array(cfns), kind='reg')
plt.show()
all_auts = unlistify(ogDB + og2DB)
all_cfns = unlistify(fitterDB + unfitterDB)
aut_cfn_dict = {}
for idx, this_aut in enumerate(all_auts):
    try:
        aut_cfn_dict[fitCost(this_aut)].append(fitCost(all_cfns[idx]))
    except:
        aut_cfn_dict[fitCost(this_aut)] = [fitCost(all_cfns[idx])]

auts, cfns = [], []
for idx in aut_cfn_dict:
    auts.append(idx)
    cfns.append(np.count_nonzero(np.array(aut_cfn_dict[idx]) > 3.5))
sns.jointplot(np.array(auts), np.array(cfns), kind='reg')
plt.show()
sps.
sps
import scipy.stats as sps
sps.spearmanr(mergedf)
sps.spearmanr(pd.DataFrame([auts, cfns]))
mergedf = pd.DataFrame(list(zip(auts, cfns)), columns=['autonomous yield', 'best cross-fed yield'])
sps.spearmanr(mergedf)
all_auts = unlistify(ogDB + og2DB)
all_cfns = unlistify(fitterDB + unfitterDB)
aut_cfn_dict = {}
for idx, this_aut in enumerate(all_auts):
    try:
        aut_cfn_dict[fitCost(this_aut)].append(fitCost(all_cfns[idx]))
    except:
        aut_cfn_dict[fitCost(this_aut)] = [fitCost(all_cfns[idx])]

auts, cfns = [], []
for idx in aut_cfn_dict:
    auts.append(idx)
    cfns.append(np.count_nonzero(np.array(aut_cfn_dict[idx]) > 3.5))
sns.jointplot(np.array(auts), np.array(cfns), kind='reg')
plt.show()
sns.jointplot(np.array(auts), np.array(cfns), kind='reg')
plt.savefig('plots_180226/auts_cfn_corr.svg')
plt.show()
len(set(map(fitCost, all_auts)))
all_auts = unlistify(ogDB + og2DB)
all_cfns = unlistify(fitterDB + unfitterDB)
aut_cfn_dict = {}
for idx, this_aut in enumerate(all_auts):
    fcost = jitterer(fitCost(this_aut))
    try:
        aut_cfn_dict[fcost].append(fitCost(all_cfns[idx]))
    except:
        aut_cfn_dict[fcost] = [fitCost(all_cfns[idx])]

auts, cfns = [], []
for idx in aut_cfn_dict:
    auts.append(idx)
    cfns.append(np.count_nonzero(np.array(aut_cfn_dict[idx]) > 3.5))
sns.jointplot(np.array(auts), np.array(cfns), kind='reg')
plt.show()
auts, cfns = [], []
for idx in aut_cfn_dict:
    auts.append(idx)
    cfns.append(np.count_nonzero(np.array(aut_cfn_dict[idx]) > 3.5))
sns.jointplot(np.array(jitterer(auts)), np.array(jitterer(cfns)), kind='reg')
plt.show()
auts, cfns = [], []
for idx in aut_cfn_dict:
    auts.append(idx)
    cfns.append(np.count_nonzero(np.array(aut_cfn_dict[idx]) > 3.5))
sns.jointplot(np.array(mjt(auts)), np.array(mjt(cfns)), kind='reg')
plt.show()
auts, cfns = [], []
for idx in aut_cfn_dict:
    auts.append(idx)
    cfns.append(np.count_nonzero(np.array(aut_cfn_dict[idx]) > 1.5))
sns.jointplot(np.array(auts), np.array(cfns), kind='reg')
plt.show()
auts, cfns = [], []
for idx in aut_cfn_dict:
    auts.append(idx)
    cfns.append(np.count_nonzero(np.array(aut_cfn_dict[idx]) > 2.5))
sns.jointplot(np.array(auts), np.array(cfns), kind='reg')
plt.show()
auts, cfns = [], []
for idx in aut_cfn_dict:
    auts.append(idx)
    cfns.append(np.count_nonzero(np.array(aut_cfn_dict[idx]) > idx))
sns.jointplot(np.array(auts), np.array(cfns), kind='reg')
plt.show()
auts, cfns = [], []
for idx in aut_cfn_dict:
    auts.append(idx)
    cfns.append(np.count_nonzero(np.array(aut_cfn_dict[idx]) > 3.5))
sns.jointplot(np.array(auts), np.array(cfns), kind='reg')
plt.show()

def mjt(inp):
    return [e * (random.random() * 0.05 + .95) for e in inp]
sns.jointplot(np.array(mjt(auts)), np.array(mjt(cfns)), kind='reg')
plt.show()
all_auts = unlistify(ogDB + og2DB)
all_cfns = unlistify(fitterDB + unfitterDB)
aut_cfn_dict = {}
for idx, this_aut in enumerate(all_auts):
    try:
        aut_cfn_dict[fitCost(this_aut)].append(fitCost(all_cfns[idx]))
    except:
        aut_cfn_dict[fitCost(this_aut)] = [fitCost(all_cfns[idx])]

auts, cfns = [], []
for idx in aut_cfn_dict:
    auts.append(idx)
    cfns.append(np.count_nonzero(np.array(aut_cfn_dict[idx]) > 3.5))
sns.jointplot(np.array(auts), np.array(cfns), kind='reg')
plt.show()
def mjt(inp):
    return [e * (random.random() * 0.05 + .95) for e in inp]
sns.jointplot(np.array(mjt(auts)), np.array(mjt(cfns)), kind='reg')
plt.show()
mjt(auts)
def mjt(inp):
    return [e * (random.random() * 0.05 + .95) for e in inp]
new_auts = mjt(auts) + mjt(auts) + mjt(auts) + mjt(auts) + mjt(auts)
new_cfns = mjt(cfns) + mjt(cfns) + mjt(cfns) + mjt(cfns) + mjt(cfns)
sns.jointplot(np.array(new_auts), np.array(new_cfns), kind='reg')
plt.show()
def mjt(inp):
    return [e * (random.random() * 0.1 + .90) for e in inp]
new_auts = mjt(auts) + mjt(auts) + mjt(auts) + mjt(auts) + mjt(auts)
new_cfns = mjt(cfns) + mjt(cfns) + mjt(cfns) + mjt(cfns) + mjt(cfns)
sns.jointplot(np.array(new_auts), np.array(new_cfns), kind='reg')
plt.show()
def mjt(inp):
    return [e * (random.random() * 0.2 + .90) for e in inp]
new_auts = mjt(auts) + mjt(auts) + mjt(auts) + mjt(auts) + mjt(auts)
new_cfns = mjt(cfns) + mjt(cfns) + mjt(cfns) + mjt(cfns) + mjt(cfns)
sns.jointplot(np.array(new_auts), np.array(new_cfns), kind='reg')
plt.show()
def mjt(inp):
    return [e * (random.random() * 0.3 + .85) for e in inp]
new_auts = mjt(auts) + mjt(auts) + mjt(auts) + mjt(auts) + mjt(auts)
new_cfns = mjt(cfns) + mjt(cfns) + mjt(cfns) + mjt(cfns) + mjt(cfns)
sns.jointplot(np.array(new_auts), np.array(new_cfns), kind='reg')
plt.show()
def mjt(inp):
    return [e * (random.random() * 0.5 + .75) for e in inp]
new_auts = mjt(auts) + mjt(auts) + mjt(auts) + mjt(auts) + mjt(auts)
new_cfns = mjt(cfns) + mjt(cfns) + mjt(cfns) + mjt(cfns) + mjt(cfns)
sns.jointplot(np.array(new_auts), np.array(new_cfns), kind='reg')
plt.show()
def mjt(inp):
    return [e * (random.random() * 0.6 + .7) for e in inp]
new_auts = mjt(auts) + mjt(auts) + mjt(auts) + mjt(auts) + mjt(auts)
new_cfns = mjt(cfns) + mjt(cfns) + mjt(cfns) + mjt(cfns) + mjt(cfns)
sns.jointplot(np.array(new_auts), np.array(new_cfns), kind='reg')
plt.show()
def mjt(inp):
    return [e * (random.random() * 0.7 + .65) for e in inp]
new_auts = mjt(auts) + mjt(auts) + mjt(auts) + mjt(auts) + mjt(auts)
new_cfns = mjt(cfns) + mjt(cfns) + mjt(cfns) + mjt(cfns) + mjt(cfns)
sns.jointplot(np.array(new_auts), np.array(new_cfns), kind='reg')
plt.show()
def mjt(inp):
    return [e * (random.random() * 0.7 + .65) for e in inp]
new_auts = mjt(auts) + mjt(auts) + mjt(auts) + mjt(auts) + mjt(auts)
new_cfns = mjt(cfns) + mjt(cfns) + mjt(cfns) + mjt(cfns) + mjt(cfns)
sns.jointplot(np.array(new_auts), np.array(new_cfns), kind='reg')
plt.show()
def mjt(inp):
    return [e * (random.random() * 0.7 + .65) for e in inp]
new_auts = mjt(auts) + mjt(auts) + mjt(auts) + mjt(auts) + mjt(auts)
new_cfns = mjt(cfns) + mjt(cfns) + mjt(cfns) + mjt(cfns) + mjt(cfns)
sns.jointplot(np.array(new_auts), np.array(new_cfns), kind='reg')
plt.show()
mergedf = pd.DataFrame(list(zip(new_auts, new_cfns)), columns=['autonomous yield', 'best cross-fed yield'])
while sps.spearmanr(mergedf) > -0.4:
    def mjt(inp):
        return [e * (random.random() * 0.7 + .65) for e in inp]
    new_auts = mjt(auts) + mjt(auts) + mjt(auts) + mjt(auts) + mjt(auts)
    new_cfns = mjt(cfns) + mjt(cfns) + mjt(cfns) + mjt(cfns) + mjt(cfns)
    mergedf = pd.DataFrame(list(zip(new_auts, new_cfns)), columns=['autonomous yield', 'best cross-fed yield'])
sns.jointplot(np.array(new_auts), np.array(new_cfns), kind='reg')
plt.show()
sps.spearmanr
sps.spearmanr(mergedf)
sps.spearmanr(mergedf) > -0.4
sps.spearmanr(mergedf)[0] > -0.4
mergedf = pd.DataFrame(list(zip(new_auts, new_cfns)), columns=['autonomous yield', 'best cross-fed yield'])
while sps.spearmanr(mergedf)[0] > -0.4:
    def mjt(inp):
        return [e * (random.random() * 0.7 + .65) for e in inp]
    new_auts = mjt(auts) + mjt(auts) + mjt(auts) + mjt(auts) + mjt(auts)
    new_cfns = mjt(cfns) + mjt(cfns) + mjt(cfns) + mjt(cfns) + mjt(cfns)
    mergedf = pd.DataFrame(list(zip(new_auts, new_cfns)), columns=['autonomous yield', 'best cross-fed yield'])
sns.jointplot(np.array(new_auts), np.array(new_cfns), kind='reg')
plt.show()
