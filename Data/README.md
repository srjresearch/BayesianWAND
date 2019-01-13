# Data

This directory contains the datasets analysed both within the paper and the supplementary materials.

## Data format

In order for the MCMC code to perform correctly the data must be set up in the correct format. The data files are required to be in space seperated value format (SSV) and contain data over n rows and K columns, where n is the number of observations (rankers) and K is the total number of entities within the data. The observations must be in preference order, that is, column 1 must contain the most preferred entity, column 2 the second most preferred entity and so on. The entities should be labelled numerically from 1 to K.

### Different ranking types

Care must be taken when handling different ranking types. First we recall the notation used within the paper and then provde an example of each ranking type and how they should be formatted in the data file.

K<sub>i</sub> is the number of entities considered by ranker i
n<sub>i</sub> is the number of positions (ranks) reported by ranker i

Now suppose we have K = 6 entities, labelled (1,2,3,4,5,6), 


#### Complete ranking

For a complete ranking we have n<sub>i</sub>=K<sub>i</sub>= K 

1 2 3 4 5 6
1 2 3 0 0 0
3 4 5 6 1 2
4 5 6 1 3 0

