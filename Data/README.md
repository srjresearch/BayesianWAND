# Data

This directory contains the datasets analysed both within the paper and the supplementary materials.

## Data format

In order for the MCMC code to perform correctly the data must be set up in the correct format. The data files are required to be in space seperated value format (SSV) and contain data over n rows and K columns, where n is the number of observations (rankers) and K is the total number of entities within the data. The observations must be in preference order, that is, column 1 must contain the most preferred entity, column 2 the second most preferred entity and so on. The entities should be labelled numerically from 1 to K.

### Different ranking types

Care must be taken when handling different ranking types. First recall that K<sub>i</sub> denotes the number of entities considered by ranker i and n<sub>i</sub> is the number of positions (ranks) reported by ranker i. Suppose we have K = 6 entities, labelled (1,2,3,4,5,6), then an example of the 4 different ranking types is as follows. 

#### Complete ranking

A complete ranking occurs when a ranker considers and ranks all possible entities and so here n<sub>i</sub> = K<sub>i</sub> = K. An entry of

1 2 3 4 5 6

within row i of the data file would specify that ranker i preferred entity 1 over entity 2, and entity 2 over entity 3 and so on.

#### Partial ranking

Partial rankings occur when a ranker considers a subset of all the entities and reports back a position for each of those considered, and so n<sub>i</sub> = K<sub>i</sub> &lt; K in this scenario.

If ranker i only considered entities (1,2,3) then  n<sub>i</sub> = K<sub>i</sub> = 3 and an entry of

1 2 3 0 0 0

within row i of the data file would specify that ranker i preferred entity 1 over entity 2, and entity 2 over entity 3. The remining entities (4,5,6) were not considered and so do not feature and zeros are used to fill the n<sub>i</sub> - K remaining positions.


3 4 5 6 1 2
4 5 6 1 3 0

