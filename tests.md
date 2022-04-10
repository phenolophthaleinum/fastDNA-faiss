# Tests fastDNA-faiss
## make-model.py
||Runtime [s]  |--filter|--dim|--length|--minn   | --maxn| --epoch|
|---|---|---|---|---|---|---|---|
|| 207.352627| family  | 20  | 125  | 10  | 10  | 1  |   
||  214.644864  |  family | 30  | 125  |  10 | 10  | 1  |
||  > 803 (undefined because of long .vec file writing)  |  family | 20  | 125  |  13 | 13  | 1  |
||  < 216 (undefined because of long .vec file writing)   |  family | 12  | 125  |  13 | 13  | 1  |
||  383 (.vec in this)  |  family | 30  | 125  |  10 | 10  | 2  |
|| 524.707605  |  family | 30  | 125  |  10 | 10  | 3 |
||   381.651000 |  family | 60  | 125  |  10 | 10  | 2 |
||   260.463109 |  family | 60  | 200  |  10 | 10  | 2 |
||   362.368333 |  family | 3  | 125  |  10 | 10  | 2 |
||   924.003704 |  genus | 30  | 125  |  10 | 10  | 2 |
||   903.002540 |  genus | 60  | 125  |  6 | 6  | 2 |
||   961.556120 |  genus | 100  | 125  |  6 | 6  | 2 |
||   1740.025691 |  genus (up to 3 representatives) | 60  | 125  |  6 | 6  | 2 |
||   2065.106571 |  none (SlimEdwards)| 60  | 125  |  7 | 7  | 2 |




## workflow.py

|Run id|Runtime [s]  |Run type|--length|--n_samples|
|---|---|---|---|---|
|run-dim_20-len_125-n_100| 115.610985| full  | 125  | 100|
|run-dim_30-len_125-n_100|  139.901045 | full  | 125  |  100| 
|run-dim_20-len_125-n_100-k_13|  [not able to run eventually] 810.893918 (undefined - computer hibernated cause i was watching f1) | full  | 125  |  100|
|run-dim_12-len_125-n_100-k_13| 1800.058877  | full  | 125  |  100|
|run-dim_30-len_125-n_100-epoch_2|  168.952101 | full  | 125  |  100|
|run-dim_30-len_125-n_200-epoch_2| 211.002438  | full  | 125  |  200|
|run-dim_30-len_125-n_100-epoch_3| 155.368846  | full  | 125  |  100|
|run-dim_60-len_125-n_100-epoch_2| 261.372551  | full  | 125  |  100|
|run-dim_60-len_200-n_100-epoch_2| 255.188641  | full  | 200  |  100|
|run-dim_3-len_125-n_100-epoch_2| 77.346754  | full  | 125 |  100|
|run-dim_30-len_125-n_100-epoch_2-genus|  146.589197 | full  | 125 |  100|
|run-dim_30-len_125-n_200-epoch_2-genus| 199.199229| full  | 125 |  200|
|run-dim_30-len_100-n_100-epoch_5-k_10 (SlimEdwards)|141.658442 | full  | 100 |  100|
|run-dim_50-len_100-n_100-epoch_20-k_10 (SlimEdwards)|192.791254 | full  | 100 |  100|
|run-dim_60-len_125-n_100-epoch_2-k_6-genus|105.597463| full  | 125 |  100|
|run-dim_50-len_100-n_100-epoch_10-k_6_13 (SlimEdwards)|86.545458| full  | 100 |  100|
|run-dim_100-len_125-n_100-epoch_2-k_6-genus |118.997620| full  | 125 |  100|
|run-dim_60-len_125-n_100-epoch_2-k_6-genus-rep_3|128.445961| full  | 125 |  100|
|run-dim_60-len_125-n_100-epoch_2-k_7-none_Slim|128.408734| full  | 125 |  100|
|run-dim_60-len_125-n_300-epoch_2-k_7-none_Slim|294.426457| full  | 125 |  300|
|run-dim_60-len_125-n_600-epoch_2-k_7-none_Slim|512.230759| full  | 125 |  600|

## workflow.py (asymmetric sampling)

|Run id|Runtime [s]  |Run type|--length|--n_vir|--n_host|
|---|---|---|---|---|---|
|run-dim_60-len_125-n_100_600-epoch_2-none_newEdwards| 261.744882| full  | 125  | 100| 600|
|run-dim_60-len_125-n_300_600-epoch_2-none_newEdwards| 261.368303| full  | 125  | 300| 600|
|run-dim_60-len_125-n_300_600-epoch_2-genus| 751.706315| full  | 125  | 300| 600|
|run-dim_60-len_125-n_300_600-epoch_2-none_Slim_newEdwards (new vector saving)| 77.261955| full  | 125  | 300| 600|
|run-dim_60-len_125-n_600_600-epoch_2-k_7-none_Slim|165.377044| full  | 125 |  600| 600|



## faiss_search.py

|Run id|Runtime [s]  |--k_nearest
|---|---|---|
|run-dim_20-len_125-n_100| 53.817716| 20  |
|run-dim_30-len_125-n_100| 74.774101| 20  |
|run-dim_12-len_125-n_100| 43.884163| 20  |
|run-dim_30-len_125-n_100-epoch_2| 89.850484| 20  |
|run-dim_30-len_125-n_200-epoch_2| 199.905919| 20  |
|run-dim_30-len_125-n_100-epoch_3| 90.767702| 20  |
|run-dim_60-len_125-n_100-epoch_2| 146.955118| 20  |
|run-dim_60-len_125-n_100-epoch_2|149.792244 | 60  |
|run-dim_30-len_125-n_100-epoch_2|95.437133| 60  |
|run-dim_60-len_200-n_100-epoch_2|135.712530| 60  |
|run-dim_3-len_125-n_100-epoch_2|39.649473| 60  |
|run-dim_30-len_125-n_100-epoch_2-genus| 84.122965| 60  |
|run-dim_30-len_125-n_200-epoch_2-genus|170.882380 | 60  |
|run-dim_30-len_100-n_100-epoch_5-k_10|78.811480|60|
|run-dim_50-len_100-n_100-epoch_20-k_10|116.156188|60|
|run-dim_60-len_125-n_100-epoch_2-k_6|129.885140|60|
|run-dim_50-len_100-n_100-epoch_10-k_6_13|118.114666|60|
|run-dim_100-len_125-n_100-epoch_2-k_6-genus|200.748213|60|
|run-dim_60-len_125-n_100-epoch_2-k_6|131.922403|30|
|run-dim_60-len_125-n_100-epoch_2-k_6-genus-rep_3|137.502430|60|
|run-dim_60-len_125-n_100-epoch_2-k_7-none_Slim|126.849045|60|
|run-dim_60-len_125-n_300-epoch_2-k_7-none_Slim|445.913480|60|
|run-dim_60-len_125-n_100-epoch_2-k_7-none_Slim (2nd)|128.676293|60|
|run-dim_60-len_125-n_600-epoch_2-k_7-none_Slim|927.781056|60|
|run-dim_60-len_125-n_100_600-epoch_2-none_newEdwards|396.819071|60|
|run-dim_60-len_125-n_300_600-epoch_2-none_newEdwards|432.643582|60|
|run-dim_60-len_125-n_300_600-epoch_2-genus|1565.948420|60|
|run-dim_60-len_125-n_300_600-epoch_2-none_Slim_newEdwards|423.511798|60|
|run-dim_60-len_125-n_600_600-epoch_2-k_7-none_Slim|999.171753|60|
|run-dim_60-len_125-n_600_600-epoch_2-k_7-none_Slim (2nd, sqrt(dist))|990.417626|60|
|run-dim_60-len_125-n_600_600-epoch_2-k_7-none_Slim (3rd, sqrt(dist), new scoring(with filter))|970.415532|60|
|run-dim_60-len_125-n_600_600-epoch_2-k_7-none_Slim (4rd, newest-correct scoring(with filter))|626.470819|60|




## evaluate.py
|Run id|species | genus | family | order | class | phylum | superkingdom
|---|---|---|---|---|---|---|---|
|run-dim_20-len_125-n_100| 4.878048780487805| 12.317073170731707  | 16.70731707317073| 24.024390243902438| 39.268292682926834| 50.36585365853659| 100|
|run-dim_30-len_125-n_100| 5.487804878048781| 13.536585365853659  | 18.536585365853657| 25.121951219512194| 36.46341463414634| 47.19512195121951| 100|
|run-dim_12-len_125-n_100-k_13|1.5853658536585367|4.634146341463414| 7.804878048780488| 11.829268292682926| 27.682926829268297| 41.09756097560975| 100|
|run-dim_30-len_125-n_100-epoch_2|9.75609756097561 | 19.390243902439025 |26.951219512195124 |36.46341463414634 |48.41463414634146 |56.46341463414634|100
|run-dim_30-len_125-n_200-epoch_2|9.146341463414634 |20.48780487804878 |26.707317073170735 |35.12195121951219 |44.87804878048781 |54.146341463414636 |100
|run-dim_30-len_125-n_100-epoch_3| 9.268292682926829|18.902439024390244 |25.365853658536587 |34.02439024390244 |47.5609756097561 |56.82926829268292 |100
|run-dim_60-len_125-n_100-epoch_2|8.536585365853659|15.365853658536585|22.195121951219512|32.5609756097561|44.02439024390244|53.04878048780488|100
|run-dim_60-len_125-n_100-epoch_2 (60 k_nearest)|9.146341463414634 |18.78048780487805 |25.853658536585368|34.8780487804878|46.82926829268293 | 56.58536585365853 | 100
|run-dim_30-len_125-n_100-epoch_2 (60 k_nearest)|11.707317073170733 |20.609756097560975|29.024390243902438|38.41463414634146|47.92682926829268|56.70731707317073|100
|run-dim_60-len_200-n_100-epoch_2 (60 k_nearest)|9.146341463414634|17.4390243902439|21.829268292682926|30.365853658536583|42.3170731707317|51.707317073170735|100
|run-dim_30-len_125-n_100-epoch_2 (60 k_nearest)|4.024390243902439|10.24390243902439|13.170731707317074|21.097560975609756|32.5609756097561|42.68292682926829|100
|run-dim_30-len_125-n_100-epoch_2 (60 k_nearest, genus)|13.90243902439024 | 24.51219512195122 | 31.951219512195124|41.21951219512195|51.951219512195124|59.14634146341463|100
|run-dim_30-len_125-n_200-epoch_2 (60 k_nearest, genus)|13.170731707317074|24.024390243902438|31.341463414634145|40.853658536585364|52.4390243902439|60.48780487804878|100
|run-dim_30-len_100-n_100-epoch_5-k_10(60 k_nearest, SlimEdwards)|12.439024390243903|25.731707317073173|34.512195121951216| 43.78048780487805|56.70731707317073|62.926829268292686|100
|run-dim_50-len_100-n_100-epoch_20-k_10(60 k_nearest, SlimEdwards)|14.02439024390244|25.365853658536587|33.048780487804876|43.41463414634146|55.365853658536594|63.78048780487805|100
|run-dim_60-len_125-n_100-epoch_2-k_6(60 k_nearest, genus)|13.292682926829269| 29.39024390243902| 38.170731707317074|44.87804878048781|57.80487804878048|66.70731707317074|100
|run-dim_50-len_100-n_100-epoch_10-k_6_13(60 k_nearest, SlimEdwards)|12.317073170731707|24.634146341463413|34.390243902439025|42.4390243902439|55.1219512195122|62.5609756097561|100
|run-dim_100-len_125-n_100-epoch_2-k_6 (60 k_nearest, genus)|8.658536585365853|16.82926829268293|18.78048780487805|24.268292682926827|29.146341463414632|32.80487804878049|100
|run-dim_60-len_125-n_100-epoch_2-k_6(30 k_nearest, genus)|12.195121951219512|26.34146341463415|34.512195121951216|41.829268292682926|54.390243902439025|63.78048780487805| 100
|run-dim_60-len_125-n_100-epoch_2-k_6-genus-rep_3(60 k_nearest, genus[reps_3])|11.951219512195122|25.48780487804878|35.24390243902439|45.24390243902439|56.95121951219512|63.65853658536585| 100
|run-dim_60-len_125-n_100-epoch_2-k_7-none_Slim (60 k_nearest, SlimEdwards)|16.21951219512195|29.878048780487802|38.78048780487805|47.682926826829|61.46341463414634|69.26829268292683|100
|run-dim_60-len_125-n_300-epoch_2-k_7-none_Slim (60 k_nearest, SlimEdwards)|19.634146341463417|33.41463414634146|43.292682926829265|49.51219512195122|58.78048780487804|65.0|100
|run-dim_60-len_125-n_100-epoch_2-k_7-none_Slim (60 k_nearest, SlimEdwards, 2nd time; gives same result as the 1st run)|16.21951219512195| 29.878048780487802|38.78048780487805|47.68292682926829|61.46341463414634|69.26829268292683|100
|run-dim_60-len_125-n_600-epoch_2-k_7-none_Slim (60 k_nearest, SlimEdwards)|21.463414634146343|36.21951219512195|44.51219512195122|52.80487804878049|61.09756097560975|68.78048780487805|100
|run-dim_60-len_125-n_100_600-epoch_2-none_newEdwards (60 k_nearest, k=7)|5.731707317073171,|23.902439024390244,|34.268292682926834,|44.75609756097561,|56.95121951219512,|64.02439024390245,|100
|run-dim_60-len_125-n_300_600-epoch_2-none_newEdwards (60 k_nearest, k=7)|8.048780487804878,|29.39024390243902,|40.0,|48.170731707317074,|58.17073170731707,|64.14634146341463,|100
|run-dim_60-len_125-n_300_600-epoch_2-genus (60 k_nearest, k=7)|11.585365853658537,| 29.024390243902438,|38.90243902439025,|46.58536585365854,|57.56097560975609|65.1219512195122|100
|run-dim_60-len_125-n_300_600-epoch_2-none_Slim_newEdwards (60 k_nearest, k=7)|8.414634146341465,| 28.41463414634146,|39.390243902439025,|49.26829268292683,|58.90243902439024,|66.09756097560975,|100
|run-dim_60-len_125-n_600_600-epoch_2-k_7-none_Slim (60 k_nearest)|21.341463414634145,|35.48780487804878,|44.26829268292683,|51.951219512195124,|60.97560975609756,|67.31707317073172,|100
|run-dim_60-len_125-n_600_600-epoch_2-k_7-none_Slim (60 k_nearest, sqrt(dist))|21.341463414634145,|35.48780487804878,|44.26829268292683,|51.951219512195124,|60.97560975609756,|67.31707317073172,|100
|run-dim_60-len_125-n_600_600-epoch_2-k_7-none_Slim (60 k_nearest, sqrt(dist) with filtering)| 24.268292682926827|34.8780487804878|41.46341463414634|48.90243902439024|58.41463414634146| 66.09756097560975|100




