# fastDNA + faiss
## Training - Step 1: fastDNA model learning by random genome sampling from a phylum
- Command (dimensions: 10; k-mer size: 10; training read length: 100; epoch: 1; rest is default):
```
./fastdna supervised -input ../edwards2016/host/random_phylum-training_fastDNA.fasta -labels ../edwards2016/host/random_phylum-training_labels.txt -output edwards_random_model -dim 10 -length 100 -minn 10 -maxn 10 -epoch 1
```
- Head of an example result (from generated `edwards_random_model.vec`)
```
524800 10
AAAAATTTTT -0.21718 -0.14016 -0.088678 0.28453 -0.076305 0.15032 -0.013353 0.04682 0.16902 0.10214
AAAATATTTT -0.069196 -0.015037 0.24694 0.28242 -0.15586 0.11387 0.024382 0.1064 0.049187 -0.066609
AAAACGTTTT -0.08012 0.03275 0.30177 0.35559 -0.13613 0.14376 0.069087 0.022798 0.10009 0.08977
AAAAGCTTTT -0.073477 -0.021923 0.088156 0.21304 -0.20992 0.11418 0.094305 -0.05679 0.076841 -0.087441
AAAAAATTTT 0.035327 0.022978 0.21736 0.36057 -0.20622 0.19239 -0.0081592 0.0045573 0.014566 0.025942
AAAACCTTTT -0.079999 0.14628 0.24882 0.55499 -0.22793 0.14801 0.13735 0.077292 0.2075 -0.01836
AAAAACTTTT -0.06527 0.050299 0.064983 0.050619 -0.03619 0.0065102 -0.070489 0.058526 0.082506 0.058876
AAAACATTTT -0.098118 0.0020106 -0.040881 0.011984 0.0008345 -0.037997 0.096516 0.092891 -0.010287 -0.039297
AAAAAGTTTT -0.088758 0.060384 0.096763 -0.089787 -0.052039 -0.07906 0.0069594 0.0042814 -0.034004 -0.036472
```
## Training - Step 2: Sampling sequences from testing dataset to be converted into embeded vectors
- Execution through `random-input-sampling.py` (with path to all host fasta files)
- Head of a generated input file (from `random_samples-training_fastDNA.fasta`)
```
>NC_000117.1 Chlamydia trachomatis D/UW-3/CX chromosome, complete genome sample_0 gen_pos:(256244:256444)
AAAGCTAGAATTGACAAAATTAACTTTTTAACAATTTTTCTTAACATCCCCTATTAGAGA
GCTCCTTTCTTGCAAAAAAAATCTCTCCTTTCTTACATTTAGTCTCAAAAGATATACCCC
CTCTCAGAAAAACTGAGGTTTTATGCGCTCACCCCTGAAACTTTTATTCTCTCCTGAAAA
AAGCAATATGATGCTTGGCT
>NC_000117.1 Chlamydia trachomatis D/UW-3/CX chromosome, complete genome sample_1 gen_pos:(463552:463752)
GCCCACTTAACAGCTGCTTTAGAAGCATCCACATGATAGACTGTTGCTCCTTGTTGAGCA
CAGAAAATAGAAGAAGCTCCAGTATAAGCAAAGAGATTCAAAACTCTACATGAAGGCTTG
GCTACAGAAGGTTGTAAATCTTTCCAGAACCCAGCATGCTCGGGAAAAATACCTACATGA
CCAAAAGGAGTCAGCTTCAA
```
- Head of a generated sample map (from `sample_map.json`)
> for debugging purposes in json; it will be binary file later
```
{
    "24759": "1163617",
    "24760": "1163617",
    "24761": "1163617",
    "24762": "1163617",
    "24763": "1163617",
    "24764": "1163617",
    "24765": "1163617",
    "24766": "1163617",
    "24767": "1163617",
```
## Training - Step 3: Converting samples to embeded vectors
- Execution through command
```commandline
./fastdna print-word-vectors edwards_random_model.bin < ../edwards2016/host/random_samples-training_fastDNA.fasta > edwards_vectors.txt
```
## Training - Step 4: Building index of host sequences for faiss
- Execution through `index_building.py`
- Generates `host_index.index`

## Example Usage - Step 1: Sampling sequences from query virus genome
- Execution through `random-input-sampling.py` (with path to single virus fasta file)

## Example Usage - Step 2: Converting samples to embeded vectors
- Execution through command
```commandline
./fastdna print-word-vectors edwards_random_model.bin < ../edwards2016/virus/random_samples_virus.fasta > virus_vectors.txt
```
## Example Usage - Step 3: Similarity search with faiss and mapping results
- Execution through `faiss_search.py`
- Example output:
```
Nearest sequences: [[ 3628 17198 15422  7129  2733 18960  7024 16364 18322 16436]
 [  466 23447 10805 17525  3625 11387  4702 17369 19358  5405]
 [24002 19986 26176 13620 11949 13078  1176  8449  6151  7853]
 [13437  9854 16159 14061 14911 14899  3125  2919  4106  4363]
 [14911  8059 18800  7949 13437  8239  3341  4182 21629  5089]
 [12423  1597 12766  2938 20496 10188 12129 10738 20386  7870]
 [19174  5884  1270  5465  3740 22173  7792 20303 22502 21337]
 [23618 22214  4177 10093 13568 14019  9164 14631  9187 15299]
 [16216 19499 15739 14214 16921  4382  2375 22244 26485  7020]
 [13710  2062 14668 20496 14935 20404 14638 14309  9553 18125]]
Classification: [['NC_008319' 'NC_017293' 'NC_016781' 'NC_010741' 'NC_007624' 'NC_017622'
  'NC_010676' 'NC_017052' 'NC_017512' 'NC_017068']
 ['NC_002655' 'NC_020888' 'NC_013941' 'NC_017317' 'NC_008319' 'NC_014219'
  'NC_009052' 'NC_017308' 'NC_017656' 'NC_009654']
 ['NC_021172' 'NC_018020' 'NC_022546' 'NC_015591' 'NC_014550' 'NC_015174'
  'NC_004129' 'NC_011985' 'NC_010103' 'NC_011353']
 ['NC_015458' 'NC_013223' 'NC_016947' 'NC_015714' 'NC_016589' 'NC_016511'
  'NC_007802' 'NC_007677' 'NC_008554' 'NC_008711']
 ['NC_016589' 'NC_011745' 'NC_017573' 'NC_011601' 'NC_015458' 'NC_011831'
  'NC_007973' 'NC_008618' 'NC_019396' 'NC_009481']
 ['NC_014623' 'NC_005835' 'NC_015052' 'NC_007626' 'NC_018290' 'NC_013526'
  'NC_014640' 'NC_013946' 'NC_018150' 'NC_011386']
 ['NC_017646' 'NC_009832' 'NC_004369' 'NC_009648' 'NC_008347' 'NC_019897'
  'NC_011369' 'NC_018143' 'NC_019973' 'NC_018700']
 ['NC_020895' 'NC_019952' 'NC_008570' 'NC_013446' 'NC_015564' 'NC_015663'
  'NC_012857' 'NC_016010' 'NC_012815' 'NC_016745']
 ['NC_016934' 'NC_017790' 'NC_016812' 'NC_015733' 'NC_017217' 'NC_008705'
  'NC_007086' 'NC_019936' 'NC_022738' 'NC_010676']
 ['NC_015596' 'NC_006576' 'NC_016024' 'NC_018290' 'NC_016516' 'NC_018268'
  'NC_016010' 'NC_015757' 'NC_013118' 'NC_017454']]
```
## Visualisation of vector clouds
- Execution through `vector_visualisation.py`
- Very big plot is generated, takes long time to read in the browser and works slow
- Example screenshots:
![Edwards vectors - whole](../master/docs/edwards_vectors.png)
![Edwards vectors - up view](../master/docs/edwards_vectors-up.png)
