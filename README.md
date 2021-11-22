# fastDNA + faiss
## Step 0: Configuration file set up - config.cfg
This needs to be set up correctly, otherwise none of the modules will not work.
> [GENERAL]

- fastdna_dir = path to fastDNA directory; (should be a version from this repository, because of its enchanced features, then path is following: `./fastDNA/`)
- active_model = path to a model file (.bin) to be used 

> [HOST]

- host_genomes = path to all host genomes fasta files (.fna)

> [VIRUS]

- virus_genomes = path to all virus genomes fasta files (.fna)

## Step 1: Model creation - make-model.py module
- Samples all given host genomes according to chosen taxonomy level or takes whole genome set to create model
- Created model (.bin) is for fastDNA vector embedding
- Optionally saves model in readable form (.vec)
- Usage:
```commandline
python make-model.py -h

usage: make-model.py [-h] -i INPUT_DIR -o OUTPUT -f {phylum,class,order,family,genus,species,none} [-r REPS] -d DIM -l
                     LENGTH --minn MINN --maxn MAXN -e EPOCH -t THREAD [--rm] [--saveVec]

fastDNA model creation

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_DIR, --input_dir INPUT_DIR
                        Directory with host genomes.
  -o OUTPUT, --output OUTPUT
                        Path to result FASTA file, labels file and model file.
  -f {phylum,class,order,family,genus,species,none}, --filter {phylum,class,order,family,genus,species,none}
                        Taxonomy level to which genomes should be filtered. Choosing 'none' implies no taxonomy
                        filtering.
  -r REPS, --reps REPS  Maximum number of representatives from the filtered group. Default value is 1.
  -d DIM, --dim DIM     Dimensionality of vectors
  -l LENGTH, --length LENGTH
                        Length of sequences
  --minn MINN           Minimum k-mer size
  --maxn MAXN           Maximum k-mer size (max k=15, otherwise fastDNA fails)
  -e EPOCH, --epoch EPOCH
                        Number of epochs (each added epoch increases runtime significantly)
  -t THREAD, --thread THREAD
                        Number of threads to use
  --rm                  Remove potentially redundant files after model creation. Default is 'false'.
  --saveVec             Enables saving of a readable model file (.vec). Enabling this may significantly increase
                        execution time. Default is 'false'.
```
- Recommended example usage (space-saving)
```commandline
python make-model.py -i /home/hyperscroll/edwards2016/host/fasta/ -o /home/hyperscroll/edwards2016/models/epoch/ --filter phylum -d 30 --length 125 --minn 10 --maxn 10 --epoch 2 -t 16 --rm
```
## Step 2: Files preparation - workflow.py module
- Performs random sampling from host and virus genomes according to chosen parameters
- Creates host and virus embedded vectors
- Creates necessary files for _**faiss**_ 
- Whole process deploys into own container:
```
.
└─ user_defined_container_name/
    ├─ host/
    │   ├─ index/
    │   │   └─ host_index.index
    │   ├─ maps/
    │   │   └─ sample_map.json
    │   ├─ samples/
    │   │   └─ host_samples.fasta
    │   └─ vectors/
    │       └─ host_vectors.txt
    ├─ rank/
    │   └─ user_defined_rank_name.json
    └─ virus/
        ├─ samples/
        │   └─ {virus_nbci-id}_vector.txt
        └─ vectors/
            └─ {virus_accession}_samples.fasta
```
- Usage
```commandline
python workflow.py -h

usage: workflow.py [-h] -w WD [--host] [--virus] [--full] --length LENGTH --n_vir N_VIR --n_host N_HOST -t THREAD

fastDNA+faiss virus-host interaction analysis

optional arguments:
  -h, --help            show this help message and exit
  -w WD, --wd WD        Working directory, where all files will be deployed
  --host                Host mode: every available host genome is randomly sampled according to a given criteria, then
                        a cloud of host vectors is generated and compiled into a faiss index
  --virus               Virus mode: every available virus genome is randomly sampled according to a given criteria,
                        then a cloud of virus vectors is generated which is compared with host cloud and results are
                        generated in a form of a rank of virus-host pairs.
  --full                Full mode: Combines host and virus mode in one go.
  --length LENGTH       Length of the samples
  --n_vir N_VIR         Number of samples to take from a virus genome
  --n_host N_HOST       Number of samples to take from a host genome
  -t THREAD, --thread THREAD
                        Number of threads to use
```
## Step 3: Search and rank interactions - faiss_search.py module
- Runs _**faiss**_ to find each virus k-nearest host sequences 
- Rank is saved to a .json file
- Usage:
```commandline
python faiss_search.py -h

usage: faiss_search.py [-h] -i INPUT_DIR -o OUTPUT -k K_NEAREST -f FAISS_INDEX -m MAP

faiss search deployer

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_DIR, --input_dir INPUT_DIR
                        Directory with virus samples of vectors.
  -o OUTPUT, --output OUTPUT
                        Path to result file with rankings.
  -k K_NEAREST, --k_nearest K_NEAREST
                        How many (k) nearest samples should be found
  -f FAISS_INDEX, --faiss_index FAISS_INDEX
                        Path to faiss index
  -m MAP, --map MAP     Path to sample map
```
- Output example:
```json
{
    "NC_000866": [
        [
            "NC_017353",
            21.320540770000004
        ],
        [
            "NC_017481",
            19.493801010000006
        ],
        [
            "NC_017333",
            19.45372254
        ],
        [
            "NC_007350",
            19.394375320000005
        ],
        [
            "NC_009879",
            19.02006546
        ],
        [
            "NC_017347",
            18.12450329
        ],
        [
            "NC_008710",
            17.945988370000002
        ],
    ...
```
## Step 4: Prediction evaluation - evaluate.py module
- Returns percent of correct predictions for each taxonomy level from Edwards2016 dataset.
- Example output:
```json
{'species': 11.951219512195122, 'genus': 25.48780487804878, 'family': 35.24390243902439, 'order': 45.24390243902439, 'class': 56.95121951219512, 'phylum': 6
3.65853658536585, 'superkingdom': 100.0}
```
- More results in `tests.md` file.
## Visualisation of vector clouds
- Execution through `vector_visualisation.py` (maual code editing to filter and manipulate the data frame)
- Data frame exported to .csv
- Csv from `vector_visualisation.py` can visualised in [EVCVTool](https://mega.nz/folder/yr52VKSL#55mQ4PSO3C1eUqMAJecemA), a custom 3D viewer for virus-host embedded vectors.
  - Correct csv form:
  ```csv
  ,x,y,z,ncbi_id,sample_name,full_name,vir_x (optional),vir_y (optional),vir_z (optional),vir_ncbi_id (optional),vir_sample_name (optional),vir_full_name (optional)
  ```
  - Csv files needs to be put into following directory `EVCVTool_DX11-x.x\Embedded Vector Cloud Visualisation Tool_DX11_Data\Resources\DataClouds\`
  - Example csv files are located in the `./visualisation/`
- Example screenshots:
![EVCV - host cluster](../master/docs/evcv_1.png)
![EVCV - host-virus cluster](../master/docs/evcv_2.png)
![EVCV - host cluster with organism names](../master/docs/evcv_3.png)
