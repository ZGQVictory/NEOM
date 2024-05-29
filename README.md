# EvaST

EvaST is a novel neoantigen maturation framework encompassing five distinct modules: “policy”, “structure”, “evaluation”, “selection” and “filter”. 

It is Flexible, Interpretable, Precise and Fast.

Hope you enjoy it.

---
# Requirements

EvaST requires [PANDORA](https://github.com/X-lab-3D/PANDORA/), python and some python libraries. The following installations are required before EvaST installation:

* Python >=3.7
* conda
* pip3

1. Generate the environment

      ```bash
      conda env create -n EvaST python=3.9
      conda activate EvaST
      ```

2. Install the [PANDORA](https://github.com/X-lab-3D/PANDORA/)

      In this paper, we use PANDORA-v1.0.0. You should follow the [ReadME.md](https://github.com/X-lab-3D/PANDORA/blob/master/README.md) step by step.

      When constructing the PANDORA database, it's better to set its directory in the './data/pmhc/', which means that there will be a 'pandora_Database.pkl' in the './data/pmhc/'. Or, change the PANDORA database directory in the './code/main.py' to fit your directory. (In the main.py, I highlighted three places that you can change.)

      After PANDORA installation, please use python to make a validation:

      ```python
      ## import requested modules
      from PANDORA import Target
      from PANDORA import Pandora
      from PANDORA import Database
      
      ## A. Load local Database
      db = Database.load()
      
      ## B. Create Target object
      target = Target(id = 'myTestCase',
          allele_type = 'HLA-A*0201',
          peptide = 'LLFGYPVYV',
          anchors = [2,9])
      
      ## C. Perform modelling
      case = Pandora.Pandora(target, db)
      case.model()

      ```

4. Install EvaST

    Clone the repository:
   ```bash
   git clone https://github.com/ZGQVictory/EvaST.git
   cd EvaST
   ```

# Command options

```python
python main.py -h
```

You will see a bunch of options. Below are the details of each argument that can be used:

- `-h, --help`: Show this help message and exit.
- `-l L`: The length of the peptide. **DEFAULT = 9**.
- `-mhc MHC`: The MHC type. **DEFAULT = HLA-A*0201**.
- `-start_pep START_PEP`: Indicates we generate peptide using our start frequency. If you want a specific start point, please enter it here. **DEFAULT = None**.
- `-start_fre START_FRE`: We generate peptide using our default frequency (IEDB+BLOSUM). To specify another frequency, enter 1 in the argument and put the frequency matrix in `process.py`. **DEFAULT = None**.
- `-frequency_weight FREQUENCY_WEIGHT`: Our default start frequency map consists of two parts, 1:BLOSUM62, 2:IEDB_0201, with default weights of [1,1]. **DEFAULT = 1_1**.
- `-rest_mut_pos REST_MUT_POS`: Select specific positions of the neo-antigen, especially the hydrophobic residues in the middle, to maintain hydrophobicity. Example: 3_5_7. **DEFAULT = None**.
- `-fre_mut_pos FRE_MUT_POS`: Select the position you want to mutate using our frequency matrix. **DEFAULT = 1_2_3_4_5_6_7_8_9**.
- `-o O`: Your output direction. **DEFAULT = ./**.
- `--total_structures`: To save all files, especially the PANDORA results. 
- `--verbose`: For more details in the Details file.
- `--extract_pdbfile`: To extract pdbfile to a LOSS directory.
- `--tolerance TOLERANCE`: The tolerance on the loss sliding window for terminating the MCMC trajectory early. **DEFAULT = None**.
- `-mutation_rate MUTATION_RATE`: Number of mutations at each MCMC step (start-finish, stepped linear decay), likely scaled with protomer length. **DEFAULT = 3-1**.
- `--T_init T_INIT`: Starting temperature for simulated annealing, decayed exponentially. **DEFAULT = 25**.
- `--half_life HALF_LIFE`: Half-life for temperature decay during simulated annealing. **DEFAULT = 1000**.
- `--metropolis_hasting`: Use the Metropolis Hasting algorithm, checking the reference probability in tcr-specific--peptide.txt.
- `--steps STEPS`: Number of steps for the MCMC trajectory. **DEFAULT = 10**.
- `--add_booster`: To use the booster generating peptides similar to the TCR file, checking tcr-specific--peptide.txt.
- `--booster_num BOOSTER_NUM`: The position number to keep as one sequence of the TCR file. **DEFAULT = 3**.
- `--frequence_change_rate FREQUENCE_CHANGE_RATE`: Rate of change for the frequency matrix. **DEFAULT = 0.3**.
- `--mutant_times MUTANT_TIMES`: Mutation times in one step. **DEFAULT = 10**.
- `--step_model_num STEP_MODEL_NUM`: Number of PANDORA models in one step. **DEFAULT = 5**.
- `--num_cores NUM_CORES`: Number of CPUs to use, see data using `lscpu`. **DEFAULT = 16**.
- `--TCR_loss_pos TCR_LOSS_POS`: Exert the TCR_loss on specific positions. **DEFAULT = 1_2_3_4_5_6_7_8_9**.
- `--cdr_sequence CDR_SEQUENCE`: Specify the cdr loop sequence in the new loss, e.g., extracted from 5NMG. **DEFAULT = None**.
- `--weight_cdr_dis WEIGHT_CDR_DIS`: The weight of cdr discrete loss. **DEFAULT = 100**.
- `--weight_cdr WEIGHT_CDR`: The weight of cdr loss, needing a sequence. **DEFAULT = 2**.
- `--weight_iedb WEIGHT_IEDB`: The weight of iedb discrete loss. **DEFAULT = 50**.
- `--nomemory`: To try this process with a changing frequency but no memory.
- `--freqnotchange`: To try this process with a fixed frequency.

---
# Test

Test example, 
1. run 1000 steps maturation from 'VMNILLQYV' with 12 cpu cores (--steps 1000, -start_pep VMNILLQYV, --num_cores 12)
2. mutate 10 times in one step (--mutant_times 10)
3. set starting temperature to be 100 (to change the acceptance rate, --T_init 100)
4. exert the TCR loss on the peptide position 3~8 (--TCR_loss_pos 3_4_5_6_7_8)
5. extract the pdb files that PANDORA generated (--extract_pdbfile)
6. add the booster for the mutations (--add_booster)
7. set the output directory (-o './output_data')

then run the following command to generate an Expanded peptide pool:

```bash
python ./code/main.py --steps 1000  -start_pep VMNILLQYV --mutant_times 10  --T_init 100  --TCR_loss_pos 3_4_5_6_7_8  --extract_pdbfile --add_booster -o './output_data'
```

# Filter Options

1. [NetMHC-4.0](https://services.healthtech.dtu.dk/services/NetMHC-4.0/)
2. [IEDB-Consensus](https://nextgen-tools.iedb.org/pipeline)
3. [Episcan predictor](https://www.episcan-predictor.com)
4. Free Energy Perturbation
5. Stability MD validation




