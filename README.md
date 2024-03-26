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

3. Install EvaST

    Clone the repository:
   ```bash
   git clone https://github.com/ZGQVictory/EvaST.git
   cd EvaST
   ```

# Command options

```python
python main.py -h
```

You will see a bunch of options:

- `-h`, `--help`: Show this help message and exit.
- `-l L`: The length of the peptide.
- `-mhc MHC`: The MHC type.
- `-start_pep START_PEP`: DEFAULT = None, which means we generate peptide using our start frequency. If you want a specific start point, please enter it here.
- `-start_fre START_FRE`: DEFAULT = None, which means we generate peptide using our default frequency (IEDB+BLOSUM). If you want a specific start point, please enter 1 in the argument and put the frequency matrix in the code.
- `-frequency_weight FREQUENCY_WEIGHT`: Our start frequency_map consists of two parts, 1:BLOSUM62, 2:IEDB_0201, default weights are [1,1].
- `-rest_mut_pos REST_MUT_POS`: Select some specific positions of the neo-antigen, especially the hydrophobic residues in the middle, to maintain the hydrophobicity. Example: 3_5_7 (default: None).
- `-fre_mut_pos FRE_MUT_POS`: Select the position that you want to mutate using our frequency matrix (default: 1_2_3_4_5_6_7_8_9).
- `-o O`: Your output direction.
- `--total_structures`: If you want to save all the files especially the PANDORA results, use --total_structures.
- `--verbose`: If you want to see more details in the Details file, use --verbose.
- `--extract_pdbfile`: If you want to extract pdbfile to a LOSS directory, use --extract_pdbfile.
- `--tolerance TOLERANCE`: The tolerance on the loss sliding window for terminating the MCMC trajectory early (default: None).
- `-mutation_rate MUTATION_RATE`: Number of mutations at each MCMC step (start-finish, stepped linear decay). Should probably be scaled with protomer length (default: 3-1).
- `--T_init T_INIT`: Starting temperature for simulated annealing. Temperature is decayed exponentially (default: 25).
- `--half_life HALF_LIFE`: Half-life for the temperature decay during simulated annealing (default: 1000).
- `--metropolis_hasting`: If you want to use the Metropolis Hasting algorithm, remember to check the reference probability in tcr-specific--peptide.txt, use --metropolis_hasting.
- `--steps STEPS`: Number for steps for the MCMC trajectory (default: 10).
- `--add_booster`: If you want to use the booster that generate the peptides similar to the TCR file, remember to check the file of tcr-specific--peptide.txt, using --add_booster.
- `--booster_num BOOSTER_NUM`: The position number that should be kept as same as one sequence of the TCR file (default: 3).
- `--frequence_change_rate FREQUENCE_CHANGE_RATE`: Change the frequency matrix using loss function (default: 0.3).
- `--mutant_times MUTANT_TIMES`: In one step, we mutate 10 times and then decide acceptance, DEFAULT = 10.
- `--step_model_num STEP_MODEL_NUM`: Number about the model in a pandora generation loop, DEFAULT = 5.
- `--num_cores NUM_CORES`: Number of CPUs, you can look at the data using command -> lscpu, DEFAULT = 16.
- `--TCR_loss_pos TCR_LOSS_POS`: Exert the TCR_loss on the specific positions, DEFAULT = 1_2_3_4_5_6_7_8_9.
- `--cdr_sequence CDR_SEQUENCE`: Put the cdr loop sequence in the new loss, DEFAULT=None. For example, extracted from 5NMG, like 'YSDRGSQSFFWMFIYSNGDKAVRTNSGYALNCSPKQGHDTVSIFQYYEEEERQRGDTVSYEQY'.
- `--weight_cdr_dis WEIGHT_CDR_DIS`: The weight of cdr discrete loss, DEFAULT = 100.
- `--weight_cdr WEIGHT_CDR`: The weight of cdr loss, which you need a sequence, DEFAULT = 2.
- `--weight_iedb WEIGHT_IEDB`: The weight of iedb discrete loss, DEFAULT = 50.
- `--nomemory`: If you want to try this process with a changing frequency but no memory, use --nomemory.
- `--freqnotchange`: If you want to try this process with a fixed frequency, use --freqnotchange.






