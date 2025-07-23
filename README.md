# NEOM

NEOM is a neoantigen maturation framework with five key modules: **policy · structure · evaluation · selection · filter**.  
We now provide **two runnable implementations**:

- **Part 1 – Structure-driven version (root `main.py`)**  
  Original pipeline that relies on structural modeling (PANDORA required).

- **Part 2 – Neural-network version (`NeuralNetwork-code/main.py`)**  
  End-to-end learning pipeline **without PANDORA**.  
  **Recommended** for quick use and less setup time.

---

## Table of Contents

- [NEOM](#neom)
  - [Table of Contents](#table-of-contents)
  - [Repository Layout](#repository-layout)
  - [Environment Setup](#environment-setup)
    - [Create Conda env (example)](#create-conda-env-example)
    - [PANDORA (Only for Part 1)](#pandora-only-for-part-1)
  - [Part 1: Structure-driven NEOM (root `main.py`)](#part-1-structure-driven-neom-root-mainpy)
  - [Part 2: Neural-network NEOM (`NeuralNetwork-code/main.py`)](#part-2-neural-network-neom-neuralnetwork-codemainpy)
  - [Post-processing Logs](#post-processing-logs)
  - [Filtering \& Downstream Tools](#filtering--downstream-tools)
  - [Citation \& License](#citation--license)

---

## Repository Layout

```

NEOM/
├── analyse_script/                 # (optional) analysis helpers
├── data/
├── example/
│   ├── random_Details_logfile.txt
│   ├── random_Details_logfile_step_loss.log
│   ├── random_Details_logfile_step_loss_accept.log
│   └── random_Frequencymap.txt
├── NeuralNetwork-code/
│   ├── database/                   # NN-specific resources
│   ├── datanet/
│   ├── example/
│   ├── Training_scripts/
│   ├── environment.yml
│   ├── main.py                     # NN version entry point
│   └── ...
├── addition_func.py
├── arg_parse.py
├── Loss.py
├── main.py                         # structure version entry point
├── new_loss.py
├── Policy.py
├── process.py
├── Structure_generation.py
└── README.md

````

---

## Environment Setup

> If you only plan to run the **Neural-network version**, you can skip PANDORA completely.

### Create Conda env (example)

```bash
conda env create -n NEOM python=3.9
conda activate NEOM
````

(Or adapt from your own `environment.yml` if provided.)

### PANDORA (Only for Part 1)

Follow PANDORA’s README: [https://github.com/X-lab-3D/PANDORA](https://github.com/X-lab-3D/PANDORA) (paper used v1.0.0).
Place `pandora_Database.pkl` under `./data/pmhc/` **or** adjust the path in `main.py` (`###Change 3###`).

Quick validation after installing PANDORA:

```python
from PANDORA import Target, Pandora, Database
db = Database.load()
target = Target(id='myTestCase', allele_type='HLA-A*0201',
                peptide='LLFGYPVYV', anchors=[2,9])
case = Pandora.Pandora(target, db)
case.model()
```

---

## Part 1: Structure-driven NEOM (root `main.py`)

Run:

```bash
python main.py [options]
```

Common arguments (subset):

| Arg                  | Meaning                           | Default                   |
| -------------------- | --------------------------------- | ------------------------- |
| `-l`                 | peptide length                    | 9                         |
| `-mhc`               | MHC type                          | HLA-A\*0201               |
| `-start_pep`         | start peptide sequence            | None                      |
| `-start_fre`         | use custom start frequency        | None                      |
| `-frequency_weight`  | weights for BLOSUM62 & IEDB\_0201 | 1\_1                      |
| `-rest_mut_pos`      | positions to keep hydrophobic     | None                      |
| `-fre_mut_pos`       | positions allowed to mutate       | 1\_2\_3\_4\_5\_6\_7\_8\_9 |
| `-o`                 | output directory                  | ./                        |
| `--total_structures` | save all (esp. PANDORA outputs)   | False                     |
| `--extract_pdbfile`  | dump PDBs to LOSS dir             | False                     |
| `--steps`            | MCMC steps                        | 10                        |
| `-mutation_rate`     | per-step mutations (e.g. `3-1`)   | 3-1                       |
| `--T_init`           | initial temperature               | 25                        |
| `--half_life`        | temp half-life                    | 1000                      |
| `--mutant_times`     | mutations per step                | 10                        |
| `--num_cores`        | CPU cores                         | 16                        |
| `--TCR_loss_pos`     | positions for TCR loss            | 1\_2\_3\_4\_5\_6\_7\_8\_9 |
| `--add_booster`      | booster mode                      | False                     |

**Example** (1000 steps from `VMNILLQYV`, 12 CPUs, save PDB, add booster, etc.):

```bash
python main.py --steps 1000 -start_pep VMNILLQYV \
  --mutant_times 10 --T_init 100 --TCR_loss_pos 3_4_5_6_7_8 \
  --extract_pdbfile --add_booster -o ./output_data --num_cores 12
```

---

## Part 2: Neural-network NEOM (`NeuralNetwork-code/main.py`)

> **No PANDORA needed.** This version is faster to set up.

1. Enter the folder:

   ```bash
   cd NeuralNetwork-code
   ```

2. (Optional) build env from its `environment.yml`:

   ```bash
   conda env create -f environment.yml -n NEOM_NN
   conda activate NEOM_NN
   ```

3. Run:

   ```bash
   python main.py [options]
   ```

Parameters are analogous but specific to the NN pipeline (see `arg_parse.py` inside this folder). Adjust dataset / database paths as needed under `database/`, `datanet/`, etc.

---

## Post-processing Logs

After **either** version finishes:

```bash

cp analyse_script/ extract_step_loss_sequence.py .
python extract_step_loss_sequence.py ./example random_Details_logfile.txt
```

This generates:

* `random_Details_logfile_step_loss.log`
* `random_Details_logfile_step_loss_accept.log`  ← **all accepted peptides**

Both appear under `./example/` by default (change path if needed).

---

## Filtering & Downstream Tools

To further prioritize binders/immunogenic peptides, you may use:

1. NetMHC-4.0
2. IEDB Consensus
3. Episcan Predictor

(Links can be found in the original README; keep them in your workflow as needed.)

---

## Citation & License

If you use NEOM in academic work, please cite our paper.
https://www.biorxiv.org/content/10.1101/2024.08.14.607669v1

---

**Tip:** Prefer the **Neural-network version** for a quicker start, then switch to the structure-driven version if you need explicit structural modeling and interpretability.

Enjoy!


