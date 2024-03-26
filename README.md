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

