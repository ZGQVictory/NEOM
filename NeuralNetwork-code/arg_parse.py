import argparse


def get_args():
    """Parse input arguments"""
    parser = argparse.ArgumentParser(
        description="This is zgq demo about hallucinating peptide on the MHC",
        usage="input peptide length (-l) and MHC type (-mhc); output target sequence of peptide",
    )

    # Basic arguments
    parser.add_argument("-l", type=int, default=9, help="the length of the peptide")
    parser.add_argument("-mhc", type=str, default="HLA-A*0201", help="the MHC type")
    parser.add_argument(
        "-start_pep",
        type=str,
        default=None,
        help="DEFAULT = None, which means we generate peptide using our start frequency. If you want a specific start point, please enter it here.",
    )
    parser.add_argument(
        "-start_fre",
        type=list,
        default=None,
        help="DEFAULT = None, which means we generate peptide using our default frequency (IEDB+BLOSUM). If you want a specific start point, please enter 1 in the argument and put the frequency matrix in the code.",
    )
    parser.add_argument(
        "-frequency_weight",
        default="1_1",
        help="our start frequency_map consists of two parts, 1:BLOSUM62, 2:IEDB_0201, default weights are [1,1]",
    )
    parser.add_argument(
        "-rest_mut_pos",
        default=None,
        type=str,
        help="Select some specific positions of the neo-antigen, especially the hydrophobic residues in the middle, to maintain the hydrophobicity. Example: 3_5_7 (default: %(default)s).",
    )
    parser.add_argument(
        "-fre_mut_pos",
        default="1_2_3_4_5_6_7_8_9",
        type=str,
        help="Select the position that you want to mutate using our frequency matrix (default: %(default)s).",
    )
    parser.add_argument(
        "-o",
        default="./example",
        type=str,
        help="your output direction",
    )
    parser.add_argument(
        "-verbose",
        action="store_true",
        help="if you want to see more details in the Details file, use -verbose.",
    )
    # MCMC arguments
    parser.add_argument(
        "-tolerance",
        default=None,
        action="store",
        type=float,
        help="the tolerance on the loss sliding window for terminating the MCMC trajectory early (default: %(default)s).",
    )

    parser.add_argument(
        "-mutation_rate",
        default="3_1",
        action="store",
        help="number of mutations at each MCMC step (start-finish, stepped linear decay). Should probably be scaled with protomer length (default: %(default)s)",
    )

    parser.add_argument(
        "-T_init",
        default=0.025,
        action="store",
        type=float,
        help="starting temperature for simulated annealing. Temperature is decayed exponentially (default: %(default)s).",
    )

    parser.add_argument(
        "-half_life",
        default=1000,
        action="store",
        type=float,
        help="half-life for the temperature decay during simulated annealing (default: %(default)s).",
    )

    parser.add_argument(
        "-metropolis_hasting",
        action="store_true",
        help="if you want to use the Metropolis Hasting algorithm, remeber to check the reference probability in tcr-specific--peptide.txt, use -metropolis_hasting.",
    )

    parser.add_argument(
        "-steps",
        default=10,
        action="store",
        type=int,
        help="number for steps for the MCMC trajectory (default: %(default)s).",
    )

    # Our policy
    parser.add_argument(
        "-add_booster",
        action="store_true",
        help="if you want to use the booster that generate the peptdies similar to the TCR file, remeber to check the file of tcr-specific--peptide.txt, using -add_booster.",
    )

    parser.add_argument(
        "-booster_num",
        default=3,
        type=int,
        help="the position number that should be kept as same as one sequence of the TCR file (default: %(default)s).",
    )

    parser.add_argument(
        "-frequence_change_rate",
        default=0.3,
        type=float,
        help="change the frequency matrix using loss function (default: %(default)s)",
    )

    parser.add_argument(
        "-mutant_times",
        default=10,
        type=int,
        help="In one step, we mutate 10 times and then decide acceptance, DEFAULT = %(default)s.",
    )

    # Loss arguments
    parser.add_argument(
        "-TCR_loss_pos",
        default="3_4_5_6_7_8",
        type=str,
        help="exert the TCR_loss on the specific positions, DEFAULT = %(default)s.",
    )

    parser.add_argument(
        "-cdr_sequence",
        default=None,
        type=str,
        help="put the cdr loop sequence in the new loss, DEFAULT=None. For example, extracted from 5NMG, like 'YSDRGSQSFFWMFIYSNGDKAVRTNSGYALNCSPKQGHDTVSIFQYYEEEERQRGDTVSYEQY'.",
    )

    parser.add_argument(
        "-weight_cdr_dis",
        default=0.1,
        type=float,
        help="the weight of cdr discrete loss, DEFAULT = %(default)s.",
    )
    parser.add_argument(
        "-weight_cdr",
        default=2,
        type=float,
        help="the weight of cdr loss, which you need a sequence, DEFAULT = %(default)s.",
    )
    parser.add_argument(
        "-weight_iedb",
        default=0.05,
        type=float,
        help="the weight of iedb discrete loss, DEFAULT = %(default)s.",
    )
    ## Other parts
    parser.add_argument(
        "-nomemory",
        action="store_true",
        help="if you want to try this process with a changing frequency but no memory, use -nomemory.",
    )
    parser.add_argument(
        "-freqnotchange",
        action="store_true",
        help="if you want to try this process with a fixed frequency, use -freqnotchange.",
    )

    args = parser.parse_args()
    return args
