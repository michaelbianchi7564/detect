def validate_config(cfg):
    if "strandedness" not in cfg or cfg["strandedness"] is None:
        raise ValueError("strandedness must be set to 0, 1, or 2")

    if cfg["strandedness"] not in [0, 1, 2]:
        raise ValueError("Strandedness must be 0, 1, or 2")


# PARAMTER HELPER FUNCTION
def dict_to_cli(cfg, exclude="extra",):                            # Define a dictionary for parameters
    cli = []                                      # Create an empty list of parameters
    for p, v in cfg.items():                     # Iterate al parameters
        if v is None or v is False or v == "":    # For empty parameters
            continue
        elif v is True:                           # For True/False parameters
            cli.append(f"--{p}")                  # Append the parameter
        else:
            cli.append(f"--{p} {v}")              # For parameters with values
    return " ".join(cli)

# EXTRA PARAMETER HELPER FUNCTION
def add_extra(base, cfg):
    return base + (" " + cfg["extra"] if cfg.get("extra") else "")    # return the resul of dict_to_cli function and add extra parameters

# HISAT2 PARAMETER HELPER FUNCTION
def build_hisat2_cli(cfg):
    args = []

    if cfg.get("preset"):                                # Specify the alignment algorithm
        args.append(f"--{cfg['preset'].lstrip('-')}")    # Append the preset parameter 
    if cfg.get("dta"):                                   # Specify the presence/absence of dta parameter
        args.append("--dta")                             # Append dta if it is true
    input_format = cfg.get("input_format", "fastq")
    if input_format == "fastq":                          # Append -q if fastq input is specified
        args.append("-q")
    elif input_format  == "fasta":                       # Append -f if fasta input is specified
        args.append("-f")

    if cfg.get("extra"):                                 # Append extra parameters if present
        args.append(cfg["extra"])

    return " ".join(args)

# FEATURECOUNTS HELPER FUNCTION
def featurecounts_cli(cfg):
    args = []

    if cfg.get("paired_end"):                       # Add -p parameter if samples are paired-end
        args.append("-p")

    args +=["-s", str(cfg.get("strandness", 0))]    # set strandness, if the parameter is empty, add 0

    if cfg.get("feature_type"):
        args.append(f"-t {cfg['feature_type']}")

    if cfg.get("attribute"):
        args.append(f"-g {cfg['attribute']}")

    return " ".join(args)
