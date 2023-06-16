#!/usr/bin/env python

import os.path

# Get config file
configfile: "config.yaml"

# Dependencies
os.environ['PATH'] += ':' + os.path.abspath("./scripts")
os.environ['PATH'] += ':' + os.path.abspath("./bin")

# Opts
threads = config["threads"]
out_dir = config["outdir"]


# Include rule files
include: "rules/arima_hic_mapping.smk"
include: "rules/yahs.smk"

# Main
rule all:
    input:
        # Hi-C Assembly
        expand(os.path.join(out_dir, "assembly", "hic_contact_map", "{type}_scaffolds_final.hic"), type = types)
