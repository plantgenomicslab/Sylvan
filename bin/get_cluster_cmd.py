#!/usr/bin/env python3
"""Extract cluster_cmd from a Sylvan config YAML."""
import sys
import yaml

cfg = yaml.safe_load(open(sys.argv[1]))
cmd = cfg.get("__default__", {}).get("cluster_cmd", "")
if not cmd:
    print("ERROR: cluster_cmd not found in __default__ section of " + sys.argv[1], file=sys.stderr)
    sys.exit(1)
print(cmd)
