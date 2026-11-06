"""
Run the analysis from raw Visium files to tables, figures and the run manifest.

    python scripts/run_pipeline.py                  # all stages
    python scripts/run_pipeline.py --stage markers  # one stage (earlier outputs must exist)
"""

from __future__ import annotations

import argparse

from prostate_markers import io, pipeline, repro


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--stage", choices=[*pipeline.STAGES, "all"], default="all")
    parser.add_argument("--skip-report", action="store_true",
                        help="with --stage all: do not execute the notebooks")
    args = parser.parse_args()

    config = io.load_config()
    repro.seed_everything(config["seed"])
    stages = list(pipeline.STAGES) if args.stage == "all" else [args.stage]
    if args.skip_report and "report" in stages:
        stages.remove("report")
    for name in stages:
        print(f"== {name}")
        pipeline.STAGES[name](config)
    print(f"manifest: {repro.write_manifest(config, stages)}")


if __name__ == "__main__":
    main()
