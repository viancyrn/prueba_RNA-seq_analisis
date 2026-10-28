# Command documentation

This directory contains explanations of the Bash scripts used in the RNA-seq analysis workflow.

Each document describes the purpose of a script, its inputs and outputs, the commands and parameters used, and the scientific or technical reasons behind the workflow choices.

## Common Bash conventions

### Bash interpreter

All workflow scripts begin with:

`#!/usr/bin/env bash`

This shebang asks the operating system to locate Bash through the current environment rather than relying on a fixed Bash installation path.

### Strict error handling

The scripts use:

`set -euo pipefail`

- `-e`: exit when an unhandled command fails.
- `-u`: report an error when an undefined variable is used.
- `-o pipefail`: report failures occurring within pipelines, not only failures from the final command.

These options help detect problems early, although they do not replace validation of scientific results.

### Project-relative paths

The scripts use:

`cd "$(dirname "$0")/../.."`

This changes the working directory to the project root, assuming the script is stored in `code/bash/`. Consequently, file paths can be written relative to the project root.

## Workflow documentation

Individual Markdown files document the analysis steps, from read acquisition and quality control to alignment, quantification, and MultiQC reporting.

The alignment-based and Salmon-based workflows share the same final trimmed FASTQ files but use different quantification approaches.
