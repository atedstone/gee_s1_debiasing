# README

Developed by Andrew Tedstone, 2021-2024.

## Dependencies

`atedstone/gee_tools.git`

## Usage

Requires a job configuration TOML file - see the example included.

Main use case involves doing debiasing first, then working with the outputs.

One way of doing this is in a Jupyter Notebook:

```python
%run core_debiasing_workflow.py <confg_file.toml>
```

And then work with the output variables, e.g. `refs_all`, `s1_debiased`.

For a practical use case see scripts in `paper_rlim_retention_repo.git`.