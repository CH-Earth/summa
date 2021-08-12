# Case study: Reynolds Mountain East
This experiment compares the impact of two different albedo decay parametrizations on simulations of snow depth for a specific mountain side. Before starting, ensure you have a compiled SUMMA executable. Then generate the simulations by running SUMMA twice, with different settings:

```
summa.exe -m fileManager_reynoldsConstantDecayRate.txt
summa.exe -m fileManager_reynoldsVariableDecayRate.txt
```

Evaluation data and a sample (Python) script to evaluate these simulations are provided in sub-folder `evaluation` and `output` respectively.