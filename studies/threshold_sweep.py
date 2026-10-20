import numpy as np
from src.boys import boys_array
from tests.boys_reference import boys_hyp1f1
from pathlib import Path

with open(Path(__file__).resolve().parent / "threshold_sweep.txt", "w") as file:
    file.write(f"{'s':>20} {'l':>20} {'worst x':>20} {'max err':>20}\n")

    small_x_grill = np.linspace(0.1, 5, 100)
    larger_x_grill = np.linspace(6, 100, 100)
    x_grill = np.logspace(-15, 5, 1000)
    for s in small_x_grill:
        for l in larger_x_grill:
            worst_x = 0
            max_e = 0
            for x in x_grill:
                result = boys_array(4, x, s, l)
                oracle = [boys_hyp1f1(n,x) for n in range(5)]
                err = np.max(abs(result - oracle)/oracle)
                if err > max_e:
                    worst_x = x
                    max_e = err
            file.write(f"{s:20.4f} {l:20} {worst_x:20} {max_e:20.3e}\n")


