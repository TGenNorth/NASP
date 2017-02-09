````
Hardware Overview:

  Model Name:   MacBook Pro
  Model Identifier: MacBookPro9,2
  Processor Name:   Intel Core i5
  Processor Speed:  2.5 GHz
  Number of Processors: 1
  Total Number of Cores:    2
  L2 Cache (per Core):  256 KB
  L3 Cache: 3 MB
  Memory:   8 GB
  Boot ROM Version: MBP91.00D3.B09
  SMC Version (system): 2.2f44
  Serial Number (system):   X
  Hardware UUID:    X
  Sudden Motion Sensor:
  State:    Enabled
````

1 reference
254 Samples ~4-5MB each
27373 contigs

````
import os
import glob
from concurrent.futures import ProcessPoolExecutor

def indexContigs(filepath):
    with open(filepath) as handle:
        for line in handle:
            if line.startswith('>'):
                pass

with ProcessPoolExecutor() as executor:
    for x in executor.map(printContigs, glob.glob("./benchmarks/Ecoli/external/*.frankenfasta")):
        pass
````

`time python3 indexFasta.py`

real    0m4.195s
user    0m11.078s
sys 0m1.099s


````
if __name__ == "__main__":
    import timeit
    print(timeit.repeat("main()", setup="from __main__ import main", repeat=5, number=1))
````

[3.5490191690041684, 3.5807349779934157, 3.489380874001654, 3.4999573140012217, 3.5834933549995185]

Without ProcessPoolExecutor:
[8.092141925997566, 6.849999589001527, 6.9470494570050505, 7.495032171995263, 8.729224843998963]
