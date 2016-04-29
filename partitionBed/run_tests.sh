#/bin/bash

./partitionBed.py test_data/a.bed test_data/b.bed > test_data/observed.bed

echo "OK IF YOU SEE NOTHING:"
diff test_data/observed.bed test_data/expected-1.bed