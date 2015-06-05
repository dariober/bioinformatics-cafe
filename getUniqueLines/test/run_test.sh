#!/bin/bash

echo -e "\n--- No duplicates. INPUT == OUTPUT"
cat dups-1.txt | ../getUniqueLines.py

echo -e "\n--- Only duplicates. EMPTY OUTPUT"
sort dups-1.txt dups-1.txt | ../getUniqueLines.py

echo -e "\n--- 1 line only in input -> 1 line output"
head -n 1 dups-1.txt | ../getUniqueLines.py

echo -e "\n--- 1st line duplicate: Skip aaa line only"
head -n 1 dups-1.txt | cat - dups-1.txt | ../getUniqueLines.py

echo -e "\n--- Only last line duplicate:"
tail -1 dups-1.txt | cat dups-1.txt - | ../getUniqueLines.py

echo -e "\n--- Mix of dups and unique"
cat dups-2.txt | ../getUniqueLines.py

echo -e "\n--- Catch exception: 3 lines equal:"
sort dups-1.txt dups-1.txt dups-1.txt | ../getUniqueLines.py