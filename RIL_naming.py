import os
import re

matches = set()
for filename in os.listdir("."):
    sep = filename.split("_")
    if len(sep) > 1:
        matches.add(sep[0])

if matches:
    with open("RIL_list.txt", "w") as file1:
        file1.writelines(' '.join(matches))
