import re

def extract_info(filename):
# Split file names into components
    parts = filename.split("/")
    filename1 = parts[-1]
    file1_without_extension = filename1.split(".txt")[0]
    components = file1_without_extension.split("_")
    # Extract information from components
    if len(components) == 8:
        name, d1, d2, d3, err, d1_s, d2_s, d3_s = components[0], int(components[1]), int(components[2]),int(components[3]), \
        float(components[4]),int(components[5]), int(components[6]),int(components[7])

    else:
        print("Invalid file name:", filename1)
    print(name, d1, d2, d3, err, d1_s, d2_s, d3_s)
    return name, d1, d2, d3, err, d1_s, d2_s, d3_s

#extract_info("GCLDLWP_26_1800_3600_0.1_13_900_1800.txt")
# extract_info("GCLDLWP_26_1800_3600_0.1_6400_stat.txt")
