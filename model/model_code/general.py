import re
# use 32_32_64.log code and 1024_1024_1024.log test
def remove_symbols(line):
    if "GPU activities:" in line:
        line = line.replace("GPU activities:", " ")
    line = re.sub(r'\s+', ' ', line)
    line = re.sub(r'(\d+(\.\d+)?)ms', lambda m: str(float(m.group(1)) * 1), line)
    line = re.sub(r'(\d+(\.\d+)?)us', lambda m: str(float(m.group(1)) / 1000), line)
    line = re.sub(r'(\d+(\.\d+)?)s', lambda m: str(float(m.group(1)) * 1000), line)
    return line

def remove_mgard(line):
    line = re.sub(r'\s+', ' ', line)
    line = re.sub(r'(\d+(\.\d+)?) ms', lambda m: str(float(m.group(1)) * 1), line)
    line = re.sub(r'(\d+(\.\d+)?) us', lambda m: str(float(m.group(1)) / 1000), line)
    line = re.sub(r'(\d+(\.\d+)?) s', lambda m: str(float(m.group(1)) * 1000), line)
    # print(line)
    return line
