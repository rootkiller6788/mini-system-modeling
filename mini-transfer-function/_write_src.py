
import os

def w(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)
    lc = text.count(chr(10))
    print(f"  wrote {path}: {lc} lines, {len(text)} bytes")
    return lc

# Read and write tf_analysis.c
c3 = open("_tf_analysis.c.tmp", "r", encoding="utf-8").read()
w("src/tf_analysis.c", c3)
os.remove("_tf_analysis.c.tmp")
print("Done writing tf_analysis.c")
