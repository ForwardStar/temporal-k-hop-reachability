import os
import shutil
import subprocess

if __name__ == "__main__":
    if not os.path.exists("datasets"):
        print("No datasets detected!")
        exit()
    if not os.path.exists("main"):
        print("No executable file detected!")
        exit()
    num_dataset = len(os.listdir("datasets"))
    for i in range(num_dataset + 1):
        os.system("python3 graph-gen.py %s" % (i))
        if os.path.exists("./results/%s" % (i)):
            shutil.rmtree("./results/%s" % (i))
        os.makedirs("./results/%s" % (i))
        for k in [2, 4, 6]:
            os.makedirs("./results/%s/k=%s" % (i, k))
            os.system("python3 query-gen.py 1000 %s -1" % (k))
            try:
                output = subprocess.check_output(["./main", "graph.txt", "query.txt", "./results/%s/k=%s/output_online1.txt" % (i, k), "Online1", str(k), "Config"], timeout=259200)
                with open("./results/%s/k=%s/online1.txt" % (i, k), "wb") as f:
                    f.write(output)
            except:
                print("Timed out!")
            try:
                output = subprocess.check_output(["./main", "graph.txt", "query.txt", "./results/%s/k=%s/output_online2.txt" % (i, k), "Online2", str(k), "Config"], timeout=259200)
                with open("./results/%s/k=%s/online2.txt" % (i, k), "wb") as f:
                    f.write(output)
            except:
                print("Timed out!")
            try:
                output = subprocess.check_output(["./main", "graph.txt", "query.txt", "./results/%s/k=%s/output_MP.txt" % (i, k), "MP", str(k), "Config"], timeout=259200)
                with open("./results/%s/k=%s/MP.txt" % (i, k), "wb") as f:
                    f.write(output)
            except:
                print("Timed out!")
            try:
                output = subprocess.check_output(["./main", "graph.txt", "query.txt", "./results/%s/k=%s/output_T2H.txt" % (i, k), "T2H", str(k), "Config"], timeout=259200)
                with open("./results/%s/k=%s/T2H.txt" % (i, k), "wb") as f:
                    f.write(output)
            except:
                print("Timed out!")
            os.system("cmp ./results/%s/k=%s/output_online2.txt ./results/%s/k=%s/output_online1.txt" % (i, k, i, k))
            os.system("cmp ./results/%s/k=%s/output_MP.txt ./results/%s/k=%s/output_online1.txt" % (i, k, i, k))
            os.system("cmp ./results/%s/k=%s/output_T2H.txt ./results/%s/k=%s/output_online1.txt" % (i, k, i, k))