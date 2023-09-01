import random

def get_tmax_and_n(filename):
    tmax, n = None, None
    lines = open(filename, "r").readlines()
    for line in lines:
        line = line.split()
        if tmax is None:
            tmax = int(line[2])
        else:
            tmax = max(tmax, int(line[2]))
        if n is None:
            n = max(int(line[0]), int(line[1]))
        else:
            n = max(int(line[0]), int(line[1]), n)
    return tmax, n

if __name__ == "__main__":
    tmax, n = get_tmax_and_n("graph.txt")
    contents = ""
    num_of_queries = input("Enter the number of queries to be generated: ")
    try:
        num_of_queries = int(num_of_queries)
    except:
        print("Invalid input! Program terminated.")
        exit()
    length = int(tmax * float(input("Enter the length of the query windows (0 < x < 1): ")))
    for i in range(num_of_queries):
        u, v = random.randint(1, n), random.randint(1, n)
        ts = random.randint(0, tmax - length)
        contents += str(u) + " " + str(v) + " " + str(ts) + " " + str(ts + length) + " " + "\n"
    open("query.txt", "w").write(contents)